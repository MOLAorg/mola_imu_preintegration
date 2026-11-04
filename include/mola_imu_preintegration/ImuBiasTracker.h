/*               _
 _ __ ___   ___ | | __ _
| '_ ` _ \ / _ \| |/ _` | Modular Optimization framework for
| | | | | | (_) | | (_| | Localization and mApping (MOLA)
|_| |_| |_|\___/|_|\__,_| https://github.com/MOLAorg/mola

 Copyright (C) 2018-2026 Jose Luis Blanco, University of Almeria,
                         and individual contributors.
 SPDX-License-Identifier: GPL-3.0
 See LICENSE for full license information.
 Closed-source licenses available upon request, for this odometry package
 alone or in combination with the complete SLAM system.
*/
#pragma once

#include <mola_imu_preintegration/types.h>
#include <mrpt/containers/yaml.h>
#include <mrpt/math/CMatrixFixed.h>
#include <mrpt/math/TPoint3D.h>
#include <mrpt/obs/CObservationIMU.h>

#include <cstddef>
#include <deque>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace mola::imu
{
/** Tracks the gyroscope and accelerometer biases of one or more IMUs over time,
 *  one estimate per sensor label, as slowly varying random walks.
 *
 *  Two kinds of measurements update each estimate (a Kalman filter per IMU):
 *
 *  - **Standstill windows**, detected by the tracker itself from the raw readings:
 *    over a sliding window of `standstill_window` seconds, the gyroscope and
 *    accelerometer readings stay close to their means (vibration and noise allowed),
 *    the mean angular velocity is consistent with the current bias estimate and its
 *    uncertainty, and the mean specific force has the magnitude of gravity. The
 *    window mean of the gyroscope is then a measurement of its bias. For the
 *    accelerometer, only the component along gravity is observable without knowing
 *    the attitude, so it updates just that one: `|mean(a)| - g = u^T b_acc`, with
 *    `u` the measured gravity direction (stops at different tilts make the rest
 *    observable). That component is relative to `gravity_magnitude`: an error in it, or in the
 *    accelerometer scale, ends up in the bias. Accepted windows do not overlap.
 *    Unlike ImuInitialCalibrator, it keeps doing this at every stop, not only at
 *    startup.
 *  - **External measurements** (add_measurement()), with their covariance, in the
 *    sensor frame. E.g. ImuInitialCalibrator results, or MapGravityEstimator biases
 *    (rotate them from the vehicle to the sensor frame, and fuse only solves over
 *    non-overlapping data, or they are counted more than once).
 *
 *  Between measurements, each covariance grows with the random-walk parameters, so
 *  an old estimate is reported as less certain.
 *
 *  ### Limitations of the standstill detector
 *
 *  From the IMU alone, a smooth motion with a nearly constant angular velocity and
 *  little vibration cannot be told apart from a standstill: a slow turn in place,
 *  a wide curve at constant speed, slow pitching over undulating ground. Its
 *  rotation would be taken as bias. Feed motion hints from another source (wheel
 *  or LiDAR odometry) with add_motion_hint(): any hint saying "moving" during a
 *  window vetoes it, and `require_static_hint` makes a "static" hint mandatory.
 *  Once a hint is received, readings are processed only when the hints reach their
 *  time (or `max_hint_latency` after it), so late hints still count.
 *
 *  The Earth rotation (up to 7.3e-5 rad/s) is part of the measured gyroscope bias.
 *
 *  ### Frames and units
 *
 *  Biases are in the IMU sensor frame, as raw readings and
 *  ImuIntegrationParams::bias_gyro / bias_acc: rad/s and m/s^2. Feed raw
 *  readings, not ImuTransformer outputs. Time stamps are seconds (TimeStamp).
 *  IMUs slower than `standstill_min_samples / standstill_window` Hz never pass
 *  the standstill test.
 *
 *  Thread-safety: none. Callers serialize access externally.
 *
 *  \ingroup mola_imu_preintegration_grp
 */
class ImuBiasTracker
{
   public:
    ImuBiasTracker() = default;

    struct Parameters
    {
        Parameters() = default;

        /** @name Standstill detection
         * @{ */
        /// Length of the sliding standstill window [s].
        double standstill_window = 1.0;
        /// Minimum readings in a window.
        std::size_t standstill_min_samples = 20;
        /// A gap between consecutive readings longer than this restarts the window [s].
        double max_sample_gap = 0.25;
        /// Maximum RMS deviation of the angular velocity from its window mean
        /// [rad/s]: noise and vibration of a vehicle standing with its motors on.
        double standstill_max_gyro_rms = 0.02;
        /// The window mean angular velocity is accepted as a bias measurement if
        /// it differs from the current estimate by less than the largest of this
        /// [rad/s] and `standstill_gyro_offset_sigmas` times the sigma of that
        /// difference (estimate propagated to the window, plus the window's own).
        /// Rejects steady rotations once the estimate has converged, while a large
        /// bias is still accepted at first. An offset within the gate but beyond the
        /// estimate uncertainty inflates the covariance, so the estimate follows it.
        double standstill_max_gyro_mean_offset = 0.01;
        double standstill_gyro_offset_sigmas   = 3.0;
        /// Windows rejected by the gate above but agreeing with each other, for this
        /// long in total [s] (over interruptions: stops, driving, gaps), mean that the
        /// bias has changed (e.g. after a power cycle): its covariance is reset to
        /// the prior. Twice `standstill_window` suffices if a "static" motion hint
        /// touches them. Offsets beyond what the prior allows for a bias
        /// (`standstill_gyro_offset_sigmas` prior sigmas) are always rotations.
        /// A smooth constant-rate rotation within that bound, lasting this long
        /// without motion hints, would be taken as a bias change too.
        double gyro_reacquire_time = 30.0;
        /// Maximum RMS deviation of the specific force from its window mean [m/s^2].
        double standstill_max_acc_rms = 0.3;
        /// Maximum difference between |mean specific force| and `gravity_magnitude` [m/s^2].
        double standstill_max_acc_norm_error = 0.3;
        /// Gravity magnitude [m/s^2].
        double gravity_magnitude = 9.81;
        /** @} */

        /** @name Motion hints
         * @{ */
        /// If true, a window is only taken as standstill if a motion hint saying
        /// "static" covers part of it (and none says "moving").
        bool require_static_hint = false;
        /// A hint holds until the next one, for at most this long [s].
        double motion_hint_hold = 1.0;
        /// Once hints are in use, a reading waits for a hint at or after its time,
        /// at most this long after newer readings arrived [s].
        double max_hint_latency = 0.5;
        /** @} */

        /** @name Bias model
         * @{ */
        /// Bias random walk [rad/s/sqrt(s)] and [m/s^2/sqrt(s)].
        double bias_gyro_random_walk = 2e-5;
        double bias_acc_random_walk  = 2e-4;
        /// Sigma of the zero-mean prior, before any measurement [rad/s] and [m/s^2].
        double bias_gyro_prior_sigma = 0.02;
        double bias_acc_prior_sigma  = 0.3;
        /// Floors of each standstill window measurement sigma [rad/s] and [m/s^2].
        double standstill_gyro_min_sigma = 2e-4;
        double standstill_acc_min_sigma  = 0.02;
        /// Floors of the estimate sigmas, per axis [rad/s] and [m/s^2]: errors that
        /// do not average out over many windows (temperature, vibration
        /// rectification, gravity and scale errors).
        double bias_gyro_min_sigma = 1e-4;
        double bias_acc_min_sigma  = 0.01;
        /** @} */

        void load_from(const mrpt::containers::yaml& cfg);
        void save_to(mrpt::containers::yaml& cfg) const;
    };

    Parameters parameters;

    /// The bias estimate of one IMU.
    struct Estimate
    {
        Estimate() = default;

        AngularVelocity             bias_gyro{0, 0, 0};  //!< [rad/s], sensor frame
        LinearAcceleration          bias_acc{0, 0, 0};  //!< [m/s^2], sensor frame
        mrpt::math::CMatrixDouble33 cov_gyro;  //!< [(rad/s)^2]
        mrpt::math::CMatrixDouble33 cov_acc;  //!< [(m/s^2)^2]

        /// Time of the covariances above.
        TimeStamp stamp = 0;
        /// Time of the newest measurement (0: none yet, the prior).
        TimeStamp last_update = 0;

        std::size_t standstill_updates = 0;
        std::size_t external_updates   = 0;

        [[nodiscard]] bool has_measurements() const
        {
            return standstill_updates + external_updates > 0;
        }
    };

    /** Feeds one raw reading (sensor frame) of the IMU named by its sensorLabel.
     *  Readings must be in time order per IMU: duplicated or slightly older stamps
     *  are dropped, and a jump back in time of more than 1 s (e.g. a dataset
     *  restart) restarts that IMU, keeping its bias estimate (its `stamp` and
     *  `last_update` move back to the new time). Readings without angular
     *  velocity or a valid time stamp, or with non-finite or implausible values,
     *  are ignored (only the accelerometer part, if just that one is bad).
     *  \return how many times the estimate of that IMU was updated. Once motion
     *  hints are in use, most updates happen in add_motion_hint() instead.
     */
    std::size_t add(const mrpt::obs::CObservationIMU& imu);

    /** Tells whether the vehicle (all its IMUs) is static at time `t`, from another
     *  source, in the same clock as the IMU readings.
     *
     *  Hints are point events that hold forward, until the next one (at most
     *  `motion_hint_hold`): to describe a motion increment over [t_a, t_b], give
     *  the same hint at both ends. Hints of all sources are merged ("moving"
     *  wins at equal stamps), and the newest of any source releases the
     *  readings waiting for hints.
     *  \return the labels of the IMUs whose estimate was updated, as readings
     *          waiting for this hint are processed.
     */
    std::vector<std::string> add_motion_hint(TimeStamp t, bool isStatic);

    /** Processes all the readings still waiting for motion hints, e.g. at the end
     *  of a dataset. \return the labels of the IMUs updated. */
    std::vector<std::string> flush();

    /** Fuses an external bias measurement of one IMU (sensor frame), with its
     *  covariance. Either part may be omitted. A measurement older than the
     *  estimate is fused with its covariance grown by the random walk since then.
     */
    void add_measurement(
        const std::string& sensorLabel, TimeStamp t,
        const std::optional<std::pair<AngularVelocity, mrpt::math::CMatrixDouble33>>&    gyro,
        const std::optional<std::pair<LinearAcceleration, mrpt::math::CMatrixDouble33>>& acc);

    /** The current estimate of one IMU, or nothing if that label was never seen.
     *  Covariances are propagated to `t` if given (and newer than the estimate),
     *  else to the newest reading of that IMU.
     */
    [[nodiscard]] std::optional<Estimate> estimate(
        const std::string& sensorLabel, const std::optional<TimeStamp>& t = std::nullopt) const;

    /// Labels of all the IMUs seen so far.
    [[nodiscard]] std::vector<std::string> sensor_labels() const;

    /// Forgets everything (estimates, windows and hints), keeping the parameters.
    void reset();

   private:
    struct Sample
    {
        TimeStamp                            t = 0;
        mrpt::math::TVector3D                w{0, 0, 0};
        std::optional<mrpt::math::TVector3D> a;
    };

    /// A window of readings with running sums, so each new one costs O(1).
    struct Window
    {
        std::deque<Sample>          samples;
        mrpt::math::TVector3D       sum_w{0, 0, 0};
        mrpt::math::TVector3D       sum_ww{0, 0, 0};  //!< per axis squares
        std::size_t                 na = 0;  //!< samples with accelerometer
        mrpt::math::TVector3D       sum_a{0, 0, 0};
        mrpt::math::CMatrixDouble33 sum_aa;  //!< sum of a * a^T

        void push(const Sample& s);
        void pop();
        void clear() { *this = Window(); }
    };

    struct PerImu
    {
        Estimate                 est;
        std::deque<Sample>       input;  //!< readings waiting for motion hints
        Window                   window;
        std::optional<TimeStamp> newest;  //!< newest accepted reading time

        /// Re-acquisition of a changed bias: mean of the first window rejected by
        /// the offset gate, and the time covered by those agreeing with it.
        struct Reacquire
        {
            mrpt::math::TVector3D mean{0, 0, 0};
            double                time    = 0;  //!< [s]
            TimeStamp             last_t1 = 0;
        };
        std::optional<Reacquire> reacquire;
    };

    enum class HintVerdict
    {
        Unknown,
        Static,
        Moving
    };

    std::map<std::string, PerImu> imus_;
    std::map<TimeStamp, bool>     hints_;
    bool                          hints_used_  = false;
    TimeStamp                     newest_hint_ = 0;

    PerImu& imu_state(const std::string& sensorLabel, TimeStamp t);
    /// Propagates the covariances of an estimate to time t (random walk).
    void predict(Estimate& e, TimeStamp t) const;
    /// Applies the posterior sigma floors.
    void apply_floors(Estimate& e) const;
    /// Processes the waiting readings that motion hints allow (all, if `all`).
    /// Returns the number of updates.
    std::size_t drain(PerImu& s, bool all);
    /// Adds one reading to the window and checks it. True if updated.
    bool process_sample(PerImu& s, const Sample& smp);
    /// Checks the current window; updates the estimate if it is a standstill.
    bool                      try_update(PerImu& s);
    [[nodiscard]] HintVerdict hints_verdict(TimeStamp from, TimeStamp to) const;
    void                      prune_hints();
};

}  // namespace mola::imu

/** Feature macro: ImuBiasTracker is available. */
#define MOLA_IMU_PREINTEGRATION_HAS_BIAS_TRACKER 1
