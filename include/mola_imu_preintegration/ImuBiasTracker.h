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
#include <mrpt/poses/CPose3D.h>

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
 *    over a window of `standstill_window` seconds, the gyroscope and accelerometer
 *    readings stay close to their means (vibration and noise allowed), the mean
 *    angular velocity stays close to the current bias, and the mean specific force
 *    has the magnitude of gravity. The window mean of the gyroscope is then a
 *    measurement of its bias. For the accelerometer, only the component along
 *    gravity is observable without knowing the attitude, so it updates just that
 *    one: `|mean(a)| - g = u^T b_acc`, with `u` the measured gravity direction.
 *    Unlike ImuInitialCalibrator, it keeps doing this at every stop, not only at
 *    startup.
 *  - **External measurements** (add_measurement()), e.g. the bias estimated by
 *    MapGravityEstimator while moving, with their covariance.
 *
 *  Between measurements, each covariance grows with the random-walk parameters, so
 *  an old estimate is reported as less certain.
 *
 *  ### Limitation of the standstill detector
 *
 *  From the IMU alone, a perfectly smooth slow rotation (slower than
 *  `standstill_max_gyro_mean_offset`, without vibration) cannot be told apart from
 *  a bias, and would be taken as one. On vehicles that can do that (e.g. a slow
 *  turn in place), feed motion hints from another source (wheel or LiDAR
 *  odometry) with add_motion_hint(): any hint saying "moving" inside a window
 *  vetoes it, and `require_static_hint` makes a "static" hint mandatory.
 *
 *  ### Frames and units
 *
 *  Biases are in the IMU sensor frame, as raw readings and
 *  ImuIntegrationParams::bias_gyro / bias_acc: rad/s and m/s^2. Feed raw
 *  readings, not ImuTransformer outputs. Time stamps are seconds (TimeStamp).
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
        /// Length of each (non-overlapping) standstill window [s].
        double standstill_window = 1.0;
        /// Minimum readings in a window.
        std::size_t standstill_min_samples = 20;
        /// A gap between consecutive readings longer than this restarts the window [s].
        double max_sample_gap = 0.25;
        /// Maximum RMS deviation of the angular velocity from its window mean
        /// [rad/s]: noise and vibration of a vehicle standing with its motors on.
        double standstill_max_gyro_rms = 0.02;
        /// Maximum difference between the window mean angular velocity and the
        /// current bias estimate [rad/s]: rejects steady rotations. Must be larger
        /// than the expected bias error before the first standstill.
        double standstill_max_gyro_mean_offset = 0.01;
        /// Maximum RMS deviation of the specific force from its window mean [m/s^2].
        double standstill_max_acc_rms = 0.3;
        /// Maximum difference between |mean specific force| and `gravity` [m/s^2].
        double standstill_max_acc_norm_error = 0.3;
        /// Gravity magnitude [m/s^2].
        double gravity = 9.81;
        /// If true, a window is only taken as standstill if a motion hint saying
        /// "static" falls inside it (and none says "moving").
        bool require_static_hint = false;
        /** @} */

        /** @name Bias model
         * @{ */
        /// Bias random walk [rad/s/sqrt(s)] and [m/s^2/sqrt(s)].
        double gyro_bias_random_walk = 2e-5;
        double acc_bias_random_walk  = 2e-4;
        /// Sigma of the zero-mean prior, before any measurement [rad/s] and [m/s^2].
        double gyro_bias_prior_sigma = 0.02;
        double acc_bias_prior_sigma  = 0.3;
        /// Floors of the standstill measurement sigmas [rad/s] and [m/s^2], for
        /// what the sample statistics cannot see (temperature, vibration
        /// rectification, a not quite still vehicle).
        double gyro_measurement_min_sigma = 2e-4;
        double acc_measurement_min_sigma  = 0.02;
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

        /// Time of the covariances above: the newest measurement, or the query time.
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
     *  Readings must be in time order per IMU. Readings without angular velocity
     *  are ignored.
     *  \return true if it closed a standstill window that updated the estimate.
     */
    bool add(const mrpt::obs::CObservationIMU& imu);

    /** Tells whether the vehicle (all its IMUs) is static at time `t`, from another
     *  source. Hints are kept for a short time (a few windows).
     */
    void add_motion_hint(TimeStamp t, bool isStatic);

    /** Fuses an external bias measurement of one IMU, e.g. from MapGravityEstimator.
     *  Either part may be omitted.
     */
    void add_measurement(
        const std::string& sensorLabel, TimeStamp t,
        const std::optional<std::pair<AngularVelocity, mrpt::math::CMatrixDouble33>>&    gyro,
        const std::optional<std::pair<LinearAcceleration, mrpt::math::CMatrixDouble33>>& acc);

    /** The current estimate of one IMU, or nothing if that label was never seen.
     *  If `t` is given (and newer than the last update), the covariances are
     *  propagated to that time.
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

    struct PerImu
    {
        Estimate           est;
        std::deque<Sample> window;
    };

    std::map<std::string, PerImu> imus_;
    std::map<TimeStamp, bool>     hints_;

    PerImu& imu_state(const std::string& sensorLabel, TimeStamp t);
    /// Propagates the covariances of an estimate to time t (random walk).
    void predict(Estimate& e, TimeStamp t) const;
    /// Checks the current window; updates the estimate if it is a standstill.
    bool               process_window(PerImu& s);
    [[nodiscard]] bool hints_allow(TimeStamp from, TimeStamp to) const;
};

}  // namespace mola::imu

/** Feature macro: ImuBiasTracker is available. */
#define MOLA_IMU_PREINTEGRATION_HAS_BIAS_TRACKER 1
