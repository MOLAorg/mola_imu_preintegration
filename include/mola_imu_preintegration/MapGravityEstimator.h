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

/**
 * @file   MapGravityEstimator.h
 * @brief  Estimates the gravity vector in the odometry/map frame, from
 *         preintegrated IMU plus externally-supplied poses and velocities.
 * @author Jose Luis Blanco Claraco
 * @date   Jul 22, 2026
 */
#pragma once

#include <mola_imu_preintegration/ImuPreintegrator.h>
#include <mrpt/containers/yaml.h>
#include <mrpt/core/optional_ref.h>

#include <cstddef>
#include <deque>
#include <optional>
#include <string>

namespace mola::imu
{
/** Estimates the GRAVITY VECTOR IN THE MAP FRAME (plus the IMU biases) from
 *  preintegrated IMU measurements, using the attitudes and velocities of an
 *  external odometry source (e.g. LiDAR odometry) over a sliding window.
 *
 *  ### The idea
 *
 *  The preintegrated velocity delta satisfies, exactly,
 *
 *      R_i * dV_ij  =  v_j - v_i - g * dt_ij
 *
 *  so, rearranged,
 *
 *      g  =  ( v_j - v_i - R_i * dV_ij ) / dt_ij
 *
 *  Every interval is therefore a direct measurement of gravity expressed in
 *  whatever frame `R` and `v` live in, i.e. the odometry/map frame. Body
 *  acceleration is carried by `dV_ij` and cancels out, which is exactly what
 *  averaging raw accelerometer samples cannot do: that only works while the
 *  vehicle is quasi-static, whereas this works under arbitrary motion.
 *
 *  How far the estimated `g` has departed from `(0, 0, -g)` IS the tilt the
 *  map has accumulated, and `correction()` is the rotation that removes it.
 *
 *  ### What is estimated, and what is trusted
 *
 *  Unknowns are just `(g, bias_acc, bias_gyro)`: NINE scalars, regardless of
 *  the window length. There are deliberately no per-keyframe attitude
 *  variables:
 *  - the external odometry's RELATIVE rotations over a short window are
 *    accurate (that is the strength of scan matching); it is the ABSOLUTE
 *    tilt that drifts, and that is precisely what `g` absorbs;
 *  - with no absolute-attitude variables there is no gauge freedom to pin and
 *    no frame-mismatch between the IMU's gravity-aligned world and the
 *    odometry's map frame.
 *
 *  The gyroscope bias stays observable through a second residual comparing the
 *  preintegrated rotation against the external odometry's relative rotation.
 *
 *  ### Observability caveat (do not remove the bias priors)
 *
 *  The accelerometer bias is only separable from a small tilt given ROTATION
 *  EXCITATION: on a perfectly straight, level, constant-attitude run, a bias
 *  and a tilt are indistinguishable. The zero-anchored bias priors are
 *  therefore load-bearing, not decorative. Likewise, the gravity-magnitude
 *  soft constraint is what keeps the direction well conditioned when the
 *  velocity data is weak.
 *
 *  ### Quality is REPORTED, not gated
 *
 *  This class deliberately has no "the estimate is good enough" threshold.
 *  `Result::converged` means only that the solve produced a usable estimate
 *  (see its docs); how much to trust it is answered by the EARNED
 *  `pitch_sigma` / `roll_sigma`, which the consumer should use as the weight
 *  of whatever constraint it builds. An internal quality threshold would
 *  merely hide informative-but-uncertain estimates behind a boolean, and a
 *  threshold tight enough to look reassuring silently disables the estimator
 *  on real MEMS hardware, where honest per-window tilt sigmas of several
 *  degrees are normal under vehicle motion.
 *
 *  Thread-safety: none. Callers serialize access externally.
 *
 *  \ingroup mola_imu_preintegration_grp
 */
class MapGravityEstimator
{
   public:
    MapGravityEstimator() = default;

    struct Parameters
    {
        Parameters() = default;

        /// Number of intervals kept in the sliding window.
        std::size_t window_size = 20;

        /// Nominal local gravity magnitude [m/s^2].
        double gravity_magnitude = 9.81;

        /// Sigma of the soft constraint on |g| [m/s^2]. Keeps the direction
        /// well conditioned; loose enough to absorb a real local-gravity or
        /// scale-factor discrepancy.
        double gravity_magnitude_sigma = 0.05;

        /// An interval whose integrated time spans less than this fraction of
        /// its wall-clock duration is REJECTED: a partially covered interval
        /// yields a delta that does not span the gap, which would silently
        /// corrupt the constraint.
        double min_interval_coverage = 0.95;

        /// Intervals outside this duration range [s] are rejected. Too short
        /// makes `g = (...)/dt` ill-conditioned; too long breaks the
        /// constant-bias and first-order-bias-update assumptions.
        double min_interval_duration = 0.05;
        double max_interval_duration = 30.0;

        /// Assumed sigma of the external odometry's RELATIVE rotation over one
        /// interval [degrees]. Only affects how strongly the gyro bias is
        /// pulled toward agreeing with the odometry.
        double odometry_relative_rotation_sigma_deg = 0.5;

        /// Sigma of the zero-anchored accelerometer bias prior [m/s^2].
        double bias_acc_prior_sigma = 0.15;

        /// Sigma of the zero-anchored gyroscope bias prior [rad/s].
        double bias_gyro_prior_sigma = 0.01;

        /// Huber threshold applied to the whitened gravity residuals, so a bad
        /// odometry velocity or a jolt cannot drag the estimate.
        double huber_threshold = 2.0;

        std::size_t max_iterations = 10;

        /// Gauss-Newton stops when the increment norm falls below this.
        double tolerance = 1e-10;

        /// Minimum number of intervals in the window before the result is
        /// reported as usable.
        std::size_t min_intervals_for_convergence = 4;

        void load_from(const mrpt::containers::yaml& cfg);
        void save_to(mrpt::containers::yaml& cfg) const;
    };

    Parameters parameters;

    /** One inter-keyframe interval: the preintegrated IMU delta plus the
     *  external odometry's attitude and velocity at both ends, all expressed
     *  in the same (map/odometry) frame.
     */
    struct Interval
    {
        Interval() = default;

        TimeStamp t_from = 0;
        TimeStamp t_to   = 0;

        /// Preintegrated over `(t_from, t_to]`.
        PreintegratedImuMeasurements delta;

        SO3 R_from = SO3::Identity();  //!< odometry attitude (map <- body) at t_from
        SO3 R_to   = SO3::Identity();  //!< odometry attitude (map <- body) at t_to

        LinearVelocity v_from{0, 0, 0};  //!< map-frame velocity at t_from
        LinearVelocity v_to{0, 0, 0};  //!< map-frame velocity at t_to

        /// Map-frame covariances of the velocities above. Defaults are
        /// deliberately loose; supply the real ones when available.
        mrpt::math::CMatrixDouble33 cov_v_from = defaultVelocityCov();
        mrpt::math::CMatrixDouble33 cov_v_to   = defaultVelocityCov();

        static mrpt::math::CMatrixDouble33 defaultVelocityCov();
    };

    /// The estimator's output.
    struct Result
    {
        Result() = default;

        TimeStamp timestamp = 0;  //!< stamp of the newest interval's end

        /// The estimated gravity vector in the map frame [m/s^2]. Points
        /// DOWN, i.e. close to `(0, 0, -9.81)` while the map is level.
        mrpt::math::TVector3D gravity_in_map{0, 0, -9.81};

        /// Rotation that takes the map frame to a gravity-aligned one, with no
        /// yaw component. Left-multiply an odometry attitude by this to obtain
        /// the gravity-corrected attitude:  `R_corrected = correction * R_odom`.
        SO3 correction = SO3::Identity();

        /// The pitch/roll of `correction` [rad], i.e. how tilted the map has
        /// become, and their "earned" sigmas propagated from the solved
        /// information matrix. Consumers should weight the estimate with these
        /// instead of a hand-tuned constant. Note that these are the ONLY
        /// measure of quality: a large sigma is a weak estimate, not an
        /// invalid one, and is not filtered out here (see the class docs).
        double pitch_correction = 0;
        double roll_correction  = 0;
        double pitch_sigma      = 0;
        double roll_sigma       = 0;

        /// Total tilt angle of `correction` [rad].
        double tilt = 0;

        LinearAcceleration bias_acc{0, 0, 0};
        AngularVelocity    bias_gyro{0, 0, 0};

        std::size_t num_intervals = 0;
        std::size_t iterations    = 0;

        /// Norm of the last Gauss-Newton increment, for diagnostics: it tells
        /// whether the iteration settled or ran out of iterations.
        double increment_norm = 0;

        /// True when the solve yielded a USABLE estimate: the window held at
        /// least `min_intervals_for_convergence` intervals and the information
        /// matrix could be inverted into finite sigmas. It is NOT a statement
        /// about accuracy: that is `pitch_sigma` / `roll_sigma`.
        bool converged = false;
    };

    /** Minimal rotation (no yaw) taking the "up" direction implied by a
     *  gravity vector onto `+Z`, i.e. the correction that levels a frame in
     *  which gravity is measured as `gravityInMap`.
     *  Returns the identity for a degenerate (near-zero) input.
     */
    static SO3 correction_from_gravity(const mrpt::math::TVector3D& gravityInMap);

    /** @name Residuals
     *  The two residuals the solver minimizes, with their analytic Jacobians.
     *  They are public so that the Jacobians can be checked against numerical
     *  differentiation, and so that solver and tests share ONE definition of
     *  the math instead of two that can drift apart.
     *  All Jacobians are plain Euclidean derivatives (the unknowns are three
     *  vectors in R^3), so numerical differentiation of these very functions
     *  is a valid reference.
     *  @{ */

    /** Gravity residual of one interval: how much the odometry velocity change
     *  disagrees with the preintegrated one under a candidate gravity vector,
     *
     *      r = (v_j - v_i) - g * dt_ij - R_i * dV_ij(b_a, b_g)
     *
     *  where `dV_ij` is bias-corrected to first order around the interval's
     *  own linearization point. This is the residual that makes the gravity
     *  DIRECTION observable under arbitrary motion.
     */
    static mrpt::math::TVector3D gravity_residual(
        const Interval& interval, const mrpt::math::TVector3D& gravityInMap,
        const LinearAcceleration& biasAcc, const AngularVelocity& biasGyro,
        const mrpt::optional_ref<mrpt::math::CMatrixDouble33>& J_gravity   = std::nullopt,
        const mrpt::optional_ref<mrpt::math::CMatrixDouble33>& J_bias_acc  = std::nullopt,
        const mrpt::optional_ref<mrpt::math::CMatrixDouble33>& J_bias_gyro = std::nullopt);

    /** Relative-rotation residual of one interval, comparing the preintegrated
     *  rotation against the external odometry's relative rotation,
     *
     *      r = Log( dR_ij(b_g)^T * R_i^T * R_j )
     *
     *  This is what makes the gyroscope bias observable.
     */
    static mrpt::math::TVector3D rotation_residual(
        const Interval& interval, const AngularVelocity& biasGyro,
        const mrpt::optional_ref<mrpt::math::CMatrixDouble33>& J_bias_gyro = std::nullopt);

    /** @} */

    /** Sets the center of the bias priors, e.g. from a calibration at rest
     *  (ImuInitialCalibrator). Defaults to zero.
     */
    void set_bias_prior(const LinearAcceleration& biasAcc, const AngularVelocity& biasGyro);

    /** Appends one interval to the window, dropping the oldest when full.
     *
     *  @return false if the interval was REJECTED by the input guards
     *          (insufficient IMU coverage, empty, non-finite, out of order,
     *          implausible duration). Rejecting is always safe here: the
     *          unknowns are global, so a skipped interval leaves nothing
     *          unconstrained.
     */
    bool add_interval(const Interval& interval);

    /// Reason the last add_interval() returned false (for logging).
    const std::string& last_rejection_reason() const { return last_rejection_reason_; }

    /** Runs Gauss-Newton over the current window.
     *  @return false if there is not enough data, or the solve failed.
     */
    bool solve();

    /// The result of the last successful solve(), if any.
    const std::optional<Result>& latest_result() const { return latest_result_; }

    /// Number of intervals currently in the window.
    std::size_t window_length() const { return intervals_.size(); }

    /// Clears the window, keeping parameters and bias priors.
    void reset();

   private:
    std::deque<Interval> intervals_;

    mrpt::math::TVector3D gravity_{0, 0, -9.81};
    LinearAcceleration    bias_acc_{0, 0, 0};
    AngularVelocity       bias_gyro_{0, 0, 0};

    LinearAcceleration bias_acc_prior_{0, 0, 0};
    AngularVelocity    bias_gyro_prior_{0, 0, 0};

    bool gravity_initialized_ = false;

    std::optional<Result> latest_result_;
    std::string           last_rejection_reason_;

    bool reject(const std::string& why);
};

}  // namespace mola::imu

/** Feature macro: the package provides MapGravityEstimator, which estimates
 *  the gravity vector in the odometry/map frame (and the IMU biases) from
 *  preintegrated IMU measurements aided by external poses and velocities.
 */
#define MOLA_IMU_PREINTEGRATION_HAS_MAP_GRAVITY_ESTIMATOR 1

/** Feature macro: `Result::converged` reports USABILITY only. The former
 *  `Parameters::max_tilt_sigma_deg` quality gate is gone: consumers judge the
 *  estimate by the earned `pitch_sigma` / `roll_sigma` instead of receiving a
 *  pre-filtered boolean.
 */
#define MOLA_IMU_PREINTEGRATION_MAP_GRAVITY_UNGATED_CONVERGENCE 1
