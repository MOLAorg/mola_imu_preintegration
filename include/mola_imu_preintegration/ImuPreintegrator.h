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
 * @file   ImuPreintegrator.h
 * @brief  Full IMU preintegration on manifold (rotation + velocity + position).
 * @author Jose Luis Blanco Claraco
 * @date   Jul 22, 2026
 */
#pragma once

#include <mola_imu_preintegration/ImuIntegrationParams.h>
#include <mola_imu_preintegration/so3_jacobians.h>
#include <mola_imu_preintegration/types.h>
#include <mrpt/math/CMatrixFixed.h>

#include <cstddef>

namespace mola::imu
{
/** Covariance of the preintegrated measurement, ordered as
 *  `[theta(3) , v(3) , p(3)]`.
 */
using PreintegrationCovariance = mrpt::math::CMatrixFixed<double, 9, 9>;

/** The result of preintegrating a set of IMU readings over an interval (i,j],
 *  around a fixed bias linearization point.
 *
 *  The deltas are defined, as usual (Forster et al. 2015), in the body frame
 *  of the interval start `i`:
 *
 *      dR_ij = R_i^T R_j
 *      dV_ij = R_i^T (v_j - v_i - g * dt_ij)
 *      dP_ij = R_i^T (p_j - p_i - v_i*dt_ij - 0.5*g*dt_ij^2)
 *
 *  so that gravity never enters the integration itself: it is introduced only
 *  when the delta is compared against a pair of states.
 *
 *  \ingroup mola_imu_preintegration_grp
 */
struct PreintegratedImuMeasurements
{
    PreintegratedImuMeasurements() = default;

    /// Time span covered by the integrated samples [s]. NOTE: this is the sum
    /// of the actual sample dt's, which is NOT necessarily the wall-clock
    /// duration of the interval (see the coverage guard in
    /// MapGravityEstimator).
    double deltaTij = 0;

    SO3                   deltaRij = SO3::Identity();  //!< Preintegrated rotation
    LinearVelocity        deltaVij{0, 0, 0};  //!< Preintegrated velocity
    mrpt::math::TVector3D deltaPij{0, 0, 0};  //!< Preintegrated position

    /** \name First-order bias-update Jacobians
     *  Defined such that, for a bias increment `db` w.r.t. the linearization
     *  point, `dR(b) ~= dR_bar * Exp(dR_dbg * db_g)`,
     *  `dV(b) ~= dV_bar + dV_dba * db_a + dV_dbg * db_g`, and likewise for dP.
     *  @{ */
    mrpt::math::CMatrixDouble33 dR_dbg = mrpt::math::CMatrixDouble33::Zero();
    mrpt::math::CMatrixDouble33 dV_dba = mrpt::math::CMatrixDouble33::Zero();
    mrpt::math::CMatrixDouble33 dV_dbg = mrpt::math::CMatrixDouble33::Zero();
    mrpt::math::CMatrixDouble33 dP_dba = mrpt::math::CMatrixDouble33::Zero();
    mrpt::math::CMatrixDouble33 dP_dbg = mrpt::math::CMatrixDouble33::Zero();
    /** @} */

    /// Propagated covariance of the preintegrated measurement, ordered
    /// `[theta, v, p]`. Only the noise of the raw readings is propagated here;
    /// bias uncertainty is handled by the estimator via the Jacobians above.
    PreintegrationCovariance cov = PreintegrationCovariance::Zero();

    /// The bias linearization point these deltas and Jacobians refer to.
    LinearAcceleration bias_acc{0, 0, 0};
    AngularVelocity    bias_gyro{0, 0, 0};

    /// Number of raw samples integrated.
    std::size_t num_measurements = 0;

    /// Bias-corrected deltas, to first order, for a new bias estimate.
    struct BiasCorrectedDelta
    {
        SO3                   deltaRij = SO3::Identity();
        LinearVelocity        deltaVij{0, 0, 0};
        mrpt::math::TVector3D deltaPij{0, 0, 0};
    };

    /** Applies the first-order bias update to obtain the deltas that would
     *  have resulted from integrating around `(newBiasAcc, newBiasGyro)`,
     *  WITHOUT re-integrating the raw samples. This is the whole point of
     *  preintegration.
     */
    BiasCorrectedDelta bias_corrected(
        const LinearAcceleration& newBiasAcc, const AngularVelocity& newBiasGyro) const;
};

/** Integrates accelerometer and gyroscope readings into a
 *  PreintegratedImuMeasurements, following:
 *
 *  Forster, C., Carlone, L., Dellaert, F., & Scaramuzza, D. (2017). On-manifold
 *  preintegration for real-time visual-inertial odometry. IEEE Transactions on
 *  Robotics, 33(1), 1-21.
 *
 *  Unlike ImuIntegrator (which is rotation-only), this class propagates the
 *  full `(dR, dV, dP)` triplet together with its covariance and the
 *  first-order bias Jacobians.
 *
 *  \note Readings MUST already be expressed in the vehicle/body frame,
 *        including any lever-arm compensation: use ImuTransformer for that.
 *        `parameters.sensor_pose` is therefore IGNORED by this class (an
 *        assertion enforces it is unset or the identity).
 *
 *  \note The accelerometer reading is the proper acceleration (specific
 *        force): at rest with the body Z axis up it reads `(0,0,+9.81)`.
 *        Gravity is deliberately NOT subtracted here, see
 *        PreintegratedImuMeasurements.
 *
 *  \note Discretization convention: each reading is held CONSTANT over the
 *        FOLLOWING `dt`, and the rotation used to rotate the acceleration is
 *        the one at the START of the step. This is the standard Euler-forward
 *        form of Forster et al. (and of GTSAM's implementation), which keeps
 *        the state propagation, the bias Jacobians and the covariance exactly
 *        consistent with one another. Its `dV`/`dP` discretization error is
 *        first order in `dt`, which at any realistic IMU rate is far below the
 *        sensor's own noise and bias errors. `dR` is unaffected whenever the
 *        angular rate is (piecewise) constant.
 *
 *  Usage:
 *  - (1) set `parameters` (noise densities + the bias linearization point),
 *  - (2) reset_integration(),
 *  - (3) integrate_measurement() for every sample of the interval,
 *  - (4) take current_state(), then go back to (2) for the next interval.
 *
 *  \ingroup mola_imu_preintegration_grp
 */
class ImuPreintegrator
{
   public:
    ImuPreintegrator() = default;
    explicit ImuPreintegrator(const ImuIntegrationParams& p) : parameters(p) {}

    /// All public parameters, including the bias linearization point and the
    /// continuous-time noise densities `cov_gyro` / `cov_acc`.
    ImuIntegrationParams parameters;

    /** Resets the accumulated state, keeping the bias linearization point
     *  currently held in `parameters`.
     */
    void reset_integration();

    /** Resets the accumulated state and sets a new bias linearization point.
     */
    void reset_integration(const LinearAcceleration& biasAcc, const AngularVelocity& biasGyro);

    /** Accumulates one IMU sample, held constant over `dt` seconds.
     *  @param a  Proper acceleration in the body frame [m/s^2].
     *  @param w  Angular velocity in the body frame [rad/s].
     *  @param dt Time step [s]; must be > 0.
     */
    void integrate_measurement(const LinearAcceleration& a, const AngularVelocity& w, double dt);

    const PreintegratedImuMeasurements& current_state() const { return state_; }

   private:
    PreintegratedImuMeasurements state_;
};

}  // namespace mola::imu

/** Feature macro: the package provides full on-manifold IMU preintegration
 *  (rotation + velocity + position, with covariance propagation and
 *  first-order bias Jacobians) via ImuPreintegrator, plus the SO(3) right
 *  Jacobian helpers in so3_jacobians.h. Downstream packages in separate repos
 *  should guard usage with
 *  `#if defined(MOLA_IMU_PREINTEGRATION_HAS_FULL_PREINTEGRATION)` to remain
 *  buildable against older checkouts.
 */
#define MOLA_IMU_PREINTEGRATION_HAS_FULL_PREINTEGRATION 1
