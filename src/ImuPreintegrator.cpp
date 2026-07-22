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
 * @file   ImuPreintegrator.cpp
 * @brief  Full IMU preintegration on manifold (rotation + velocity + position).
 * @author Jose Luis Blanco Claraco
 * @date   Jul 22, 2026
 */

#include <mola_imu_preintegration/ImuPreintegrator.h>
#include <mrpt/core/exceptions.h>
#include <mrpt/math/TPose3D.h>

#include <Eigen/Dense>

using namespace mola::imu;

namespace
{
/// Converts a TVector3D into an Eigen column vector.
Eigen::Vector3d toEigen(const mrpt::math::TVector3D& v) { return {v.x, v.y, v.z}; }

/// Converts an Eigen column vector into a TVector3D.
mrpt::math::TVector3D fromEigen(const Eigen::Vector3d& v) { return {v.x(), v.y(), v.z()}; }
}  // namespace

PreintegratedImuMeasurements::BiasCorrectedDelta PreintegratedImuMeasurements::bias_corrected(
    const LinearAcceleration& newBiasAcc, const AngularVelocity& newBiasGyro) const
{
    const Eigen::Vector3d dba = toEigen(newBiasAcc - bias_acc);
    const Eigen::Vector3d dbg = toEigen(newBiasGyro - bias_gyro);

    BiasCorrectedDelta out;

    // Rotation: a right-multiplied correction on the manifold.
    const auto dPhi        = fromEigen(Eigen::Vector3d(dR_dbg.asEigen() * dbg));
    out.deltaRij.asEigen() = deltaRij.asEigen() * so3_exp(dPhi).asEigen();

    // Velocity and position: plain vector-space first-order corrections.
    out.deltaVij =
        deltaVij + fromEigen(Eigen::Vector3d(dV_dba.asEigen() * dba + dV_dbg.asEigen() * dbg));
    out.deltaPij =
        deltaPij + fromEigen(Eigen::Vector3d(dP_dba.asEigen() * dba + dP_dbg.asEigen() * dbg));

    return out;
}

void ImuPreintegrator::reset_integration()
{
    reset_integration(parameters.bias_acc, parameters.bias_gyro);
}

void ImuPreintegrator::reset_integration(
    const LinearAcceleration& biasAcc, const AngularVelocity& biasGyro)
{
    state_           = PreintegratedImuMeasurements();
    state_.bias_acc  = biasAcc;
    state_.bias_gyro = biasGyro;

    // Keep `parameters` in sync, so that a caller reading them back sees the
    // linearization point actually in use:
    parameters.bias_acc  = biasAcc;
    parameters.bias_gyro = biasGyro;
}

void ImuPreintegrator::integrate_measurement(
    const LinearAcceleration& a, const AngularVelocity& w, double dt)
{
    ASSERT_GT_(dt, 0);
    ASSERTMSG_(
        !parameters.sensor_pose.has_value() ||
            parameters.sensor_pose->asTPose() == mrpt::math::TPose3D::Identity(),
        "ImuPreintegrator expects readings already in the body frame: "
        "parameters.sensor_pose must be unset or the identity (use "
        "ImuTransformer to move raw readings to the body frame).");

    // Bias-corrected readings, at the linearization point:
    const Eigen::Vector3d acc  = toEigen(a - state_.bias_acc);
    const Eigen::Vector3d gyro = toEigen(w - state_.bias_gyro);

    const mrpt::math::TVector3D phi     = fromEigen(Eigen::Vector3d(gyro * dt));
    const auto                  incrR   = so3_exp(phi);
    const auto                  Jr      = right_jacobian_so3(phi);
    const Eigen::Matrix3d       incrR_e = incrR.asEigen();
    const Eigen::Matrix3d       Jr_e    = Jr.asEigen();

    const Eigen::Matrix3d dR_e   = state_.deltaRij.asEigen();
    const Eigen::Matrix3d accHat = skew_symmetric(fromEigen(acc)).asEigen();

    // ---- Covariance propagation (order: [theta, v, p]) --------------------
    // Must run BEFORE the state is advanced, since A and B are built at the
    // current linearization point.
    {
        Eigen::Matrix<double, 9, 9> A = Eigen::Matrix<double, 9, 9>::Identity();
        A.block<3, 3>(0, 0)           = incrR_e.transpose();
        A.block<3, 3>(3, 0)           = -dR_e * accHat * dt;
        A.block<3, 3>(6, 0)           = -0.5 * dR_e * accHat * dt * dt;
        A.block<3, 3>(6, 3)           = Eigen::Matrix3d::Identity() * dt;

        Eigen::Matrix<double, 9, 3> Bg = Eigen::Matrix<double, 9, 3>::Zero();
        Bg.block<3, 3>(0, 0)           = Jr_e * dt;

        Eigen::Matrix<double, 9, 3> Ba = Eigen::Matrix<double, 9, 3>::Zero();
        Ba.block<3, 3>(3, 0)           = dR_e * dt;
        Ba.block<3, 3>(6, 0)           = 0.5 * dR_e * dt * dt;

        // cov_gyro / cov_acc are continuous-time noise DENSITIES, so the
        // discrete-time covariance over this step is density/dt.
        const Eigen::Matrix3d covG = parameters.cov_gyro.asEigen() / dt;
        const Eigen::Matrix3d covA = parameters.cov_acc.asEigen() / dt;

        state_.cov.asEigen() = A * state_.cov.asEigen() * A.transpose() +
                               Bg * covG * Bg.transpose() + Ba * covA * Ba.transpose();
    }

    // ---- Bias Jacobians ---------------------------------------------------
    // Also before advancing dR, for the same reason.
    {
        const Eigen::Matrix3d dR_dbg_prev = state_.dR_dbg.asEigen();
        const Eigen::Matrix3d dV_dba_prev = state_.dV_dba.asEigen();
        const Eigen::Matrix3d dV_dbg_prev = state_.dV_dbg.asEigen();

        state_.dP_dba.asEigen() += dV_dba_prev * dt - 0.5 * dR_e * dt * dt;
        state_.dP_dbg.asEigen() += dV_dbg_prev * dt - 0.5 * dR_e * accHat * dR_dbg_prev * dt * dt;

        state_.dV_dba.asEigen() -= dR_e * dt;
        state_.dV_dbg.asEigen() -= dR_e * accHat * dR_dbg_prev * dt;

        state_.dR_dbg.asEigen() = incrR_e.transpose() * dR_dbg_prev - Jr_e * dt;
    }

    // ---- State propagation ------------------------------------------------
    // Position first (it uses the pre-update velocity and rotation), then
    // velocity, then rotation.
    state_.deltaPij = state_.deltaPij + state_.deltaVij * dt +
                      fromEigen(Eigen::Vector3d(0.5 * dR_e * acc)) * dt * dt;
    state_.deltaVij           = state_.deltaVij + fromEigen(Eigen::Vector3d(dR_e * acc)) * dt;
    state_.deltaRij.asEigen() = dR_e * incrR_e;

    state_.deltaTij += dt;
    state_.num_measurements++;
}
