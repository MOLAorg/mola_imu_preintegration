/*               _
 _ __ ___   ___ | | __ _
| '_ ` _ \ / _ \| |/ _` | Modular Optimization framework for
| | | | | | (_) | | (_| | Localization and mApping (MOLA)
|_| |_| |_|\___/|_|\__,_| https://github.com/MOLAorg/mola

 Copyright (C) 2018-2025 Jose Luis Blanco, University of Almeria,
                         and individual contributors.
 SPDX-License-Identifier: GPL-3.0
 See LICENSE for full license information.
 Closed-source licenses available upon request, for this odometry package
 alone or in combination with the complete SLAM system.
*/

/**
 * @file   ImuIntegrationParams.h
 * @brief  Parameters for IMU preintegration.
 * @author Jose Luis Blanco Claraco
 * @date   Sep 19, 2021
 */

#pragma once

#include <mola_imu_preintegration/types.h>
#include <mrpt/containers/yaml.h>
#include <mrpt/math/CMatrixFixed.h>
#include <mrpt/math/TPoint3D.h>
#include <mrpt/poses/CPose3D.h>

#include <optional>

namespace mola::imu
{
/** Parameters needed by IMU integration routines to integrate acceleration and rotation.
 *
 * Refer to:
 *  - Crassidis, J. L. (2006). Sigma-point Kalman filtering for integrated GPS and inertial
 *    navigation. IEEE Transactions on Aerospace and Electronic Systems, 42(2), 750-756.
 *  - Forster, C., Carlone, L., Dellaert, F., & Scaramuzza, D. (2015).
 *    IMU preintegration on manifold for efficient visual-inertial maximum-a-posteriori estimation.
 *  - Nikolic, J. (2016). Characterisation, calibration, and design of visual-inertial sensor
 *    systems for robot navigation (Doctoral dissertation, ETH Zurich).
 *
 * \ingroup mola_imu_preintegration_grp
 */
class ImuIntegrationParams
{
   public:
    ImuIntegrationParams() = default;

    /// Loads all parameters from a YAML map node.
    void load_from(const mrpt::containers::yaml& cfg);

    /// Saves all parameters as YAML.
    void save_to(mrpt::containers::yaml& cfg) const;

    /// Gravity vector (units are m/s²), in the global gravity-aligned frame of coordinates.
    LinearAcceleration gravity_vector = {0, 0, -9.81};

    /// Gyroscope (initial or constant) bias, in the local IMU frame of reference (units: rad/s).
    AngularVelocity bias_gyro = {.0, .0, .0};

    /// Default gyroscope noise density [rad/s/√Hz], representative of an
    /// industrial-grade MEMS IMU (~0.01 deg/s/√Hz, e.g. ADIS16448 / Xsens MTi
    /// class). \sa cov_gyro
    constexpr static double DEFAULT_GYRO_NOISE_DENSITY = 1.7e-4;

    /// Default accelerometer noise density [m/s²/√Hz], representative of an
    /// industrial-grade MEMS IMU (~200 µg/√Hz). \sa cov_acc
    constexpr static double DEFAULT_ACCEL_NOISE_DENSITY = 2.0e-3;

    /// Builds the isotropic covariance matrix for a given noise density.
    [[nodiscard]] static mrpt::math::CMatrixDouble33 isotropic_cov(double sigma)
    {
        mrpt::math::CMatrixDouble33 m;
        m.setDiagonal(sigma * sigma);
        return m;
    }

    /** Gyroscope covariance (units of sigma are rad/s/√Hz).
     *
     *  Defaults to DEFAULT_GYRO_NOISE_DENSITY squared, i.e. a typical MEMS
     *  datasheet figure. Set the real value for your sensor: field deployments
     *  commonly inflate the datasheet density several-fold to absorb vibration
     *  and unmodeled dynamics.
     *
     *  \note This is a continuous-time noise DENSITY, not a per-sample
     *        covariance: consumers scale it by the sample period.
     */
    mrpt::math::CMatrixDouble33 cov_gyro = isotropic_cov(DEFAULT_GYRO_NOISE_DENSITY);

    /// Accelerometer (initial or constant) bias, in the local IMU frame of reference (units: m/s²).
    LinearAcceleration bias_acc = {.0, .0, .0};

    /** Accelerometer covariance (units of sigma are m/s²/√Hz).
     *  Defaults to DEFAULT_ACCEL_NOISE_DENSITY squared. See cov_gyro for the
     *  units convention and the caveat about inflating datasheet figures.
     */
    mrpt::math::CMatrixDouble33 cov_acc = isotropic_cov(DEFAULT_ACCEL_NOISE_DENSITY);

    /// Integration covariance: jerk, that is, how much acceleration can change over time:
    // mrpt::math::CMatrixDouble33 cov_integration = mrpt::math::CMatrixDouble33::Identity();

    /** If provided, defines an IMU placed at a pose different than the vehicle origin of
     * coordinates.
     * Default: IMU used as reference of the vehicle frame, i.e. sensor_pose = SE(3) identity
     * I_{4x4}).
     */
    std::optional<mrpt::poses::CPose3D> sensor_pose;
};

}  // namespace mola::imu
