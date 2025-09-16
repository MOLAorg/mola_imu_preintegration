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

#include <mrpt/containers/yaml.h>
#include <mrpt/math/CMatrixFixed.h>
#include <mrpt/math/TPoint3D.h>
#include <mrpt/poses/CPose3D.h>

#include <optional>

namespace mola
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

    /// Gyroscope (initial or constant) bias, in the local IMU frame of reference (units: rad/s).
    mrpt::math::TVector3D gyroBias = {.0, .0, .0};

    /// Gyroscope covariance (units of sigma are rad/s/√Hz )
    mrpt::math::CMatrixDouble33 gyroCov = mrpt::math::CMatrixDouble33::Identity();

    /// Accelerometer (initial or constant) bias, in the local IMU frame of reference (units: m/s²).
    mrpt::math::TVector3D accelerometerBias = {.0, .0, .0};

    /// Gravity vector (units are m/s²), in the global frame of coordinates.
    mrpt::math::TVector3D gravityVector = {0, 0, -9.81};

    /// Accelerometer covariance (units of sigma are m/s²/√Hz )
    mrpt::math::CMatrixDouble33 accCov = mrpt::math::CMatrixDouble33::Identity();

    /// Integration covariance: jerk, that is, how much acceleration can change over time:
    mrpt::math::CMatrixDouble33 integrationCov = mrpt::math::CMatrixDouble33::Identity();

    /** If provided, defines an IMU placed at a pose different than the vehicle origin of
     * coordinates.
     * Default: IMU used as reference of the vehicle frame, i.e. sensorPose = SE(3) identity
     * I_{4x4}).
     */
    std::optional<mrpt::poses::CPose3D> sensorPose;
};

}  // namespace mola
