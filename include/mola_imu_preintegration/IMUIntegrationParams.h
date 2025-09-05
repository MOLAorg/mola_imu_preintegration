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
 * @file   IMUIntegrationParams.h
 * @brief  Parameters for IMU preintegration.
 * @author Jose Luis Blanco Claraco
 * @date   Sep 19, 2021
 */

#pragma once

#include <mola_imu_preintegration/RotationIntegrationParams.h>
#include <mrpt/containers/yaml.h>
#include <mrpt/math/CMatrixFixed.h>
#include <mrpt/math/TPoint3D.h>

namespace mola
{
/** Parameters needed by IMU preintegration classes, when integrating both,
 * acceleration and rotation.
 *
 *  Refer to:
 *  - Crassidis, J. L. (2006). Sigma-point Kalman filtering for integrated GPS
 * and inertial navigation. IEEE Transactions on Aerospace and Electronic
 * Systems, 42(2), 750-756.
 *  - Forster, C., Carlone, L., Dellaert, F., & Scaramuzza, D. (2015). IMU
 * preintegration on manifold for efficient visual-inertial maximum-a-posteriori
 * estimation. Georgia Institute of Technology.
 *  - Nikolic, J. (2016). Characterisation, calibration, and design of
 * visual-inertial sensor systems for robot navigation (Doctoral dissertation,
 * ETH Zurich).
 *
 * \ingroup mola_imu_preintegration_grp
 */
class IMUIntegrationParams
{
   public:
    IMUIntegrationParams() = default;

    /// Loads all parameters from a YAML map node.
    void loadFrom(const mrpt::containers::yaml& cfg);

    /// Parameters for gyroscope integration:
    RotationIntegrationParams rotationParams;

    /// Gravity vector (units are m/s²), in the global frame of coordinates.
    mrpt::math::TVector3D gravityVector = {0, 0, -9.81};

    /// Accelerometer covariance (units of sigma are m/s²/√Hz )
    mrpt::math::CMatrixDouble33 accCov = mrpt::math::CMatrixDouble33::Identity();

    /// Integration covariance: jerk, that is, how much acceleration can change
    /// over time:
    mrpt::math::CMatrixDouble33 integrationCov = mrpt::math::CMatrixDouble33::Identity();
};

}  // namespace mola
