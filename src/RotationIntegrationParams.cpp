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
 * @file   RotationIntegrationParams.cpp
 * @brief  Parameters for angular velocity integration.
 * @author Jose Luis Blanco Claraco
 * @date   Sep 20, 2021
 */

#include <mola_imu_preintegration/RotationIntegrationParams.h>

using namespace mola;

void RotationIntegrationParams::load_from(const mrpt::containers::yaml& cfg)
{
    gyroBias = mrpt::math::TVector3D::FromVector(cfg["gyroBias"].toStdVector<double>());

    const auto poseQuat  = cfg["sensorLocationInVehicle"]["quaternion"].toStdVector<double>();
    const auto poseTrans = cfg["sensorLocationInVehicle"]["translation"].toStdVector<double>();
    ASSERT_EQUAL_(poseQuat.size(), 4U);
    ASSERT_EQUAL_(poseTrans.size(), 3U);

    auto pose = mrpt::poses::CPose3D::FromQuaternionAndTranslation(
        mrpt::math::CQuaternionDouble(poseQuat[3], poseQuat[0], poseQuat[1], poseQuat[2]),
        mrpt::math::TPoint3D::FromVector(poseTrans));

    if (pose != mrpt::poses::CPose3D::Identity())
    {
        // Store:
        sensorPose = pose;
    }
    else
    {
        // Leave as unasigned to reflect it's just I_{4,4}
    }
}
