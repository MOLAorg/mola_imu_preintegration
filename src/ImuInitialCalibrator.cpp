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
 * @file   ImuInitialCalibrator.h
 * @brief  Provides a rough initial calibration and attitude for IMUs without
 *         proper bias calibration
 * @author Jose Luis Blanco Claraco
 * @date   Sep 5, 2025
 */

#include <mola_imu_preintegration/ImuInitialCalibrator.h>
#include <mrpt/core/Clock.h>
#include <mrpt/math/TPoint3D.h>
#include <mrpt/math/TTwist3D.h>
#include <mrpt/obs/CObservationIMU.h>

using namespace mola;

void ImuInitialCalibrator::add(const std::shared_ptr<const mrpt::obs::CObservationIMU>& obs)
{
    // Add:
    ASSERT_(obs);
    samples_[mrpt::Clock::toDouble(obs->timestamp)] = obs;

    // Remove old samples:
    while (!samples_.empty() &&
           samples_.begin()->first < samples_.rbegin()->first - max_samples_age_)
    {
        samples_.erase(samples_.begin());
    }
}

bool ImuInitialCalibrator::isReady() const
{
    // samples enough?
    return (samples_.size() >= required_samples_);
}

std::tuple<double, double> ImuInitialCalibrator::getPitchRoll() const
{
    mrpt::math::TVector3D avr_accel(0, 0, 0);
    std::size_t           count = 0;

    for (const auto& [_, imu] : samples_)
    {
        ASSERT_(imu);
        // TODO: Check for direct IMU-provided pitch/roll values?

        const auto accel_sensor = mrpt::math::TTwist3D(  //
      imu->get(mrpt::obs::IMU_X_ACC),                //
      imu->get(mrpt::obs::IMU_Y_ACC),                //
      imu->get(mrpt::obs::IMU_Z_ACC),                //
      0, 0, 0);

        // TODO: Minimum sanity check for the acceleration vector?

        const auto accel_base_link = accel_sensor.rotated(imu->sensorPose.asTPose());

        // Accumulate:
        avr_accel +=
            mrpt::math::TVector3D(accel_base_link.vx, accel_base_link.vy, accel_base_link.vz);
        count++;
    }

    // Average:
    if (count > 0) { avr_accel *= 1.0 / static_cast<double>(count); }

    // Remove accelerometer bias:
    // avr_accel-= bias_accel;

    // Compute pitch & roll from the XYZ acceleration vector:
    const auto up_vector = avr_accel.unitarize();

    const double pitch = -std::asin(up_vector.x);
    const double roll  = -std::asin(up_vector.y);
    return {pitch, roll};
}
