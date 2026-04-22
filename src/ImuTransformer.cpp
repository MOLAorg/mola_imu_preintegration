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
 * @file   ImuTransformer.cpp
 * @brief  Transforms IMU readings from an arbitrary sensor pose to the vehicle reference base_link
 * @author Jose Luis Blanco Claraco
 * @date   Sep 5, 2025
 */

#include <mola_imu_preintegration/ImuTransformer.h>
#include <mrpt/math/TPoint3D.h>
#include <mrpt/math/TTwist3D.h>
#include <mrpt/math/geometry.h>
#include <mrpt/math/ops_containers.h>
#include <mrpt/obs/CObservationIMU.h>

using namespace mola::imu;

mrpt::obs::CObservationIMU ImuTransformer::process(const mrpt::obs::CObservationIMU& raw_imu)
{
    mrpt::obs::CObservationIMU imu = raw_imu;

    const bool hasAngVel = raw_imu.has(mrpt::obs::IMU_WX) && raw_imu.has(mrpt::obs::IMU_WY) &&
                           raw_imu.has(mrpt::obs::IMU_WZ);

    mrpt::math::TVector3D ang_vel_body = {0, 0, 0};
    // ang_acc is initialised to the current EMA state so that, if hasAngVel is
    // false, the last filtered estimate is reused rather than forcing zero.
    mrpt::math::TVector3D ang_acc = filtered_ang_acc_;

    if (hasAngVel)
    {
        // Transform angular velocity:
        const auto ang_vel_sensor = mrpt::math::TVector3D(  //
            raw_imu.get(mrpt::obs::IMU_WX),  //
            raw_imu.get(mrpt::obs::IMU_WY),  //
            raw_imu.get(mrpt::obs::IMU_WZ));

        ang_vel_body = imu.sensorPose.rotateVector(ang_vel_sensor);
        imu.set(mrpt::obs::IMU_WX, ang_vel_body.x);
        imu.set(mrpt::obs::IMU_WY, ang_vel_body.y);
        imu.set(mrpt::obs::IMU_WZ, ang_vel_body.z);

        // Estimate angular acceleration via finite difference, then apply an
        // exponential moving average (EMA) low-pass filter to suppress the
        // amplified gyro noise that finite differencing introduces.
        const auto this_stamp = mrpt::Clock::toDouble(raw_imu.timestamp);
        double     dt         = this_stamp - last_stamp_;
        if (dt <= 0 || dt > 1.0)
        {
            // It's either the first reading, an error, or data flow stopped and resumed.
            // Then use default rate:
            dt = 1.0 / 100.0;
        }

        if (first_sample_)
        {
            // On the very first call last_ang_vel_body_ is zero, so the raw
            // finite difference would produce a large spurious spike equal to
            // ang_vel_body / dt.  Bootstrap the filter state to zero instead
            // and skip the lever-arm correction for this sample.
            filtered_ang_acc_ = {0, 0, 0};
            filtered_ang_vel_ = ang_vel_body;  // bootstrap to current reading; no prior exists
            first_sample_     = false;
        }
        else
        {
            // Raw finite-difference angular acceleration:
            const auto raw_ang_acc = (ang_vel_body - last_ang_vel_body_) / dt;

            // EMA low-pass filter for angular acceleration:
            //   y[n] = alpha * x[n] + (1 - alpha) * y[n-1]
            // alpha close to 0: heavy smoothing; alpha = 1: no filtering.
            const double alpha_acc = parameters.ang_acc_lpf_alpha;
            filtered_ang_acc_ = raw_ang_acc * alpha_acc + filtered_ang_acc_ * (1.0 - alpha_acc);

            // EMA low-pass filter for angular velocity (used only in the
            // centripetal lever-arm term; the output channels stay unfiltered):
            const double alpha_vel = parameters.ang_vel_lpf_alpha;
            filtered_ang_vel_ = ang_vel_body * alpha_vel + filtered_ang_vel_ * (1.0 - alpha_vel);
        }

        ang_acc = filtered_ang_acc_;

        last_ang_vel_body_ = ang_vel_body;
        last_stamp_        = this_stamp;
    }

    // Transform acceleration:
    if (raw_imu.has(mrpt::obs::IMU_X_ACC) && raw_imu.has(mrpt::obs::IMU_Y_ACC) &&
        raw_imu.has(mrpt::obs::IMU_Z_ACC))
    {
        const auto raw_accel_sensor = mrpt::math::TVector3D(  //
            raw_imu.get(mrpt::obs::IMU_X_ACC),  //
            raw_imu.get(mrpt::obs::IMU_Y_ACC),  //
            raw_imu.get(mrpt::obs::IMU_Z_ACC));

        // Rotate acceleration from sensor frame to body frame:
        const auto accel_body_rotated = imu.sensorPose.rotateVector(raw_accel_sensor);

        //  a_imu  = a_body + α×t + ω×(ω×t)  ==>
        //  a_body = R * a_imu - α×t - ω×(ω×t)
        const auto t = raw_imu.sensorPose.translation();

        // On the first sample ang_acc is zero (see above), so the lever-arm
        // Euler term contributes nothing and no spike is introduced.
        // filtered_ang_vel_ is used for the centripetal term instead of the
        // raw ang_vel_body to suppress quadratic amplification of gyro noise;
        // the output IMU_WX/WY/WZ channels remain the unfiltered rotated reading.
        const auto accel_body =
            accel_body_rotated - mrpt::math::crossProduct3D(ang_acc, t) -
            mrpt::math::crossProduct3D(
                filtered_ang_vel_, mrpt::math::crossProduct3D(filtered_ang_vel_, t));

        imu.set(mrpt::obs::IMU_X_ACC, accel_body.x);
        imu.set(mrpt::obs::IMU_Y_ACC, accel_body.y);
        imu.set(mrpt::obs::IMU_Z_ACC, accel_body.z);
    }

    // Mark the new reference frame:
    imu.sensorPose = mrpt::poses::CPose3D::Identity();

    return imu;
}