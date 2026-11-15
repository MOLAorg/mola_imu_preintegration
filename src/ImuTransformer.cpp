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
#include <mrpt/core/exceptions.h>
#include <mrpt/math/TPoint3D.h>
#include <mrpt/math/TTwist3D.h>
#include <mrpt/math/geometry.h>
#include <mrpt/math/ops_containers.h>
#include <mrpt/obs/CObservationIMU.h>

#include <algorithm>
#include <cmath>

using namespace mola::imu;

namespace
{
// Span [s] of the first accelerometer readings whose mean norm decides their
// units, the mean norm below which they are taken in units of g, and g.
constexpr double ACCEL_UNITS_DETECTION_PERIOD = 0.5;
constexpr double ACCEL_NORM_MIN_FOR_SI_UNITS  = 3.0;
constexpr double ACCEL_NORM_MIN_VALID         = 0.3;
constexpr double STANDARD_GRAVITY             = 9.80665;  // [m/s²]

// Per-sample factor of a first-order EMA low-pass filter: a legacy fixed
// factor if given (> 0), otherwise the one matching the cut-off frequency for
// this time step. A non-positive cut-off disables the filter.
double ema_alpha(const double legacy_alpha, const double cutoff_hz, const double dt)
{
    if (legacy_alpha > 0)
    {
        return std::min(legacy_alpha, 1.0);
    }
    if (cutoff_hz <= 0)
    {
        return 1.0;
    }
    return 1.0 - std::exp(-2.0 * M_PI * cutoff_hz * dt);
}
}  // namespace

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
        const auto   this_stamp = mrpt::Clock::toDouble(raw_imu.timestamp);
        const double dt         = this_stamp - last_stamp_;

        // Sane range for the finite-difference step. The lower bound rejects
        // near-duplicate timestamps (some IMU drivers emit pairs of messages a
        // few microseconds apart), while accepting IMU rates up to 10 kHz.
        constexpr double MIN_SANE_DT = 1e-4;
        constexpr double MAX_SANE_DT = 1.0;

        if (first_sample_ || dt < 0 || dt > MAX_SANE_DT)
        {
            // First reading, time going backwards, or data flow stopped and
            // resumed: there is no usable previous sample, so a finite
            // difference would produce a spurious spike. Restart the filter
            // state from the current reading instead.
            filtered_ang_acc_  = {0, 0, 0};
            filtered_ang_vel_  = ang_vel_body;
            first_sample_      = false;
            last_ang_vel_body_ = ang_vel_body;
            last_stamp_        = this_stamp;
        }
        else if (dt < MIN_SANE_DT)
        {
            // Near-duplicate timestamp: dividing by such a tiny dt would turn a
            // negligible angular-velocity change into a huge acceleration
            // spike. Keep both the filter state and the reference sample, so
            // the next difference spans a sane interval.
        }
        else
        {
            // Raw finite-difference angular acceleration:
            const auto raw_ang_acc = (ang_vel_body - last_ang_vel_body_) / dt;

            // EMA low-pass filter for angular acceleration:
            //   y[n] = alpha * x[n] + (1 - alpha) * y[n-1]
            // alpha close to 0: heavy smoothing; alpha = 1: no filtering.
            const double alpha_acc =
                ema_alpha(parameters.ang_acc_lpf_alpha, parameters.ang_acc_lpf_cutoff_hz, dt);
            filtered_ang_acc_ = raw_ang_acc * alpha_acc + filtered_ang_acc_ * (1.0 - alpha_acc);

            // EMA low-pass filter for angular velocity (used only in the
            // centripetal lever-arm term; the output channels stay unfiltered):
            const double alpha_vel =
                ema_alpha(parameters.ang_vel_lpf_alpha, parameters.ang_vel_lpf_cutoff_hz, dt);
            filtered_ang_vel_ = ang_vel_body * alpha_vel + filtered_ang_vel_ * (1.0 - alpha_vel);

            last_ang_vel_body_ = ang_vel_body;
            last_stamp_        = this_stamp;
        }

        ang_acc = filtered_ang_acc_;
    }

    // Transform acceleration:
    if (raw_imu.has(mrpt::obs::IMU_X_ACC) && raw_imu.has(mrpt::obs::IMU_Y_ACC) &&
        raw_imu.has(mrpt::obs::IMU_Z_ACC))
    {
        const auto raw_accel_sensor = mrpt::math::TVector3D(  //
            raw_imu.get(mrpt::obs::IMU_X_ACC),  //
            raw_imu.get(mrpt::obs::IMU_Y_ACC),  //
            raw_imu.get(mrpt::obs::IMU_Z_ACC));

        // Units: some drivers publish the acceleration in units of g. If not
        // given, they are decided from the mean norm of the valid readings over
        // their first span: the decision so far applies within it, and is then
        // kept for good.
        ASSERT_GE_(parameters.accel_scale, 0.0);
        double accel_scale = parameters.accel_scale;
        if (accel_scale == 0)
        {
            const double stamp = mrpt::Clock::toDouble(raw_imu.timestamp);
            const double norm  = raw_accel_sensor.norm();
            // Readings without information (e.g. zeros while a driver starts)
            // are left out:
            const bool valid = std::isfinite(norm) && norm > ACCEL_NORM_MIN_VALID;
            if (!accel_units_decided_ && valid &&
                (accel_norm_count_ == 0 || stamp >= accel_first_stamp_))
            {
                if (accel_norm_count_ == 0)
                {
                    accel_first_stamp_ = stamp;
                }
                accel_norm_sum_ += norm;
                accel_norm_count_++;
                const double meanNorm = accel_norm_sum_ / static_cast<double>(accel_norm_count_);
                accel_auto_scale_ = meanNorm < ACCEL_NORM_MIN_FOR_SI_UNITS ? STANDARD_GRAVITY : 1.0;
                accel_units_decided_ = stamp - accel_first_stamp_ >= ACCEL_UNITS_DETECTION_PERIOD;
            }
            accel_scale = accel_auto_scale_;
        }

        // Rotate acceleration from sensor frame to body frame:
        const auto accel_body_rotated = imu.sensorPose.rotateVector(raw_accel_sensor * accel_scale);

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