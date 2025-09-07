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
 * @file   ImuInitialCalibrator.cpp
 * @brief  Provides a rough initial calibration and attitude for IMUs without
 *         proper bias calibration
 * @author Jose Luis Blanco Claraco
 * @date   Sep 5, 2025
 */

#include <mola_imu_preintegration/ImuInitialCalibrator.h>
#include <mrpt/core/Clock.h>
#include <mrpt/core/exceptions.h>  // ASSERT_()
#include <mrpt/math/TPoint3D.h>
#include <mrpt/math/TTwist3D.h>

#include <sstream>

using namespace mola;

void ImuInitialCalibrator::add(const mrpt::obs::CObservationIMU::Ptr& obs)
{
    ASSERT_(obs);
    ASSERT_(parameters.required_samples > 2);
    ASSERT_(parameters.max_samples_age > 0);

    // Add IMU reading, after rotating it so it is body frame-referenced:
    samples_.emplace(
        mrpt::Clock::toDouble(obs->timestamp), imu_transformers_[obs->sensorLabel].process(*obs));

    // Remove old samples:
    while (!samples_.empty() &&
           samples_.begin()->first < samples_.rbegin()->first - parameters.max_samples_age)
    {
        samples_.erase(samples_.begin());
    }
}

bool ImuInitialCalibrator::isReady() const
{
    // samples enough?
    return (samples_.size() >= parameters.required_samples);
}

namespace
{
mrpt::math::TVector3D stddev_from_variances(
    const mrpt::math::TVector3D& variances, const double n_inverse)
{
    return {
        std::sqrt(variances.x * n_inverse),  //
        std::sqrt(variances.y * n_inverse),  //
        std::sqrt(variances.z * n_inverse)};
}
}  // namespace

std::optional<ImuInitialCalibrator::Results> ImuInitialCalibrator::getCalibration() const
{
    if (!isReady()) { return {}; }

    auto forEachAcc = [this](const std::function<void(const mrpt::math::TVector3D& acc)>& f)
    {
        for (const auto& [_, imu] : samples_)
        {
            const auto accel_base_link = mrpt::math::TVector3D(
                imu.get(mrpt::obs::IMU_X_ACC), imu.get(mrpt::obs::IMU_Y_ACC),
                imu.get(mrpt::obs::IMU_Z_ACC));
            f(accel_base_link);
        }
    };

    auto forEachGyro = [this](const std::function<void(const mrpt::math::TVector3D& omega)>& f)
    {
        for (const auto& [_, imu] : samples_)
        {
            const auto ang_vel_base_link = mrpt::math::TVector3D(
                imu.get(mrpt::obs::IMU_WX), imu.get(mrpt::obs::IMU_WY), imu.get(mrpt::obs::IMU_WZ));
            f(ang_vel_base_link);
        }
    };

    mrpt::math::TVector3D average_accel(0, 0, 0);
    forEachAcc([&](const auto& acc) { average_accel += acc; });

    mrpt::math::TVector3D average_gyro(0, 0, 0);
    forEachGyro([&](const auto& omega) { average_gyro += omega; });

    const std::size_t count = samples_.size();

    // Average:
    const auto n_1 = [count]() { return count ? 1.0 / static_cast<double>(count) : .0; }();
    average_accel *= n_1;
    average_gyro *= n_1;

    Results results;

    // Compute unbiased estimation of std. deviation (additive noise):
    const auto n_1b = [count]() { return count > 1 ? 1.0 / static_cast<double>(count - 1) : .0; }();

    {
        mrpt::math::TVector3D accel_variance(0, 0, 0);
        forEachAcc(
            [&](const mrpt::math::TVector3D& acc)
            {
                const auto err = acc - average_accel;
                accel_variance += {mrpt::square(err.x), mrpt::square(err.y), mrpt::square(err.z)};
            });
        results.noise_stddev_acc = stddev_from_variances(accel_variance, n_1b);
    }

    {
        mrpt::math::TVector3D gyro_variance(0, 0, 0);
        forEachGyro(
            [&](const mrpt::math::TVector3D& acc)
            {
                const auto err = acc - average_gyro;
                gyro_variance += {mrpt::square(err.x), mrpt::square(err.y), mrpt::square(err.z)};
            });
        results.noise_stddev_gyro = stddev_from_variances(gyro_variance, n_1b);
    }

    // Compute pitch & roll from the XYZ acceleration vector:
    const auto up_vector = average_accel.unitarize();

    results.pitch = -std::asin(up_vector.x);
    results.roll  = -std::asin(up_vector.y);

    // Accelerometer bias:
    const auto nominal_gravity_vector = mrpt::math::TVector3D(0., 0., parameters.gravity);

    results.bias_acc_b = (up_vector * parameters.gravity) - nominal_gravity_vector;

    // Gyroscope bias:
    results.bias_gyro = average_gyro;

    return results;
}

std::string ImuInitialCalibrator::Results::asString() const
{
    std::ostringstream oss;
    oss << "IMU Initial Calibration Results:\n";
    oss << "  Gyro bias: " << bias_gyro.asString() << "\n";
    oss << "  Acc bias (base_link): " << bias_acc_b.asString() << "\n";
    oss << "  Gyro noise stddev: " << noise_stddev_gyro.asString() << "\n";
    oss << "  Acc noise stddev: " << noise_stddev_acc.asString() << "\n";
    oss << "  Pitch: " << mrpt::RAD2DEG(pitch) << " deg\n";
    oss << "  Roll: " << mrpt::RAD2DEG(roll) << " deg\n";
    return oss.str();
}