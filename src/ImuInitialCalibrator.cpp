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
#include <mrpt/math/CMatrixFixed.h>
#include <mrpt/math/CQuaternion.h>
#include <mrpt/math/TPoint3D.h>
#include <mrpt/math/TTwist3D.h>
#include <mrpt/obs/CObservationIMU.h>
#include <mrpt/poses/CPose3D.h>
#include <mrpt/poses/SO_SE_average.h>

#include <Eigen/Dense>  // required by MRPT matrix transpose()/operator* below
#include <algorithm>
#include <cmath>
#include <exception>
#include <functional>
#include <iterator>
#include <sstream>

using namespace mola::imu;

double ImuInitialCalibrator::bufferHorizon() const
{
    if (parameters.window_seconds <= 0)
    {
        return parameters.max_samples_age;
    }
    // Comfortably more than the window: pruning exactly at the window would keep the buffer
    // shorter than it, so "the window is complete" could only ever be true by an exact
    // floating-point coincidence, and at low sample rates never at all. Only the newest
    // `window_seconds` take part in the average, see windowBegin().
    return 2.0 * parameters.window_seconds;
}

ImuInitialCalibrator::const_iterator ImuInitialCalibrator::windowBegin() const
{
    if (parameters.window_seconds <= 0 || samples_.empty())
    {
        return samples_.begin();
    }
    return samples_.lower_bound(samples_.rbegin()->first - parameters.window_seconds);
}

void ImuInitialCalibrator::add(const mrpt::obs::CObservationIMU::ConstPtr& obs)
{
    ASSERT_(obs);
    ASSERT_(parameters.required_samples > 2);
    ASSERT_(bufferHorizon() > 0);

    // Rotate accel/gyro so they are body (base_link) frame-referenced:
    mrpt::obs::CObservationIMU bodyImu = imu_transformers_[obs->sensorLabel].process(*obs);

    // ImuTransformer only rotates accel/gyro; the absolute-orientation quaternion
    // (if present) is left in the SENSOR frame. Bring it to the vehicle frame too,
    // R_world_vehicle = R_world_sensor * R_vehicle_sensor^T, so the orientation-
    // derived pitch/roll in getCalibration() describe the VEHICLE and not the
    // (possibly rotated) sensor. Also DROP a PLACEHOLDER identity quaternion: some
    // IMUs that do NOT estimate attitude (e.g. the Hesai built-in IMU) still ship
    // a bit-exact identity quaternion tagged orientation_covariance[0]=0, which
    // the ROS bridge forwards as "valid". It is not a real attitude, and through a
    // rotated mount it would fabricate a bogus pitch/roll; dropped, the (already
    // body-frame) accelerometer-gravity path does the leveling instead.
    if (bodyImu.has(mrpt::obs::IMU_ORI_QUAT_W))
    {
        const auto qw = bodyImu.get(mrpt::obs::IMU_ORI_QUAT_W);
        const auto qx = bodyImu.get(mrpt::obs::IMU_ORI_QUAT_X);
        const auto qy = bodyImu.get(mrpt::obs::IMU_ORI_QUAT_Y);
        const auto qz = bodyImu.get(mrpt::obs::IMU_ORI_QUAT_Z);

        const bool isPlaceholderIdentity = std::abs(qx) < 1e-6 && std::abs(qy) < 1e-6 &&
                                           std::abs(qz) < 1e-6 &&
                                           std::abs(std::abs(qw) - 1.0) < 1e-6;

        // Some drivers that do not estimate attitude at all (e.g. a raw MEMS
        // IMU with no onboard AHRS) leave the orientation quaternion at all
        // zeros instead of a placeholder identity. That is not a rotation
        // (norm 0), and would otherwise abort below when handed to
        // CQuaternion; drop it the same way as a placeholder identity.
        // The comparison is negated (rather than "> 0.1") so a NaN qNormSq
        // (e.g. from a sensor anomaly) is also classified as degenerate:
        // NaN fails every relational comparison, so "> 0.1" alone would let
        // it slip through as if it were a valid unit quaternion.
        const double qNormSq =
            mrpt::square(qw) + mrpt::square(qx) + mrpt::square(qy) + mrpt::square(qz);
        const bool isDegenerate = !(std::abs(qNormSq - 1.0) <= 0.1);

        if (isPlaceholderIdentity || isDegenerate)
        {
            bodyImu.dataIsPresent.at(mrpt::obs::IMU_ORI_QUAT_W) = false;
            bodyImu.dataIsPresent.at(mrpt::obs::IMU_ORI_QUAT_X) = false;
            bodyImu.dataIsPresent.at(mrpt::obs::IMU_ORI_QUAT_Y) = false;
            bodyImu.dataIsPresent.at(mrpt::obs::IMU_ORI_QUAT_Z) = false;
        }
        else
        {
            // obs->sensorPose still holds R_vehicle_sensor (process() reset the
            // sensor pose on the returned COPY only, not on the input obs):
            const mrpt::math::CMatrixDouble33 R_world_sensor =
                mrpt::poses::CPose3D::FromQuaternion(mrpt::math::CQuaternionDouble(qw, qx, qy, qz))
                    .getRotationMatrix();
            const mrpt::math::CMatrixDouble33 R_vehicle_sensor =
                obs->sensorPose.getRotationMatrix();

            mrpt::math::CMatrixDouble33 R_world_vehicle;
            R_world_vehicle.asEigen() =
                R_world_sensor.asEigen() * R_vehicle_sensor.asEigen().transpose();

            mrpt::poses::CPose3D pv;
            pv.setRotationMatrix(R_world_vehicle);
            mrpt::math::CQuaternionDouble qv;
            pv.getAsQuaternion(qv);
            bodyImu.set(mrpt::obs::IMU_ORI_QUAT_W, qv.w());
            bodyImu.set(mrpt::obs::IMU_ORI_QUAT_X, qv.x());
            bodyImu.set(mrpt::obs::IMU_ORI_QUAT_Y, qv.y());
            bodyImu.set(mrpt::obs::IMU_ORI_QUAT_Z, qv.z());
        }
    }

    // Add IMU reading (now fully body frame-referenced, orientation included):
    const double t = mrpt::Clock::toDouble(obs->timestamp);
    samples_.emplace(t, std::move(bodyImu));

    if (!first_sample_time_.has_value())
    {
        first_sample_time_ = t;
    }

    // Remove old samples:
    const double horizon = bufferHorizon();
    while (!samples_.empty() && samples_.begin()->first < samples_.rbegin()->first - horizon)
    {
        samples_.erase(samples_.begin());
    }
}

std::optional<double> ImuInitialCalibrator::directionDispersion() const
{
    const auto itFirst = windowBegin();

    // Mean direction of the accelerometer readings:
    mrpt::math::TVector3D mean(0, 0, 0);
    std::size_t           count = 0;

    auto forEachDirection = [&](const std::function<void(const mrpt::math::TVector3D&)>& f)
    {
        for (auto it = itFirst; it != samples_.end(); ++it)
        {
            const auto& imu = it->second;
            const auto  acc = mrpt::math::TVector3D(
                 imu.get(mrpt::obs::IMU_X_ACC), imu.get(mrpt::obs::IMU_Y_ACC),
                 imu.get(mrpt::obs::IMU_Z_ACC));
            const double n = acc.norm();
            if (n < 1e-6)
            {
                continue;
            }
            f(acc * (1.0 / n));
        }
    };

    forEachDirection(
        [&](const mrpt::math::TVector3D& u)
        {
            mean += u;
            ++count;
        });

    if (count < 2)
    {
        return {};
    }

    const double meanNorm = mean.norm();
    if (meanNorm < 1e-6)
    {
        return {};
    }
    mean *= 1.0 / meanNorm;

    // RMS angle of each direction about the mean direction:
    double sumSqAngle = 0;
    forEachDirection(
        [&](const mrpt::math::TVector3D& u)
        {
            const double c = std::clamp(u.x * mean.x + u.y * mean.y + u.z * mean.z, -1.0, 1.0);
            sumSqAngle += mrpt::square(std::acos(c));
        });

    return std::sqrt(sumSqAngle / static_cast<double>(count));
}

ImuInitialCalibrator::Readiness ImuInitialCalibrator::readiness() const
{
    Readiness r;

    const auto itFirst = windowBegin();
    r.samples          = static_cast<std::size_t>(std::distance(itFirst, samples_.end()));

    if (samples_.empty())
    {
        r.reason = "no IMU samples yet";
        return r;
    }

    const double newest = samples_.rbegin()->first;

    // Reported span is the one that will actually be averaged; the completeness test below uses
    // the whole buffer instead. Samples inside the window are by construction no older than
    // `window_seconds`, so their span can only reach it by exact floating-point coincidence.
    r.span    = newest - itFirst->first;
    r.elapsed = newest - first_sample_time_.value_or(itFirst->first);

    // Minimum-sample sanity floor. It only rejects a degenerate handful of samples: what sets
    // the averaging time is window_seconds, since the error is dominated by platform motion
    // during the window and not by sensor noise.
    if (r.samples < parameters.required_samples)
    {
        r.reason = "not enough samples yet";
        return r;
    }

    if (parameters.window_seconds > 0 &&
        (newest - samples_.begin()->first) < parameters.window_seconds)
    {
        r.reason = "averaging window not filled yet";
        return r;
    }

    r.dispersion = directionDispersion();

    if (parameters.max_direction_dispersion > 0 && r.dispersion.has_value() &&
        *r.dispersion > parameters.max_direction_dispersion)
    {
        // The window is not measuring gravity: freezing it would bake platform motion into the
        // attitude. Defer, unless we have been deferring for too long already, in which case a
        // motion-contaminated estimate still beats never initializing.
        if (parameters.dispersion_timeout <= 0 || r.elapsed < parameters.dispersion_timeout)
        {
            r.reason = "accelerometer dispersion too high (platform not still)";
            return r;
        }
        r.timed_out = true;
    }

    r.ready = true;
    return r;
}

bool ImuInitialCalibrator::isReady() const { return readiness().ready; }

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
    if (!isReady())
    {
        return {};
    }

    // Only the samples inside the averaging window take part in the calibration: in time-window
    // mode the buffer may hold slightly more than that.
    const auto itFirst = windowBegin();

    auto forEachAcc = [&](const std::function<void(const mrpt::math::TVector3D& acc)>& f)
    {
        for (auto it = itFirst; it != samples_.end(); ++it)
        {
            const auto& imu             = it->second;
            const auto  accel_base_link = mrpt::math::TVector3D(
                 imu.get(mrpt::obs::IMU_X_ACC), imu.get(mrpt::obs::IMU_Y_ACC),
                 imu.get(mrpt::obs::IMU_Z_ACC));
            f(accel_base_link);
        }
    };

    auto forEachGyro = [&](const std::function<void(const mrpt::math::TVector3D& omega)>& f)
    {
        for (auto it = itFirst; it != samples_.end(); ++it)
        {
            const auto& imu               = it->second;
            const auto  ang_vel_base_link = mrpt::math::TVector3D(
                 imu.get(mrpt::obs::IMU_WX), imu.get(mrpt::obs::IMU_WY), imu.get(mrpt::obs::IMU_WZ));
            f(ang_vel_base_link);
        }
    };

    auto forEachOrientation = [&](const std::function<void(const mrpt::poses::CPose3D&)>& f)
    {
        for (auto it = itFirst; it != samples_.end(); ++it)
        {
            const auto& imu = it->second;
            if (!imu.has(mrpt::obs::IMU_ORI_QUAT_W))
            {
                continue;
            }
            const auto qw = imu.get(mrpt::obs::IMU_ORI_QUAT_W);
            const auto qx = imu.get(mrpt::obs::IMU_ORI_QUAT_X);
            const auto qy = imu.get(mrpt::obs::IMU_ORI_QUAT_Y);
            const auto qz = imu.get(mrpt::obs::IMU_ORI_QUAT_Z);
            const auto q_norm =
                mrpt::square(qw) + mrpt::square(qx) + mrpt::square(qy) + mrpt::square(qz);
            if (std::abs(q_norm - 1.0) > 0.1)
            {
                // Orientation data seems invalid!
                continue;
            }

            const auto p =
                mrpt::poses::CPose3D::FromQuaternion(mrpt::math::CQuaternionDouble(qw, qx, qy, qz));
            f(p);
        }
    };

    mrpt::math::TVector3D average_accel(0, 0, 0);
    forEachAcc([&](const auto& acc) { average_accel += acc; });

    mrpt::math::TVector3D average_gyro(0, 0, 0);
    forEachGyro([&](const auto& omega) { average_gyro += omega; });

    std::optional<mrpt::poses::CPose3D> avr_so3;
    if (parameters.use_imu_orientation)
    {
        mrpt::poses::SO_average<3> so3_average;
        forEachOrientation([&](const auto& pose) { so3_average.append(pose.getRotationMatrix()); });
        try
        {
            auto rot = so3_average.get_average();

            avr_so3.emplace();
            avr_so3->setRotationMatrix(rot);
        }
        catch (const std::exception& e)
        {
            // Ignore, this means we had no data.
            (void)e;
        }
    }

    const auto count = static_cast<std::size_t>(std::distance(itFirst, samples_.end()));

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

    if (avr_so3.has_value())
    {
        results.pitch = avr_so3->pitch();
        results.roll  = avr_so3->roll();
    }
    else
    {
        results.pitch = -std::asin(up_vector.x);
        results.roll  = std::atan2(up_vector.y, up_vector.z);
    }

    // Accelerometer bias:
    const auto estimated_gravity_body = up_vector * parameters.gravity;
    const auto average_accel_body     = average_accel;

    results.bias_acc_b = average_accel_body - estimated_gravity_body;

    // Gyroscope bias:
    // WARNING: this is the plain mean of the gyro, which is only a bias if the platform was
    // actually still during the window. On a moving start it is real rotation rate, so it must
    // not be fed to a preintegrator without a staticness gate first.
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