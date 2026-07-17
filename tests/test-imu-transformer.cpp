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
 * @file   test-imu-transformer.cpp
 * @brief  Unit tests for ImuTransformer (sensor pose to base_link transform)
 * @author Jose Luis Blanco Claraco
 * @date   Jun 29, 2026
 */

#include <mola_imu_preintegration/ImuTransformer.h>
#include <mrpt/core/Clock.h>
#include <mrpt/core/exceptions.h>
#include <mrpt/math/TPoint3D.h>
#include <mrpt/math/geometry.h>
#include <mrpt/obs/CObservationIMU.h>
#include <mrpt/poses/CPose3D.h>

#include <cmath>
#include <iostream>

using namespace mola::imu;
using mrpt::math::TVector3D;

namespace
{
mrpt::obs::CObservationIMU::Ptr make_imu(
    double time_sec, const TVector3D& acc, const TVector3D& gyro,
    const mrpt::poses::CPose3D& sensorPose)
{
    auto obs        = mrpt::obs::CObservationIMU::Create();
    obs->timestamp  = mrpt::Clock::fromDouble(time_sec);
    obs->sensorPose = sensorPose;

    obs->set(mrpt::obs::IMU_X_ACC, acc.x);
    obs->set(mrpt::obs::IMU_Y_ACC, acc.y);
    obs->set(mrpt::obs::IMU_Z_ACC, acc.z);

    obs->set(mrpt::obs::IMU_WX, gyro.x);
    obs->set(mrpt::obs::IMU_WY, gyro.y);
    obs->set(mrpt::obs::IMU_WZ, gyro.z);

    return obs;
}

TVector3D get_acc(const mrpt::obs::CObservationIMU& o)
{
    return {o.get(mrpt::obs::IMU_X_ACC), o.get(mrpt::obs::IMU_Y_ACC), o.get(mrpt::obs::IMU_Z_ACC)};
}

TVector3D get_gyro(const mrpt::obs::CObservationIMU& o)
{
    return {o.get(mrpt::obs::IMU_WX), o.get(mrpt::obs::IMU_WY), o.get(mrpt::obs::IMU_WZ)};
}

void check_near(const TVector3D& a, const TVector3D& b, double tol, const char* what)
{
    if (std::abs(a.x - b.x) > tol || std::abs(a.y - b.y) > tol || std::abs(a.z - b.z) > tol)
    {
        THROW_EXCEPTION_FMT(
            "%s mismatch: got %s, expected %s", what, a.asString().c_str(), b.asString().c_str());
    }
}

// 1) Angular velocity is rotated from sensor to body frame; with no lever arm,
//    acceleration is only rotated.
void test_rotation_only()
{
    // Sensor rotated +90 deg yaw, no translation.
    const auto sensorPose = mrpt::poses::CPose3D::FromYawPitchRoll(mrpt::DEG2RAD(90.0), 0, 0);

    ImuTransformer  tf;
    const TVector3D gyro_sensor = {1.0, 0.0, 0.0};
    const TVector3D acc_sensor  = {0.0, 0.0, 9.81};

    const auto out = tf.process(*make_imu(100.0, acc_sensor, gyro_sensor, sensorPose));

    // R(+90 yaw) maps x->y. So body gyro = (0,1,0); accel z stays z.
    check_near(get_gyro(out), {0.0, 1.0, 0.0}, 1e-9, "rotated gyro");
    check_near(get_acc(out), {0.0, 0.0, 9.81}, 1e-9, "rotated accel");

    // The relabeled frame must be the identity:
    ASSERT_(out.sensorPose == mrpt::poses::CPose3D::Identity());

    std::cout << "✅ test_rotation_only passed!" << std::endl;
}

// 2) Centripetal lever-arm term -w x (w x t) for constant angular velocity.
void test_centripetal_lever_arm()
{
    // Lever arm t=(0.5,0,0), identity rotation.
    const auto sensorPose = mrpt::poses::CPose3D::FromTranslation(0.5, 0.0, 0.0);

    ImuTransformer tf;
    // No filtering, so the term uses the full angular velocity immediately:
    tf.parameters.ang_acc_lpf_alpha = 1.0;
    tf.parameters.ang_vel_lpf_alpha = 1.0;

    const TVector3D w   = {0.0, 0.0, 2.0};  // constant => zero angular accel
    const TVector3D acc = {0.0, 0.0, 9.81};

    // Feed a few constant samples (first one bootstraps filter state).
    mrpt::obs::CObservationIMU out;
    for (int i = 0; i < 5; ++i)
    {
        out = tf.process(*make_imu(100.0 + 0.01 * i, acc, w, sensorPose));
    }

    // w x (w x t): w=(0,0,2), t=(0.5,0,0) => (-2,0,0). a_body = acc - that = (2,0,9.81).
    check_near(get_acc(out), {2.0, 0.0, 9.81}, 1e-9, "centripetal corrected accel");

    std::cout << "✅ test_centripetal_lever_arm passed!" << std::endl;
}

// 3) Tangential lever-arm term -alpha x t for a linear angular-velocity ramp.
void test_tangential_lever_arm()
{
    const auto      sensorPose = mrpt::poses::CPose3D::FromTranslation(0.5, 0.0, 0.0);
    const TVector3D t          = sensorPose.translation();

    ImuTransformer tf;
    tf.parameters.ang_acc_lpf_alpha = 1.0;  // raw finite difference, no lag
    tf.parameters.ang_vel_lpf_alpha = 1.0;

    const double    dt      = 0.01;
    const double    alpha_z = 4.0;  // rad/s^2
    const TVector3D acc     = {0.0, 0.0, 9.81};

    // wz(t) = alpha_z * (i*dt), so finite-difference angular accel = alpha_z.
    mrpt::obs::CObservationIMU out;
    TVector3D                  w_last = {0, 0, 0};
    for (int i = 0; i < 6; ++i)
    {
        w_last = {0.0, 0.0, alpha_z * (i * dt)};
        out    = tf.process(*make_imu(100.0 + dt * i, acc, w_last, sensorPose));
    }

    // Expected = acc - alpha x t - w x (w x t), with alpha=(0,0,alpha_z), w=w_last.
    const TVector3D alpha = {0.0, 0.0, alpha_z};
    const TVector3D expected =
        acc - mrpt::math::crossProduct3D(alpha, t) -
        mrpt::math::crossProduct3D(w_last, mrpt::math::crossProduct3D(w_last, t));

    check_near(get_acc(out), expected, 1e-9, "tangential corrected accel");

    std::cout << "✅ test_tangential_lever_arm passed!" << std::endl;
}

// 4) First sample introduces no angular-acceleration spike, and output is finite.
void test_first_sample_no_spike()
{
    const auto sensorPose = mrpt::poses::CPose3D::FromTranslation(1.0, 2.0, 3.0);

    ImuTransformer  tf;
    const TVector3D w   = {0.5, -0.3, 0.2};
    const TVector3D acc = {0.1, 0.2, 9.8};

    const auto out  = tf.process(*make_imu(100.0, acc, w, sensorPose));
    const auto aout = get_acc(out);

    ASSERT_(std::isfinite(aout.x) && std::isfinite(aout.y) && std::isfinite(aout.z));

    // On the first sample angular accel is forced to zero, so only the
    // centripetal term (using w bootstrapped to the current reading) applies.
    const TVector3D t = sensorPose.translation();
    const TVector3D expected =
        acc - mrpt::math::crossProduct3D(w, mrpt::math::crossProduct3D(w, t));
    check_near(aout, expected, 1e-9, "first-sample accel (no tangential spike)");

    std::cout << "✅ test_first_sample_no_spike passed!" << std::endl;
}

// 5) Large/negative timestamp gaps fall back to the default rate without NaNs.
void test_dt_fallback()
{
    const auto sensorPose = mrpt::poses::CPose3D::FromTranslation(0.2, 0.0, 0.0);

    ImuTransformer  tf;
    const TVector3D w   = {0.0, 0.0, 1.0};
    const TVector3D acc = {0.0, 0.0, 9.81};

    tf.process(*make_imu(100.0, acc, w, sensorPose));
    // Gap > 1 second (data flow resumed): must not throw nor produce NaN.
    const auto out  = tf.process(*make_imu(150.0, acc, w, sensorPose));
    const auto aout = get_acc(out);
    ASSERT_(std::isfinite(aout.x) && std::isfinite(aout.y) && std::isfinite(aout.z));

    std::cout << "✅ test_dt_fallback passed!" << std::endl;
}

// 6) Near-duplicate timestamps (a few microseconds apart, as emitted by some
//    IMU drivers) must not amplify a tiny angular-velocity difference into a
//    huge spurious lever-arm spike via division by a near-zero dt.
void test_near_duplicate_timestamp_no_spike()
{
    const auto sensorPose = mrpt::poses::CPose3D::FromTranslation(0.5, 0.0, 0.0);

    ImuTransformer  tf;
    const TVector3D acc = {0.0, 0.0, 9.81};

    tf.process(*make_imu(100.0, acc, {0.0, 0.0, 0.0}, sensorPose));
    // A microsecond-scale gap with a tiny gyro change, as seen from near-duplicate
    // driver messages:
    const auto out  = tf.process(*make_imu(100.0 + 5e-6, acc, {0.0, 0.0, 0.005}, sensorPose));
    const auto aout = get_acc(out);

    ASSERT_(std::isfinite(aout.x) && std::isfinite(aout.y) && std::isfinite(aout.z));
    check_near(aout, acc, 0.5, "near-duplicate-timestamp accel (no spike)");

    std::cout << "✅ test_near_duplicate_timestamp_no_spike passed!" << std::endl;
}

}  // namespace

int main([[maybe_unused]] int argc, [[maybe_unused]] char** argv)
{
    try
    {
        test_rotation_only();
        test_centripetal_lever_arm();
        test_tangential_lever_arm();
        test_first_sample_no_spike();
        test_dt_fallback();
        test_near_duplicate_timestamp_no_spike();
        std::cout << "All ImuTransformer tests passed." << std::endl;
        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Test failed: " << e.what() << std::endl;
        return 1;
    }
}
