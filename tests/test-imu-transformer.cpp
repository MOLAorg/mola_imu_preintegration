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

// 5) A large timestamp gap (data flow stopped and resumed) restarts the filters:
//    the angular-velocity change across the gap must not become a tangential spike.
void test_dt_gap_restarts_filters()
{
    const auto      sensorPose = mrpt::poses::CPose3D::FromTranslation(0.2, 0.0, 0.0);
    const TVector3D t          = sensorPose.translation();

    ImuTransformer  tf;
    const TVector3D w   = {0.0, 0.0, 1.0};
    const TVector3D acc = {0.0, 0.0, 9.81};

    tf.process(*make_imu(100.0, acc, {0.0, 0.0, 0.0}, sensorPose));
    // Gap > 1 second, with a different angular velocity afterwards:
    const auto out  = tf.process(*make_imu(150.0, acc, w, sensorPose));
    const auto aout = get_acc(out);
    ASSERT_(std::isfinite(aout.x) && std::isfinite(aout.y) && std::isfinite(aout.z));

    // Only the centripetal term, with w bootstrapped to the current reading:
    const TVector3D expected =
        acc - mrpt::math::crossProduct3D(w, mrpt::math::crossProduct3D(w, t));
    check_near(aout, expected, 1e-9, "accel after a timestamp gap");

    // Time going backwards must behave the same way:
    const auto out2 = tf.process(*make_imu(140.0, acc, w, sensorPose));
    check_near(get_acc(out2), expected, 1e-9, "accel after a backwards timestamp");

    std::cout << "✅ test_dt_gap_restarts_filters passed!" << std::endl;
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

// 7) High-rate IMUs (sample period well below 1 ms) must use their real dt in
//    the finite difference, so the tangential term matches the true angular
//    acceleration.
void test_high_rate_tangential_lever_arm()
{
    const auto      sensorPose = mrpt::poses::CPose3D::FromTranslation(0.5, 0.0, 0.0);
    const TVector3D t          = sensorPose.translation();

    ImuTransformer tf;
    tf.parameters.ang_acc_lpf_alpha = 1.0;
    tf.parameters.ang_vel_lpf_alpha = 1.0;

    const double    dt      = 1.0 / 2560.0;
    const double    alpha_z = 4.0;  // rad/s^2
    const TVector3D acc     = {0.0, 0.0, 9.81};

    mrpt::obs::CObservationIMU out;
    TVector3D                  w_last = {0, 0, 0};
    for (int i = 0; i < 20; ++i)
    {
        w_last = {0.0, 0.0, alpha_z * (i * dt)};
        out    = tf.process(*make_imu(100.0 + dt * i, acc, w_last, sensorPose));
    }

    const TVector3D alpha = {0.0, 0.0, alpha_z};
    const TVector3D expected =
        acc - mrpt::math::crossProduct3D(alpha, t) -
        mrpt::math::crossProduct3D(w_last, mrpt::math::crossProduct3D(w_last, t));

    // The tolerance covers timestamp rounding in mrpt::Clock (older MRPT
    // releases resolve only ~1.6 us through toDouble()), and is still far
    // below the error of differencing over a wrong dt.
    check_near(get_acc(out), expected, 0.05, "high-rate tangential corrected accel");

    std::cout << "✅ test_high_rate_tangential_lever_arm passed!" << std::endl;
}

// Angular acceleration estimated after a step of `alpha_z` starting at t=0, at
// time `t_end`, sampling at `rate_hz`. Read from the tangential lever-arm term.
double estimated_ang_acc_after_step(
    ImuTransformer& tf, double rate_hz, double alpha_z, double t_end)
{
    const auto      sensorPose = mrpt::poses::CPose3D::FromTranslation(0.5, 0.0, 0.0);
    const TVector3D acc        = {0.0, 0.0, 9.81};
    const double    dt         = 1.0 / rate_hz;
    const int       n          = static_cast<int>(std::lround(t_end / dt));

    mrpt::obs::CObservationIMU out;
    for (int i = 0; i <= n; ++i)
    {
        const TVector3D w = {0.0, 0.0, alpha_z * (i * dt)};
        out               = tf.process(*make_imu(100.0 + dt * i, acc, w, sensorPose));
    }
    // With w along z and the lever arm along x, the tangential term is
    // (0, alpha*0.5, 0) and the centripetal one is along x only:
    return -2.0 * (get_acc(out).y - acc.y);
}

// 8) The default low-pass filters are specified by a cut-off frequency, so their
//    response does not depend on the IMU rate.
void test_cutoff_is_rate_independent()
{
    const double alpha_z = 4.0;
    const double t_end   = 0.05;

    ImuTransformer tf_100;
    ImuTransformer tf_400;
    const double   a100 = estimated_ang_acc_after_step(tf_100, 100.0, alpha_z, t_end);
    const double   a400 = estimated_ang_acc_after_step(tf_400, 400.0, alpha_z, t_end);

    const double fc       = tf_100.parameters.ang_acc_lpf_cutoff_hz;
    const double expected = alpha_z * (1.0 - std::exp(-2.0 * M_PI * fc * t_end));

    ASSERT_NEAR_(a100, expected, 1e-6);
    ASSERT_NEAR_(a400, expected, 1e-6);

    std::cout << "✅ test_cutoff_is_rate_independent passed!" << std::endl;
}

// 9) A legacy per-sample alpha, if set, overrides the cut-off frequency.
void test_legacy_alpha_overrides_cutoff()
{
    const double alpha_z = 4.0;
    const double t_end   = 0.05;
    const double rate    = 100.0;
    const double a       = 0.02;

    ImuTransformer tf;
    tf.parameters.ang_acc_lpf_alpha = a;
    const double est                = estimated_ang_acc_after_step(tf, rate, alpha_z, t_end);

    const int    n        = static_cast<int>(std::lround(t_end * rate));
    const double expected = alpha_z * (1.0 - std::pow(1.0 - a, n));
    ASSERT_NEAR_(est, expected, 1e-6);

    std::cout << "✅ test_legacy_alpha_overrides_cutoff passed!" << std::endl;
}

// Filtered angular velocity of the centripetal lever-arm term, after a step of
// `w_z` right after the first reading, at time `t_end`, sampling at `rate_hz`.
double filtered_ang_vel_after_step(ImuTransformer& tf, double rate_hz, double w_z, double t_end)
{
    const auto      sensorPose = mrpt::poses::CPose3D::FromTranslation(0.5, 0.0, 0.0);
    const TVector3D acc        = {0.0, 0.0, 9.81};
    const TVector3D w_zero     = {0.0, 0.0, 0.0};
    const TVector3D w_step     = {0.0, 0.0, w_z};
    const double    dt         = 1.0 / rate_hz;
    const int       n          = static_cast<int>(std::lround(t_end / dt));

    auto out = tf.process(*make_imu(100.0, acc, w_zero, sensorPose));
    for (int i = 1; i <= n; ++i)
    {
        out = tf.process(*make_imu(100.0 + dt * i, acc, w_step, sensorPose));
    }
    // With w along z and the lever arm along x, the centripetal term adds
    // 0.5*w^2 along x, and the tangential one is along y only:
    return std::sqrt(2.0 * (get_acc(out).x - acc.x));
}

// 10) The centripetal term has its own cut-off frequency, which <= 0 disables,
//     and which a legacy per-sample alpha overrides.
void test_ang_vel_cutoff()
{
    const double w_z   = 2.0;
    const double t_end = 0.05;
    const double rate  = 200.0;
    const int    n     = static_cast<int>(std::lround(t_end * rate));

    // Different cut-offs, so that mixing up the two filters fails:
    ImuTransformer tf;
    tf.parameters.ang_vel_lpf_cutoff_hz = 5.0;
    tf.parameters.ang_acc_lpf_cutoff_hz = 40.0;
    ASSERT_NEAR_(
        filtered_ang_vel_after_step(tf, rate, w_z, t_end),
        w_z * (1.0 - std::exp(-2.0 * M_PI * 5.0 * t_end)), 1e-6);

    ImuTransformer tf_unfiltered;
    tf_unfiltered.parameters.ang_vel_lpf_cutoff_hz = 0;
    ASSERT_NEAR_(filtered_ang_vel_after_step(tf_unfiltered, rate, w_z, t_end), w_z, 1e-6);

    ImuTransformer tf_legacy;
    tf_legacy.parameters.ang_vel_lpf_alpha = 0.3;
    ASSERT_NEAR_(
        filtered_ang_vel_after_step(tf_legacy, rate, w_z, t_end),
        w_z * (1.0 - std::pow(1.0 - 0.3, n)), 1e-6);

    std::cout << "✅ test_ang_vel_cutoff passed!" << std::endl;
}

// 11) Accelerometer units: readings in units of g are output in m/s², those in
//     m/s² unchanged, from the first reading on. The decision ignores readings
//     without information (zeros, NaN), is kept after its first 0.5 s, and an
//     explicit accel_scale is applied as is.
void test_accel_units()
{
    const auto       sensorPose = mrpt::poses::CPose3D::Identity();
    const TVector3D  w          = {0.0, 0.0, 0.0};
    constexpr double DT         = 0.005;  // 200 Hz

    // Feeds `n` readings of acc_z from reading index `i0`, and checks that every
    // output equals `expected`:
    const auto feed =
        [&](ImuTransformer& tf, int i0, int n, double acc_z, double expected, const char* what)
    {
        for (int i = i0; i < i0 + n; ++i)
        {
            const auto out =
                tf.process(*make_imu(100.0 + DT * i, {0.0, 0.0, acc_z}, w, sensorPose));
            if (std::abs(get_acc(out).z - expected) > 1e-9)
            {
                THROW_EXCEPTION_FMT(
                    "%s: reading %d gave %f, expected %f", what, i, get_acc(out).z, expected);
            }
        }
    };

    // Units of g, then frozen: a later reading of 9.81 (in g) is still scaled.
    ImuTransformer tf_g;
    feed(tf_g, 0, 120, 1.0, 9.80665, "units of g");
    feed(tf_g, 120, 10, 9.81, 9.81 * 9.80665, "units of g, frozen");

    // m/s², then frozen: a later reading of 1.0 is not scaled.
    ImuTransformer tf_si;
    feed(tf_si, 0, 120, 9.81, 9.81, "m/s^2");
    feed(tf_si, 120, 10, 1.0, 1.0, "m/s^2, frozen");

    // Zeros while a driver starts do not count:
    ImuTransformer tf_zeros;
    feed(tf_zeros, 0, 80, 0.0, 0.0, "leading zeros");
    feed(tf_zeros, 80, 120, 9.81, 9.81, "m/s^2 after zeros");

    // Nor does a non-finite reading:
    ImuTransformer tf_nan;
    tf_nan.process(*make_imu(100.0, {0.0, 0.0, std::nan("")}, w, sensorPose));
    feed(tf_nan, 1, 120, 1.0, 9.80665, "units of g after NaN");

    // An explicit scale:
    ImuTransformer tf_fixed;
    tf_fixed.parameters.accel_scale = 1.0;
    feed(tf_fixed, 0, 10, 1.0, 1.0, "explicit scale");

    // A negative one is rejected:
    ImuTransformer tf_neg;
    tf_neg.parameters.accel_scale = -1.0;
    bool thrown                   = false;
    try
    {
        tf_neg.process(*make_imu(100.0, {0.0, 0.0, 9.81}, w, sensorPose));
    }
    catch (const std::exception&)
    {
        thrown = true;
    }
    ASSERT_(thrown);

    std::cout << "✅ test_accel_units passed!" << std::endl;
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
        test_dt_gap_restarts_filters();
        test_near_duplicate_timestamp_no_spike();
        test_high_rate_tangential_lever_arm();
        test_cutoff_is_rate_independent();
        test_legacy_alpha_overrides_cutoff();
        test_ang_vel_cutoff();
        test_accel_units();
        std::cout << "All ImuTransformer tests passed." << std::endl;
        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Test failed: " << e.what() << std::endl;
        return 1;
    }
}
