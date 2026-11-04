/*               _
 _ __ ___   ___ | | __ _
| '_ ` _ \ / _ \| |/ _` | Modular Optimization framework for
| | | | | | (_) | | (_| | Localization and mApping (MOLA)
|_| |_| |_|\___/|_|\__,_| https://github.com/MOLAorg/mola

 Copyright (C) 2018-2026 Jose Luis Blanco, University of Almeria,
                         and individual contributors.
 SPDX-License-Identifier: GPL-3.0
 See LICENSE for full license information.
 Closed-source licenses available upon request, for this odometry package
 alone or in combination with the complete SLAM system.
*/

/**
 * Tests for ImuBiasTracker, with simulated IMU streams whose noise and vibration
 * levels are those measured on a skid-steer ground robot (Clearpath Husky):
 *  - at standstill with the motors on: ~0.0035 rad/s per-sample gyro noise
 *    (100 Hz), ~0.03 m/s^2 accelerometer noise;
 *  - driving on paving: ~0.05 rad/s gyro and ~2 m/s^2 accelerometer vibration;
 *  - biases like those of a LiDAR-integrated MEMS IMU: ~-0.2 deg/s on z and
 *    ~-0.3 deg/s on y, plus a few cm/s^2 of accelerometer bias.
 */

#include <mola_imu_preintegration/ImuBiasTracker.h>
#include <mrpt/core/bits_math.h>
#include <mrpt/core/exceptions.h>
#include <mrpt/obs/CObservationIMU.h>
#include <mrpt/poses/CPose3D.h>
#include <mrpt/random/RandomGenerators.h>

#include <Eigen/Dense>
#include <cmath>
#include <iostream>

using namespace mola::imu;
using mrpt::math::TVector3D;

namespace
{
constexpr double RATE = 100.0;  // [Hz]
constexpr double G    = 9.81;

/// A simulated IMU, rigidly mounted with some tilt on a ground vehicle.
struct SimImu
{
    std::string          label = "imu";
    mrpt::poses::CPose3D mounting =
        mrpt::poses::CPose3D::FromYawPitchRoll(0, mrpt::DEG2RAD(-3.0), mrpt::DEG2RAD(2.0));
    TVector3D bias_gyro{0.0002, -0.0051, -0.0036};  // ~0.3 deg/s on y, ~0.2 deg/s on z
    TVector3D bias_acc{0.05, -0.03, 0.08};
    double    gyro_noise = 0.0035;  // per sample [rad/s]
    double    acc_noise  = 0.03;  // per sample [m/s^2]
    double    t          = 1.8e9;  // [s] UNIX-like time
    mrpt::random::CRandomGenerator rng;

    SimImu() { rng.randomize(1); }

    /// One reading; vehicle angular rate wz [rad/s] about its vertical axis and an
    /// extra vibration level (0 = standing still).
    mrpt::obs::CObservationIMU sample(double wz, double gyroVib = 0, double accVib = 0)
    {
        using namespace mrpt::obs;
        t += 1.0 / RATE;
        // Vehicle rate and gravity, rotated to the sensor frame:
        const Eigen::Matrix3d Rt  = mounting.getRotationMatrix().asEigen().transpose();
        const Eigen::Vector3d wV  = Rt * Eigen::Vector3d(0, 0, wz);
        const Eigen::Vector3d gUp = Rt * Eigen::Vector3d(0, 0, G);
        CObservationIMU       o;
        o.sensorLabel   = label;
        o.timestamp     = mrpt::Clock::fromDouble(t);
        const double wn = gyro_noise + gyroVib;
        const double an = acc_noise + accVib;
        o.set(IMU_WX, wV.x() + bias_gyro.x + rng.drawGaussian1D(0, wn));
        o.set(IMU_WY, wV.y() + bias_gyro.y + rng.drawGaussian1D(0, wn));
        o.set(IMU_WZ, wV.z() + bias_gyro.z + rng.drawGaussian1D(0, wn));
        o.set(IMU_X_ACC, gUp.x() + bias_acc.x + rng.drawGaussian1D(0, an));
        o.set(IMU_Y_ACC, gUp.y() + bias_acc.y + rng.drawGaussian1D(0, an));
        o.set(IMU_Z_ACC, gUp.z() + bias_acc.z + rng.drawGaussian1D(0, an));
        return o;
    }
};

/// Feeds `seconds` of readings; returns how many standstill updates happened.
template <typename F>
std::size_t feed(ImuBiasTracker& tr, SimImu& imu, double seconds, F&& gen)
{
    std::size_t n = 0;
    for (int i = 0; i < static_cast<int>(seconds * RATE); i++)
    {
        if (tr.add(gen(imu)))
        {
            n++;
        }
    }
    return n;
}

auto   standing = [](SimImu& s) { return s.sample(0); };
auto   driving  = [](SimImu& s) { return s.sample(0.3 * std::sin(s.t * 0.2), 0.05, 2.0); };
double norm3(const TVector3D& v) { return v.norm(); }

// ---------------------------------------------------------------------------
void test_standstill_converges()
{
    ImuBiasTracker tr;
    SimImu         imu;
    const auto     n = feed(tr, imu, 10.0, standing);
    ASSERT_GE_(n, 8U);

    const auto e = tr.estimate(imu.label);
    ASSERT_(e.has_value());
    ASSERT_(e->has_measurements());
    // 10 s of 0.0035 rad/s noise at 100 Hz: the mean has a sigma of ~1e-4 rad/s.
    const double errG = norm3(e->bias_gyro - imu.bias_gyro);
    ASSERT_LT_(errG, 4e-4);  // ~3.5 sigma of the norm of the mean error
    ASSERT_LT_(std::sqrt(e->cov_gyro(2, 2)), 3e-4);

    // Only the accelerometer bias component along gravity is observable at rest:
    const Eigen::Matrix3d Rt = imu.mounting.getRotationMatrix().asEigen().transpose();
    const Eigen::Vector3d u  = Rt * Eigen::Vector3d(0, 0, 1);
    const double trueU = u.dot(Eigen::Vector3d(imu.bias_acc.x, imu.bias_acc.y, imu.bias_acc.z));
    const double estU  = u.dot(Eigen::Vector3d(e->bias_acc.x, e->bias_acc.y, e->bias_acc.z));
    ASSERT_LT_(std::abs(estU - trueU), 0.01);
    const Eigen::Matrix3d Pa         = e->cov_acc.asEigen();
    const double          sigmaAlong = std::sqrt(u.dot(Pa * u));
    const double          sigmaX     = std::sqrt(e->cov_acc(0, 0));
    ASSERT_LT_(sigmaAlong, 0.03);
    ASSERT_GT_(sigmaX, 0.25);  // across gravity: still the prior

    std::cout << "standstill: " << n << " updates, gyro bias error " << errG
              << " rad/s, accel bias along gravity error " << std::abs(estU - trueU) << " m/s^2\n";
}

// ---------------------------------------------------------------------------
void test_driving_does_not_update()
{
    ImuBiasTracker tr;
    SimImu         imu;
    feed(tr, imu, 5.0, standing);
    const auto before = *tr.estimate(imu.label);

    const auto n = feed(tr, imu, 60.0, driving);
    ASSERT_EQUAL_(n, 0U);
    const auto after = *tr.estimate(imu.label);
    ASSERT_EQUAL_(after.standstill_updates, before.standstill_updates);
    ASSERT_LT_(norm3(after.bias_gyro - before.bias_gyro), 1e-12);
    std::cout << "driving: no updates in 60 s of driving\n";
}

// ---------------------------------------------------------------------------
void test_stop_and_go_tracks_drift()
{
    ImuBiasTracker tr;
    SimImu         imu;
    // 10 minutes: drive 55 s, stop 5 s, while the z gyro bias drifts by 0.002 rad/s
    // (e.g. warming up).
    const TVector3D b0       = imu.bias_gyro;
    const double    t0       = imu.t;
    std::size_t     stops    = 0;
    double          lastCovZ = 0;
    for (int cycle = 0; cycle < 10; cycle++)
    {
        auto drift = [&](SimImu& s, bool moving)
        {
            s.bias_gyro.z = b0.z + 0.002 * (s.t - t0) / 600.0;
            return moving ? s.sample(0.3 * std::sin(s.t * 0.2), 0.05, 2.0) : s.sample(0);
        };
        feed(tr, imu, 55.0, [&](SimImu& s) { return drift(s, true); });
        // The uncertainty grows while driving:
        const auto e = *tr.estimate(imu.label, imu.t);
        if (cycle > 0)
        {
            ASSERT_GT_(e.cov_gyro(2, 2), lastCovZ);
        }
        stops += feed(tr, imu, 5.0, [&](SimImu& s) { return drift(s, false); }) > 0 ? 1 : 0;
        lastCovZ = tr.estimate(imu.label)->cov_gyro(2, 2);
    }
    ASSERT_EQUAL_(stops, 10U);
    const auto   e   = *tr.estimate(imu.label);
    const double err = std::abs(e.bias_gyro.z - imu.bias_gyro.z);
    ASSERT_LT_(err, 3e-4);
    std::cout << "stop and go: z bias tracked within " << err
              << " rad/s after a 0.002 rad/s drift\n";
}

// ---------------------------------------------------------------------------
void test_slow_rotation_and_hints()
{
    // A smooth 0.01 rad/s turn without vibration: from the IMU alone it looks like
    // a bias (documented limitation).
    auto slowTurn = [](SimImu& s) { return s.sample(0.01); };
    {
        ImuBiasTracker tr;
        SimImu         imu;
        imu.gyro_noise = 0.001;
        ASSERT_GT_(feed(tr, imu, 3.0, slowTurn), 0U);
    }
    // "Moving" hints veto it:
    {
        ImuBiasTracker tr;
        SimImu         imu;
        imu.gyro_noise = 0.001;
        std::size_t n  = 0;
        for (int i = 0; i < 300; i++)
        {
            const auto o = slowTurn(imu);
            tr.add_motion_hint(imu.t, false);
            n += tr.add(o) ? 1 : 0;
        }
        ASSERT_EQUAL_(n, 0U);
    }
    // Requiring a static hint: none given, no updates even standing still...
    {
        ImuBiasTracker tr;
        tr.parameters.require_static_hint = true;
        SimImu imu;
        ASSERT_EQUAL_(feed(tr, imu, 3.0, standing), 0U);
        // ...and with "static" hints, it works:
        std::size_t n = 0;
        for (int i = 0; i < 300; i++)
        {
            const auto o = standing(imu);
            if (i % 10 == 0)
            {
                tr.add_motion_hint(imu.t, true);
            }
            n += tr.add(o) ? 1 : 0;
        }
        ASSERT_GT_(n, 0U);
    }
    std::cout << "slow rotation: taken as bias without hints, rejected with them\n";
}

// ---------------------------------------------------------------------------
void test_external_measurement_matches_closed_form()
{
    ImuBiasTracker tr;
    tr.parameters.gyro_bias_prior_sigma = 0.02;
    tr.parameters.gyro_bias_random_walk = 0;
    const double                t       = 100.0;
    mrpt::math::CMatrixDouble33 R1;
    R1.setDiagonal(1e-6);
    mrpt::math::CMatrixDouble33 R2;
    R2.setDiagonal(4e-6);
    tr.add_measurement("imu", t, std::make_pair(TVector3D(0.004, 0, 0), R1), std::nullopt);
    tr.add_measurement("imu", t, std::make_pair(TVector3D(0.001, 0, 0), R2), std::nullopt);
    const auto e = *tr.estimate("imu");
    // Inverse-variance fusion of the zero prior and both measurements:
    const double i0 = 1.0 / mrpt::square(0.02);
    const double i1 = 1.0 / 1e-6;
    const double i2 = 1.0 / 4e-6;
    const double x  = (i1 * 0.004 + i2 * 0.001) / (i0 + i1 + i2);
    ASSERT_NEAR_(e.bias_gyro.x, x, 1e-12);
    ASSERT_NEAR_(e.cov_gyro(0, 0), 1.0 / (i0 + i1 + i2), 1e-15);
    ASSERT_EQUAL_(e.external_updates, 2U);
    std::cout << "external measurements: Kalman fusion equals the closed form\n";
}

// ---------------------------------------------------------------------------
void test_gaps_restart_the_window()
{
    ImuBiasTracker tr;
    SimImu         imu;
    std::size_t    n = 0;
    // 0.6 s of readings, a 0.5 s gap, 0.6 s more: no 1 s window without a gap.
    n += feed(tr, imu, 0.6, standing);
    imu.t += 0.5;
    n += feed(tr, imu, 0.6, standing);
    ASSERT_EQUAL_(n, 0U);
    n += feed(tr, imu, 0.6, standing);
    ASSERT_EQUAL_(n, 1U);
    std::cout << "gaps: a window never spans a gap\n";
}

// ---------------------------------------------------------------------------
void test_independent_imus()
{
    ImuBiasTracker tr;
    SimImu         a;
    SimImu         b;
    b.label     = "imu2";
    b.bias_gyro = {-0.001, 0.002, 0.0005};
    for (int i = 0; i < 1000; i++)
    {
        tr.add(a.sample(0));
        tr.add(b.sample(0));
    }
    ASSERT_EQUAL_(tr.sensor_labels().size(), 2U);
    ASSERT_LT_(norm3(tr.estimate("imu")->bias_gyro - a.bias_gyro), 3e-4);
    ASSERT_LT_(norm3(tr.estimate("imu2")->bias_gyro - b.bias_gyro), 3e-4);
    ASSERT_(!tr.estimate("unknown").has_value());
    std::cout << "two IMUs: independent estimates\n";
}

// ---------------------------------------------------------------------------
void test_yaml_roundtrip()
{
    ImuBiasTracker::Parameters p;
    p.standstill_window      = 2.5;
    p.require_static_hint    = true;
    p.gyro_bias_random_walk  = 1e-4;
    mrpt::containers::yaml y = mrpt::containers::yaml::Map();
    p.save_to(y);
    ImuBiasTracker::Parameters q;
    q.load_from(y);
    ASSERT_EQUAL_(q.standstill_window, 2.5);
    ASSERT_(q.require_static_hint);
    ASSERT_EQUAL_(q.gyro_bias_random_walk, 1e-4);
    std::cout << "parameters: YAML round trip\n";
}

}  // namespace

int main()
{
    try
    {
        test_standstill_converges();
        test_driving_does_not_update();
        test_stop_and_go_tracks_drift();
        test_slow_rotation_and_hints();
        test_external_measurement_matches_closed_form();
        test_gaps_restart_the_window();
        test_independent_imus();
        test_yaml_roundtrip();
        std::cout << "All ImuBiasTracker tests passed.\n";
        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << mrpt::exception_to_str(e) << "\n";
        return 1;
    }
}
