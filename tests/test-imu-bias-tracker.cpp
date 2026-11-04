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
 * levels are typical of a skid-steer ground robot:
 *  - at standstill with the motors on: ~0.0035 rad/s per-sample gyro noise
 *    (100 Hz), ~0.03 m/s^2 accelerometer noise;
 *  - driving on paving: ~0.05 rad/s gyro and ~2 m/s^2 accelerometer vibration;
 *  - biases like those of a MEMS IMU integrated in a LiDAR: ~-0.2 deg/s on z and
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
#include <deque>
#include <iostream>
#include <limits>

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
    bool      with_acc   = true;
    mrpt::random::CRandomGenerator rng;

    explicit SimImu(uint32_t seed = 1) { rng.randomize(seed); }

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
        if (with_acc)
        {
            o.set(IMU_X_ACC, gUp.x() + bias_acc.x + rng.drawGaussian1D(0, an));
            o.set(IMU_Y_ACC, gUp.y() + bias_acc.y + rng.drawGaussian1D(0, an));
            o.set(IMU_Z_ACC, gUp.z() + bias_acc.z + rng.drawGaussian1D(0, an));
        }
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

Eigen::Vector3d gravity_dir(const SimImu& imu)
{
    return imu.mounting.getRotationMatrix().asEigen().transpose() * Eigen::Vector3d(0, 0, 1);
}

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
    const double    q        = mrpt::square(tr.parameters.bias_gyro_random_walk);
    for (int cycle = 0; cycle < 10; cycle++)
    {
        auto drift = [&](SimImu& s, bool moving)
        {
            s.bias_gyro.z = b0.z + 0.002 * (s.t - t0) / 600.0;
            return moving ? s.sample(0.3 * std::sin(s.t * 0.2), 0.05, 2.0) : s.sample(0);
        };
        feed(tr, imu, 55.0, [&](SimImu& s) { return drift(s, true); });
        // While driving, the uncertainty grows as the random walk:
        const auto e = *tr.estimate(imu.label);
        if (cycle > 0)
        {
            ASSERT_NEAR_(e.cov_gyro(2, 2) - lastCovZ, q * 55.0, q * 0.05);
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
void test_large_bias_accepted()
{
    // A consumer-grade gyro, ~1 deg/s off: above standstill_max_gyro_mean_offset,
    // accepted while the estimate is still uncertain.
    ImuBiasTracker tr;
    SimImu         imu;
    imu.bias_gyro = {0.012, -0.015, 0.018};
    ASSERT_GE_(feed(tr, imu, 10.0, standing), 8U);
    ASSERT_LT_(norm3(tr.estimate(imu.label)->bias_gyro - imu.bias_gyro), 4e-4);
    std::cout << "large bias: 0.026 rad/s accepted at first\n";
}

// ---------------------------------------------------------------------------
void test_bias_changes_are_reacquired()
{
    // Once converged, a step of the bias (e.g. after a power cycle) looks like a
    // rotation, and is rejected at first...
    {
        ImuBiasTracker tr;
        SimImu         imu;
        feed(tr, imu, 10.0, standing);
        imu.bias_gyro.z += 0.02;
        imu.t += 1.0;  // a new window
        ASSERT_EQUAL_(feed(tr, imu, 25.0, standing), 0U);
        // ...until it has been consistent for gyro_reacquire_time:
        ASSERT_GT_(feed(tr, imu, 15.0, standing), 0U);
        ASSERT_LT_(norm3(tr.estimate(imu.label)->bias_gyro - imu.bias_gyro), 1e-3);
    }
    // With "static" hints, the step is accepted after two windows:
    {
        ImuBiasTracker tr;
        SimImu         imu;
        tr.parameters.max_hint_latency = 0;
        feed(tr, imu, 10.0, standing);
        imu.bias_gyro.z += 0.02;
        imu.t += 1.0;
        std::size_t n = 0;
        for (int i = 0; i < 500; i++)
        {
            const auto o = standing(imu);
            if (i % 10 == 0)
            {
                tr.add_motion_hint(imu.t, true);
            }
            n += tr.add(o);
        }
        ASSERT_GE_(n, 2U);
        ASSERT_LT_(norm3(tr.estimate(imu.label)->bias_gyro - imu.bias_gyro), 1e-3);
    }
    // After a long time without updates, the gate has opened by the random walk:
    {
        ImuBiasTracker tr;
        SimImu         imu;
        feed(tr, imu, 10.0, standing);
        imu.t += 20 * 3600.0;  // e.g. a day of driving without stops
        imu.bias_gyro.z += 0.015;
        ASSERT_GE_(feed(tr, imu, 3.0, standing), 2U);
        ASSERT_LT_(norm3(tr.estimate(imu.label)->bias_gyro - imu.bias_gyro), 1e-3);
    }
    std::cout << "bias changes: re-acquired after a while, at once with static hints\n";
}

// ---------------------------------------------------------------------------
void test_reacquisition_robustness()
{
    // Stop and go: a bias step after an IMU restart is picked up over several
    // 10 s stops, the re-acquisition time adding up across the drives.
    {
        ImuBiasTracker tr;
        SimImu         imu;
        feed(tr, imu, 10.0, standing);
        imu.t += 2.0;  // the IMU restarts...
        imu.bias_gyro.z += 0.02;  // ...with another bias
        for (int cycle = 0; cycle < 6; cycle++)
        {
            feed(tr, imu, 30.0, driving);
            feed(tr, imu, 10.0, standing);
        }
        ASSERT_LT_(norm3(tr.estimate(imu.label)->bias_gyro - imu.bias_gyro), 1e-3);
    }
    // Mid-size steps (within the absolute gate, beyond the estimate sigma) are
    // followed within a few windows:
    {
        ImuBiasTracker tr;
        SimImu         imu;
        feed(tr, imu, 10.0, standing);
        imu.bias_gyro.z += 0.008;
        imu.t += 1.0;
        feed(tr, imu, 4.0, standing);
        ASSERT_LT_(norm3(tr.estimate(imu.label)->bias_gyro - imu.bias_gyro), 1e-3);
    }
    // Wrong "static" hints during a fast smooth turn: beyond any plausible bias,
    // never taken as one.
    {
        ImuBiasTracker tr;
        SimImu         imu;
        tr.parameters.max_hint_latency = 0;
        feed(tr, imu, 10.0, standing);
        const auto before = *tr.estimate(imu.label);
        for (int i = 0; i < 300; i++)
        {
            tr.add_motion_hint(imu.t, true);
            tr.add(imu.sample(0.3));
        }
        const auto after = *tr.estimate(imu.label);
        ASSERT_EQUAL_(after.standstill_updates, before.standstill_updates);
    }
    // The body settling on its suspension after a stop (a short pitch rate, with
    // "static" wheel odometry hints): one window is not enough to reset the bias.
    {
        ImuBiasTracker tr;
        SimImu         imu;
        tr.parameters.max_hint_latency = 0;
        feed(tr, imu, 10.0, standing);
        imu.t += 1.0;
        for (int i = 0; i < 115; i++)
        {
            tr.add_motion_hint(imu.t, true);
            auto o = imu.sample(0);
            o.set(mrpt::obs::IMU_WY, o.get(mrpt::obs::IMU_WY) + 0.015);
            tr.add(o);
        }
        ASSERT_LT_(norm3(tr.estimate(imu.label)->bias_gyro - imu.bias_gyro), 1e-3);
    }
    std::cout << "re-acquisition: over interruptions, mid-size steps, bounded by the prior\n";
}

// ---------------------------------------------------------------------------
void test_hint_sources_and_clock_jumps()
{
    // Two hint sources, one 1.5 s late: its stale hints do not erase the on-time
    // "moving" ones during a smooth turn.
    {
        ImuBiasTracker tr;
        SimImu         imu;
        imu.gyro_noise = 0.001;
        feed(tr, imu, 10.0, standing);
        std::deque<std::pair<double, bool>> late;
        std::size_t                         n = 0;
        for (int i = 0; i < 1000; i++)
        {
            const bool turning = i >= 200;
            const auto o       = imu.sample(turning ? 0.03 : 0.0);
            if (i % 10 == 0)
            {
                n += tr.add_motion_hint(imu.t, !turning).size();
                late.emplace_back(imu.t, !turning);
            }
            n += tr.add(o);
            while (!late.empty() && imu.t - late.front().first >= 1.5)
            {
                n += tr.add_motion_hint(late.front().first, late.front().second).size();
                late.pop_front();
            }
        }
        const auto e = *tr.estimate(imu.label);
        ASSERT_LT_(norm3(e.bias_gyro - imu.bias_gyro), 1e-3);
    }
    // A dataset replayed from its start: both clocks jump back by minutes, and
    // the tracker keeps working with the new hints.
    {
        ImuBiasTracker tr;
        SimImu         imu;
        tr.parameters.require_static_hint = true;
        tr.parameters.max_hint_latency    = 0;
        auto run                          = [&](double seconds)
        {
            std::size_t n = 0;
            for (int i = 0; i < static_cast<int>(seconds * RATE); i++)
            {
                const auto o = standing(imu);
                if (i % 10 == 0)
                {
                    tr.add_motion_hint(imu.t, true);
                }
                n += tr.add(o);
            }
            return n;
        };
        ASSERT_GE_(run(5.0), 3U);
        imu.t -= 300.0;
        ASSERT_GE_(run(5.0), 3U);
    }
    std::cout << "hints: late sources and clock jumps\n";
}

// ---------------------------------------------------------------------------
/// Feeds standstill readings, with hints given `latency` seconds late every 0.1 s
/// ("moving" only at `tMoving`). Returns the [t0, t1] of the accepted windows.
std::vector<std::pair<double, double>> run_with_late_hints(
    ImuBiasTracker& tr, SimImu& imu, double seconds, double latency, double tMoving)
{
    std::vector<std::pair<double, double>> accepted;
    std::deque<double>                     pending;
    for (int i = 0; i < static_cast<int>(seconds * RATE); i++)
    {
        const auto o = standing(imu);
        if (i % 10 == 0)
        {
            pending.push_back(imu.t);
        }
        std::size_t n = tr.add(o);
        while (!pending.empty() && imu.t - pending.front() >= latency - 1e-9)
        {
            const double th       = pending.front();
            const bool   isMoving = std::abs(th - tMoving) < 0.05;
            n += tr.add_motion_hint(th, !isMoving).size();
            pending.pop_front();
        }
        if (n > 0)
        {
            const double t1 = tr.estimate(imu.label)->last_update;
            // A window spans [standstill_window, standstill_window + 1 sample]:
            accepted.emplace_back(t1 - tr.parameters.standstill_window - 1.1 / RATE, t1);
        }
    }
    return accepted;
}

void test_slow_motion_and_late_hints()
{
    // A wide curve at constant speed, without vibration: from the IMU alone it
    // looks like a bias (documented limitation), whatever its direction.
    for (const double wz : {0.005, -0.005})
    {
        ImuBiasTracker tr;
        SimImu         imu;
        imu.gyro_noise = 0.001;
        ASSERT_GT_(feed(tr, imu, 3.0, [wz](SimImu& s) { return s.sample(wz); }), 0U);
    }

    // A single "moving" hint, 0.3 s late among 10 Hz "static" ones: readings wait
    // for it, so no accepted window spans it (it holds until the next hint).
    {
        ImuBiasTracker tr;
        SimImu         imu;
        const double   tMoving = imu.t + 5.0;
        const auto     acc     = run_with_late_hints(tr, imu, 10.0, 0.3, tMoving);
        ASSERT_GE_(acc.size(), 6U);
        for (const auto& [t0, t1] : acc)
        {
            ASSERT_(t1 < tMoving || t0 > tMoving + 0.1 - 1e-6);
        }
        // Without waiting for hints, a window does span it:
        ImuBiasTracker tr2;
        SimImu         imu2;
        tr2.parameters.max_hint_latency = 0;
        const auto acc2                 = run_with_late_hints(tr2, imu2, 10.0, 0.3, tMoving);
        bool       spans                = false;
        for (const auto& [t0, t1] : acc2)
        {
            spans = spans || (t0 < tMoving && t1 > tMoving);
        }
        ASSERT_(spans);
    }

    // Two hints at the same time: "moving" wins.
    {
        ImuBiasTracker tr;
        SimImu         imu;
        tr.parameters.max_hint_latency = 0;
        tr.add_motion_hint(imu.t + 0.5, true);
        tr.add_motion_hint(imu.t + 0.5, false);
        ASSERT_EQUAL_(feed(tr, imu, 1.4, standing), 0U);
    }

    // A single "moving" hint just after an accepted window holds for
    // motion_hint_hold (1 s): it vetoes the windows that start before it expires.
    {
        ImuBiasTracker tr;
        SimImu         imu;
        tr.parameters.max_hint_latency = 0;
        ASSERT_EQUAL_(feed(tr, imu, 1.05, standing), 1U);
        tr.add_motion_hint(imu.t + 0.005, false);
        ASSERT_EQUAL_(feed(tr, imu, 1.9, standing), 0U);
        ASSERT_EQUAL_(feed(tr, imu, 0.3, standing), 1U);
    }

    // Requiring a static hint: none given, no updates even standing still...
    {
        ImuBiasTracker tr;
        tr.parameters.require_static_hint = true;
        tr.parameters.max_hint_latency    = 0;
        SimImu imu;
        ASSERT_EQUAL_(feed(tr, imu, 3.0, standing), 0U);
        // ...and with sparse (1 Hz) "static" hints, it works:
        std::size_t n = 0;
        for (int i = 0; i < 300; i++)
        {
            const auto o = standing(imu);
            if (i % 100 == 0)
            {
                tr.add_motion_hint(imu.t, true);
            }
            n += tr.add(o);
        }
        ASSERT_GE_(n, 2U);
    }

    // Hints stalled at the end of a dataset: flush() processes the waiting readings.
    {
        ImuBiasTracker tr;
        SimImu         imu;
        tr.add_motion_hint(imu.t, true);
        ASSERT_EQUAL_(feed(tr, imu, 1.4, standing), 0U);
        ASSERT_EQUAL_(tr.flush().size(), 1U);
    }
    std::cout << "motion hints: late, simultaneous, held, required and stalled ones\n";
}

// ---------------------------------------------------------------------------
void test_short_stops_at_random_phase()
{
    // Stops of 1.5 s between drives of random length: the sliding window finds
    // the 1 s of standstill in each one, wherever it falls.
    ImuBiasTracker                 tr;
    SimImu                         imu;
    mrpt::random::CRandomGenerator rng;
    rng.randomize(7);
    std::size_t stopsWithUpdate = 0;
    const int   nStops          = 40;
    for (int i = 0; i < nStops; i++)
    {
        feed(tr, imu, 3.0 + rng.drawUniform(0.0, 1.0), driving);
        stopsWithUpdate += feed(tr, imu, 1.5, standing) > 0 ? 1 : 0;
    }
    ASSERT_GE_(stopsWithUpdate, 38U);
    std::cout << "short stops: " << stopsWithUpdate << "/" << nStops << " with an update\n";
}

// ---------------------------------------------------------------------------
void test_covariance_consistency()
{
    // Over many runs, the gyro bias error matches the reported covariance (NEES).
    double    neesSum = 0;
    const int runs    = 60;
    for (int run = 0; run < runs; run++)
    {
        ImuBiasTracker tr;
        SimImu         imu(100 + run);
        feed(tr, imu, 6.0, standing);
        const auto            e   = *tr.estimate(imu.label);
        const Eigen::Vector3d err = Eigen::Vector3d(
            e.bias_gyro.x - imu.bias_gyro.x, e.bias_gyro.y - imu.bias_gyro.y,
            e.bias_gyro.z - imu.bias_gyro.z);
        neesSum += err.dot(e.cov_gyro.asEigen().inverse() * err);
    }
    const double nees = neesSum / runs;  // 3 if consistent
    ASSERT_GT_(nees, 1.0);
    ASSERT_LT_(nees, 5.0);
    std::cout << "covariance: mean NEES " << nees << " (3 dof)\n";
}

// ---------------------------------------------------------------------------
void test_stops_at_different_tilts()
{
    // Stops on slopes of different directions make the whole accelerometer
    // bias observable, not only along gravity.
    ImuBiasTracker tr;
    SimImu         imu;
    for (const auto& [pitch, roll] :
         std::vector<std::pair<double, double>>{{0, 0}, {15, 0}, {-15, 0}, {0, 15}, {0, -15}})
    {
        imu.mounting =
            mrpt::poses::CPose3D::FromYawPitchRoll(0, mrpt::DEG2RAD(pitch), mrpt::DEG2RAD(roll));
        feed(tr, imu, 20.0, standing);
        feed(tr, imu, 5.0, driving);
    }
    const auto e = *tr.estimate(imu.label);
    for (int i = 0; i < 3; i++)
    {
        ASSERT_LT_(std::abs(e.bias_acc[i] - imu.bias_acc[i]), 0.01);
        ASSERT_LT_(std::sqrt(e.cov_acc(i, i)), 0.05);
    }
    std::cout << "tilts: full accel bias within " << norm3(e.bias_acc - imu.bias_acc) << " m/s^2\n";
}

// ---------------------------------------------------------------------------
void test_floors_hold_along_gravity()
{
    // After a long stop on a slope, the accelerometer bias along gravity is not
    // reported as more certain than the floor (errors that do not average out).
    ImuBiasTracker tr;
    SimImu         imu;
    imu.mounting =
        mrpt::poses::CPose3D::FromYawPitchRoll(0, mrpt::DEG2RAD(10.0), mrpt::DEG2RAD(5.0));
    feed(tr, imu, 600.0, standing);
    const auto            e = *tr.estimate(imu.label);
    const Eigen::Vector3d u = gravity_dir(imu);
    ASSERT_GE_(std::sqrt(u.dot(e.cov_acc.asEigen() * u)), 0.999 * tr.parameters.bias_acc_min_sigma);
    const Eigen::Vector3d z = Eigen::Vector3d::UnitZ();
    ASSERT_GE_(
        std::sqrt(z.dot(e.cov_gyro.asEigen() * z)), 0.999 * tr.parameters.bias_gyro_min_sigma);
    std::cout << "floors: hold along gravity\n";
}

// ---------------------------------------------------------------------------
void test_invalid_and_unordered_readings()
{
    // Non-finite accelerometer values: dropped, the estimate stays finite.
    {
        ImuBiasTracker tr;
        SimImu         imu;
        std::size_t    n = 0;
        for (int i = 0; i < 500; i++)
        {
            auto o = standing(imu);
            if (i % 50 == 0)
            {
                o.set(mrpt::obs::IMU_Y_ACC, std::numeric_limits<double>::quiet_NaN());
            }
            n += tr.add(o);
        }
        ASSERT_GE_(n, 4U);
        const auto e = *tr.estimate(imu.label);
        for (int i = 0; i < 3; i++)
        {
            ASSERT_(std::isfinite(e.bias_acc[i]));
            ASSERT_(std::isfinite(e.cov_acc(i, i)));
        }
        const double estU =
            gravity_dir(imu).dot(Eigen::Vector3d(e.bias_acc.x, e.bias_acc.y, e.bias_acc.z));
        const double trueU =
            gravity_dir(imu).dot(Eigen::Vector3d(imu.bias_acc.x, imu.bias_acc.y, imu.bias_acc.z));
        ASSERT_LT_(std::abs(estU - trueU), 0.02);
    }
    // Garbage values (NaN or absurd gyro, absurd accelerometer, no time stamp) are
    // ignored: the result is the same as without them.
    {
        auto run = [](bool withGarbage)
        {
            ImuBiasTracker tr;
            SimImu         imu;
            for (int i = 0; i < 500; i++)
            {
                const auto o = standing(imu);
                tr.add(o);
                if (withGarbage && i % 50 == 0)
                {
                    auto bad = o;
                    bad.timestamp += std::chrono::milliseconds(1);
                    bad.set(mrpt::obs::IMU_WZ, (i % 100 == 0) ? 1e200 : std::nan(""));
                    tr.add(bad);
                    auto noTime      = o;
                    noTime.timestamp = mrpt::Clock::time_point();
                    tr.add(noTime);
                }
            }
            // One absurd accelerometer reading, then a stop:
            auto o = standing(imu);
            if (withGarbage)
            {
                o.set(mrpt::obs::IMU_Z_ACC, 3e38);
            }
            tr.add(o);
            const auto n = feed(tr, imu, 3.0, standing);
            return std::make_pair(n, *tr.estimate(imu.label));
        };
        const auto [nClean, eClean] = run(false);
        const auto [nBad, eBad]     = run(true);
        ASSERT_GE_(nBad, 2U);
        ASSERT_EQUAL_(nBad, nClean);
        ASSERT_LT_(norm3(eBad.bias_gyro - eClean.bias_gyro), 1e-6);
        ASSERT_(std::isfinite(eBad.cov_gyro(2, 2)));
        ASSERT_LT_(std::sqrt(eBad.cov_gyro(2, 2)), 1e-3);
    }
    // Duplicated and slightly older readings (e.g. re-sent by a driver) are
    // dropped: the same result as without them.
    {
        auto run = [](bool withRepeats)
        {
            ImuBiasTracker tr;
            SimImu         imu;
            for (int i = 0; i < 500; i++)
            {
                const auto o = standing(imu);
                tr.add(o);
                if (withRepeats && i % 10 == 0)
                {
                    tr.add(o);
                    auto older = o;
                    older.timestamp -= std::chrono::milliseconds(5);
                    tr.add(older);
                }
            }
            return *tr.estimate("imu");
        };
        const auto eA = run(false);
        const auto eB = run(true);
        ASSERT_EQUAL_(eA.standstill_updates, eB.standstill_updates);
        ASSERT_EQUAL_(eA.bias_gyro.z, eB.bias_gyro.z);
    }
    // A jump back in time (a replayed dataset): the IMU restarts, keeping its bias.
    {
        ImuBiasTracker tr;
        SimImu         imu;
        feed(tr, imu, 5.0, standing);
        const auto before = *tr.estimate(imu.label);
        imu.t -= 100.0;
        const auto n     = feed(tr, imu, 3.0, standing);
        const auto after = *tr.estimate(imu.label);
        ASSERT_GE_(n, 2U);
        ASSERT_EQUAL_(after.standstill_updates, before.standstill_updates + n);
        ASSERT_LT_(after.stamp, before.stamp);
        ASSERT_LT_(norm3(after.bias_gyro - imu.bias_gyro), 4e-4);
    }
    // A gyro-only IMU: gyro bias updated, the accelerometer one stays the prior.
    {
        ImuBiasTracker tr;
        SimImu         imu;
        imu.with_acc = false;
        ASSERT_GE_(feed(tr, imu, 5.0, standing), 4U);
        const auto e = *tr.estimate(imu.label);
        ASSERT_LT_(norm3(e.bias_gyro - imu.bias_gyro), 5e-4);
        ASSERT_EQUAL_(norm3(e.bias_acc), 0.0);
    }
    // An external measurement for an IMU before its first reading does not block
    // its readings:
    {
        ImuBiasTracker              tr;
        SimImu                      imu;
        mrpt::math::CMatrixDouble33 R;
        R.setDiagonal(1e-4);
        tr.add_measurement(
            imu.label, imu.t + 0.9, std::make_pair(TVector3D(0, 0, 0), R), std::nullopt);
        ASSERT_EQUAL_(feed(tr, imu, 1.5, standing), 1U);
    }
    std::cout << "invalid, duplicated and unordered readings: handled\n";
}

// ---------------------------------------------------------------------------
void test_external_measurements()
{
    ImuBiasTracker tr;
    tr.parameters.bias_gyro_prior_sigma = 0.02;
    tr.parameters.bias_gyro_random_walk = 0;
    tr.parameters.bias_gyro_min_sigma   = 0;
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

    // A measurement one hour older than the estimate counts as less certain, by
    // the random walk since its time, and does not move last_update back:
    ImuBiasTracker tr2;
    tr2.parameters.bias_gyro_min_sigma = 0;
    const double                q      = mrpt::square(tr2.parameters.bias_gyro_random_walk);
    mrpt::math::CMatrixDouble33 R;
    R.setDiagonal(1e-8);
    tr2.add_measurement("imu", 5000.0, std::make_pair(TVector3D(0, 0, 0), R), std::nullopt);
    const double p0 = tr2.estimate("imu")->cov_gyro(2, 2);
    tr2.add_measurement("imu", 1400.0, std::make_pair(TVector3D(0, 0, 0.01), R), std::nullopt);
    const auto   e2   = *tr2.estimate("imu");
    const double rOld = 1e-8 + q * 3600.0;
    ASSERT_NEAR_(e2.cov_gyro(2, 2), 1.0 / (1.0 / p0 + 1.0 / rOld), 1e-15);
    ASSERT_NEAR_(e2.bias_gyro.z, 0.01 * p0 / (p0 + rOld), 1e-12);
    ASSERT_EQUAL_(e2.last_update, 5000.0);
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
    SimImu         a(1);
    SimImu         b(2);
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

    // Hints are about the vehicle: one "moving" hint vetoes both IMUs.
    tr.parameters.max_hint_latency = 0;
    const auto na0                 = tr.estimate("imu")->standstill_updates;
    const auto nb0                 = tr.estimate("imu2")->standstill_updates;
    for (int i = 0; i < 300; i++)
    {
        tr.add_motion_hint(a.t + 1.0 / RATE, false);
        tr.add(a.sample(0));
        tr.add(b.sample(0));
    }
    ASSERT_EQUAL_(tr.estimate("imu")->standstill_updates, na0);
    ASSERT_EQUAL_(tr.estimate("imu2")->standstill_updates, nb0);
    std::cout << "two IMUs: independent estimates, shared hints\n";
}

// ---------------------------------------------------------------------------
void test_yaml_roundtrip()
{
    ImuBiasTracker::Parameters p;
    p.standstill_window      = 2.5;
    p.require_static_hint    = true;
    p.bias_gyro_random_walk  = 1e-4;
    p.max_hint_latency       = 0.8;
    mrpt::containers::yaml y = mrpt::containers::yaml::Map();
    p.save_to(y);
    ImuBiasTracker::Parameters q;
    q.load_from(y);
    ASSERT_EQUAL_(q.standstill_window, 2.5);
    ASSERT_(q.require_static_hint);
    ASSERT_EQUAL_(q.bias_gyro_random_walk, 1e-4);
    ASSERT_EQUAL_(q.max_hint_latency, 0.8);
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
        test_large_bias_accepted();
        test_bias_changes_are_reacquired();
        test_reacquisition_robustness();
        test_hint_sources_and_clock_jumps();
        test_slow_motion_and_late_hints();
        test_short_stops_at_random_phase();
        test_covariance_consistency();
        test_stops_at_different_tilts();
        test_floors_hold_along_gravity();
        test_invalid_and_unordered_readings();
        test_external_measurements();
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
