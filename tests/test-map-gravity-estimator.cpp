/*               _
 _ __ ___   ___ | | __ _
| '_ ` _ \ / _ \| |/ _` | Modular Optimization framework for
| | | | | | (_) | | (_| | Localization and mApping (MOLA)
|_| |_| |_|\___/|_|\__,_| https://github.com/MOLAorg/mola

 Copyright (C) 2018-2026 Jose Luis Blanco, University of Almeria,
                         and individual contributors.
 SPDX-License-Identifier: GPL-3.0
 See LICENSE for full license information.
*/

/**
 * @file   test-map-gravity-estimator.cpp
 * @brief  Unit tests for MapGravityEstimator.
 * @author Jose Luis Blanco Claraco
 * @date   Jul 22, 2026
 *
 * The central case is the one that motivates the whole feature: a map frame
 * that has accumulated a tilt, on a trajectory that NEVER STOPS, with a raw
 * accel+gyro IMU and no absolute attitude reference. Averaging accelerometer
 * samples cannot recover the tilt there; preintegration can.
 */

#include <mola_imu_preintegration/MapGravityEstimator.h>
#include <mrpt/core/exceptions.h>

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <random>

using namespace mola::imu;
using mrpt::math::TVector3D;

namespace
{
constexpr double GRAVITY = 9.81;

/** A synthetic world.
 *
 *  Ground truth lives in a gravity-aligned frame {W}. The odometry ("map")
 *  frame {M} is {W} rotated by a fixed tilt, which is what the estimator has
 *  to discover: everything the estimator is fed (attitudes, velocities) is
 *  expressed in {M}, while the IMU of course senses true gravity.
 */
struct SyntheticWorld
{
    /// The tilt of the map frame: R_M_W, i.e. p_M = R_map_from_world * p_W.
    TVector3D map_tilt{0, 0, 0};  //!< rotation vector [rad]

    /// Body attitude in {W}: a constant tilt plus a steady yaw rotation, so
    /// there IS rotation excitation (needed for accel-bias observability).
    TVector3D body_tilt{0.05, 0.09, 0.0};
    double    yaw_rate = 0.30;  //!< [rad/s]

    /// Trajectory in {W}: a Lissajous curve, so the vehicle never stops and
    /// acceleration is nontrivial in all axes.
    TVector3D amplitude{3.0, 2.0, 0.5};
    TVector3D frequency{0.6, 0.9, 1.4};

    /// True IMU biases, in the body frame.
    TVector3D bias_acc{0, 0, 0};
    TVector3D bias_gyro{0, 0, 0};

    SO3 R_map_from_world() const { return so3_exp(map_tilt); }

    /// True body attitude in {W}.
    SO3 R_world_body(double t) const
    {
        SO3 out;
        out.asEigen() = so3_exp(body_tilt).asEigen() * so3_exp({0, 0, yaw_rate * t}).asEigen();
        return out;
    }

    /// Body attitude as the odometry reports it, i.e. in {M}.
    SO3 R_map_body(double t) const
    {
        SO3 out;
        out.asEigen() = R_map_from_world().asEigen() * R_world_body(t).asEigen();
        return out;
    }

    TVector3D v_world(double t) const
    {
        return {
            amplitude.x * frequency.x * std::cos(frequency.x * t),
            amplitude.y * frequency.y * std::cos(frequency.y * t),
            amplitude.z * frequency.z * std::cos(frequency.z * t)};
    }

    /// Velocity as the odometry reports it, i.e. in {M}.
    TVector3D v_map(double t) const
    {
        const Eigen::Vector3d v = R_map_from_world().asEigen() * toE(v_world(t));
        return {v.x(), v.y(), v.z()};
    }

    TVector3D a_world(double t) const
    {
        return {
            -amplitude.x * frequency.x * frequency.x * std::sin(frequency.x * t),
            -amplitude.y * frequency.y * frequency.y * std::sin(frequency.y * t),
            -amplitude.z * frequency.z * frequency.z * std::sin(frequency.z * t)};
    }

    /// Ideal accelerometer reading (proper acceleration in body frame), plus
    /// the true bias.
    TVector3D accel_reading(double t) const
    {
        const TVector3D       aw = a_world(t);
        const Eigen::Vector3d d(aw.x, aw.y, aw.z + GRAVITY);
        const Eigen::Vector3d b = R_world_body(t).asEigen().transpose() * d;
        return TVector3D{b.x(), b.y(), b.z()} + bias_acc;
    }

    /// Ideal gyroscope reading, plus the true bias.
    /// With `R(t) = Exp(body_tilt) * Rz(yaw_rate*t)`, we have
    /// `R^T dR/dt = [ (0,0,yaw_rate) ]_x` exactly, i.e. the body-frame angular
    /// rate is constant and equal to (0,0,yaw_rate), independently of the
    /// constant tilt. (Rotating a world-frame rate into the body instead would
    /// make the synthetic gyro inconsistent with the synthetic attitude.)
    TVector3D gyro_reading(double /*t*/) const { return TVector3D{0, 0, yaw_rate} + bias_gyro; }

    /// The true gravity vector expressed in {M}: what the estimator must find.
    TVector3D gravity_in_map() const
    {
        const Eigen::Vector3d g = R_map_from_world().asEigen() * Eigen::Vector3d(0, 0, -GRAVITY);
        return {g.x(), g.y(), g.z()};
    }

    static Eigen::Vector3d toE(const TVector3D& v) { return {v.x, v.y, v.z}; }
};

/// Builds one interval of the given duration, preintegrating at `rate` Hz.
MapGravityEstimator::Interval makeInterval(
    const SyntheticWorld& w, double t0, double t1, double rate,
    const ImuIntegrationParams& imuParams, std::mt19937* rng = nullptr, double accelNoise = 0,
    double gyroNoise = 0)
{
    ImuPreintegrator pim(imuParams);
    pim.reset_integration();

    const double dt = 1.0 / rate;
    const int    n  = static_cast<int>(std::llround((t1 - t0) / dt));

    std::normal_distribution<double> nd(0.0, 1.0);

    for (int k = 0; k < n; k++)
    {
        const double t = t0 + k * dt;
        TVector3D    a = w.accel_reading(t);
        TVector3D    g = w.gyro_reading(t);
        if (rng != nullptr)
        {
            a = a + TVector3D{nd(*rng) * accelNoise, nd(*rng) * accelNoise, nd(*rng) * accelNoise};
            g = g + TVector3D{nd(*rng) * gyroNoise, nd(*rng) * gyroNoise, nd(*rng) * gyroNoise};
        }
        pim.integrate_measurement(a, g, dt);
    }

    MapGravityEstimator::Interval itv;
    itv.t_from = t0;
    // The interval is delimited by the KEYFRAME stamps, as in the real
    // wiring; whether the IMU samples actually cover it is what the coverage
    // guard is there to check.
    itv.t_to   = t1;
    itv.delta  = pim.current_state();
    itv.R_from = w.R_map_body(itv.t_from);
    itv.R_to   = w.R_map_body(itv.t_to);
    itv.v_from = w.v_map(itv.t_from);
    itv.v_to   = w.v_map(itv.t_to);
    return itv;
}

/// Angle between two vectors [deg].
double angleBetweenDeg(const TVector3D& a, const TVector3D& b)
{
    const double na = a.norm();
    const double nb = b.norm();
    ASSERT_GT_(na, 0);
    ASSERT_GT_(nb, 0);
    const double c = std::clamp((a.x * b.x + a.y * b.y + a.z * b.z) / (na * nb), -1.0, 1.0);
    return mrpt::RAD2DEG(std::acos(c));
}

/// Runs a full sequence of intervals through the estimator.
MapGravityEstimator::Result runSequence(
    MapGravityEstimator& est, const SyntheticWorld& w, const ImuIntegrationParams& imuParams,
    int numIntervals, double intervalDuration = 0.5, double rate = 200.0,
    std::mt19937* rng = nullptr, double accelNoise = 0, double gyroNoise = 0)
{
    for (int k = 0; k < numIntervals; k++)
    {
        const double t0 = k * intervalDuration;
        const auto   itv =
            makeInterval(w, t0, t0 + intervalDuration, rate, imuParams, rng, accelNoise, gyroNoise);

        const bool ok = est.add_interval(itv);
        ASSERTMSG_(ok, std::string("interval rejected: ") + est.last_rejection_reason());
    }
    ASSERT_(est.solve());
    ASSERT_(est.latest_result().has_value());
    return *est.latest_result();
}

// -----------------------------------------------------------------------
// THE key test: a tilted map, a trajectory that never stops, raw IMU only.
// -----------------------------------------------------------------------
void test_recovers_map_tilt_while_always_moving()
{
    SyntheticWorld w;
    w.map_tilt = {mrpt::DEG2RAD(2.5), mrpt::DEG2RAD(-1.8), 0.0};  // ~3.1 deg of tilt

    MapGravityEstimator est;
    const auto          res = runSequence(est, w, ImuIntegrationParams(), 20);

    const double dirErrDeg = angleBetweenDeg(res.gravity_in_map, w.gravity_in_map());

    ASSERTMSG_(
        dirErrDeg < 0.05,
        mrpt::format(
            "gravity direction error too large: %.4f deg (estimated %s, true %s)", dirErrDeg,
            res.gravity_in_map.asString().c_str(), w.gravity_in_map().asString().c_str()));

    ASSERT_NEAR_(res.gravity_in_map.norm(), GRAVITY, 0.02);
    ASSERT_(res.converged);

    // The correction must actually level the map: applying it to the measured
    // gravity has to yield (0,0,-g).
    const Eigen::Vector3d leveled =
        res.correction.asEigen() *
        Eigen::Vector3d(res.gravity_in_map.x, res.gravity_in_map.y, res.gravity_in_map.z);
    ASSERT_NEAR_(leveled.x(), 0.0, 1e-6);
    ASSERT_NEAR_(leveled.y(), 0.0, 1e-6);
    ASSERT_LT_(leveled.z(), 0.0);

    std::cout << "test_recovers_map_tilt_while_always_moving passed (dir err " << dirErrDeg
              << " deg, tilt " << mrpt::RAD2DEG(res.tilt) << " deg, " << res.iterations
              << " iters)." << std::endl;
}

// A level map must be reported as level: no phantom tilt.
void test_level_map_yields_no_correction()
{
    SyntheticWorld w;  // map_tilt = 0

    MapGravityEstimator est;
    const auto          res = runSequence(est, w, ImuIntegrationParams(), 20);

    ASSERTMSG_(
        mrpt::RAD2DEG(res.tilt) < 0.05,
        mrpt::format("phantom tilt on a level map: %.4f deg", mrpt::RAD2DEG(res.tilt)));

    std::cout << "test_level_map_yields_no_correction passed (tilt " << mrpt::RAD2DEG(res.tilt)
              << " deg)." << std::endl;
}

// -----------------------------------------------------------------------
// Bias observability.
// -----------------------------------------------------------------------
void test_gyro_bias_observability()
{
    SyntheticWorld w;
    w.map_tilt  = {mrpt::DEG2RAD(1.0), mrpt::DEG2RAD(1.5), 0};
    w.bias_gyro = {0.004, -0.003, 0.002};  // rad/s

    // (a) With a DELIBERATELY LOOSE prior, the bias must be recovered almost
    //     exactly. This is what proves the residual and its Jacobian are
    //     correctly scaled, as opposed to merely "close".
    {
        MapGravityEstimator est;
        est.parameters.bias_gyro_prior_sigma = 1.0;  // effectively no prior

        const auto   res = runSequence(est, w, ImuIntegrationParams(), 20);
        const double err = (res.bias_gyro - w.bias_gyro).norm();

        ASSERTMSG_(
            err < 5e-5, mrpt::format(
                            "gyro bias not recovered with a loose prior: estimated %s vs true %s "
                            "(err %g rad/s)",
                            res.bias_gyro.asString().c_str(), w.bias_gyro.asString().c_str(), err));

        std::cout << "  loose prior: bias err " << err << " rad/s" << std::endl;
    }

    // (b) With the DEFAULT (tight, zero-anchored) prior, the estimate is
    //     deliberately shrunk toward zero: that is the prior working as
    //     designed, and is what keeps the bias bounded when the data is weak.
    //     It must still recover most of the true bias, and must not overshoot.
    {
        MapGravityEstimator est;
        const auto          res = runSequence(est, w, ImuIntegrationParams(), 20);

        const double trueNorm = w.bias_gyro.norm();
        const double recovered =
            (res.bias_gyro.x * w.bias_gyro.x + res.bias_gyro.y * w.bias_gyro.y +
             res.bias_gyro.z * w.bias_gyro.z) /
            (trueNorm * trueNorm);

        ASSERTMSG_(
            recovered > 0.8 && recovered <= 1.0,
            mrpt::format(
                "gyro bias shrinkage out of the expected range: recovered fraction %.3f "
                "(estimated %s vs true %s)",
                recovered, res.bias_gyro.asString().c_str(), w.bias_gyro.asString().c_str()));

        // and the tilt must still come out right despite the bias:
        const double dirErrDeg = angleBetweenDeg(res.gravity_in_map, w.gravity_in_map());
        ASSERT_LT_(dirErrDeg, 0.1);

        std::cout << "  default prior: recovered fraction " << recovered << ", dir err "
                  << dirErrDeg << " deg" << std::endl;
    }

    std::cout << "test_gyro_bias_observability passed." << std::endl;
}

void test_accel_bias_observability_with_rotation_excitation()
{
    SyntheticWorld w;
    w.map_tilt = {mrpt::DEG2RAD(2.0), mrpt::DEG2RAD(-1.0), 0};
    w.bias_acc = {0.05, -0.04, 0.03};  // m/s^2
    w.yaw_rate = 0.6;  // plenty of rotation excitation

    // (a) Loose prior: rotation excitation makes the accelerometer bias
    //     separable from gravity (the bias rotates with the body while gravity
    //     stays fixed in the map frame), so it must be recovered accurately.
    {
        MapGravityEstimator est;
        est.parameters.bias_acc_prior_sigma = 10.0;  // effectively no prior

        const auto   res = runSequence(est, w, ImuIntegrationParams(), 30);
        const double err = (res.bias_acc - w.bias_acc).norm();

        // The floor here is the first-order discretization error of the
        // Euler-forward preintegration at 200 Hz, which biases dV slightly and
        // is partly absorbed by the accel bias. It is ~1e-3 m/s^2, far below
        // any real accelerometer's bias instability.
        ASSERTMSG_(
            err < 5e-3, mrpt::format(
                            "accel bias not recovered with a loose prior: estimated %s vs true %s "
                            "(err %g m/s^2)",
                            res.bias_acc.asString().c_str(), w.bias_acc.asString().c_str(), err));

        const double dirErrDeg = angleBetweenDeg(res.gravity_in_map, w.gravity_in_map());
        ASSERT_LT_(dirErrDeg, 0.05);

        std::cout << "  loose prior: bias err " << err << " m/s^2, dir err " << dirErrDeg << " deg"
                  << std::endl;
    }

    // (b) Default prior: the accelerometer bias is shrunk noticeably harder
    //     than the gyro one, because it trades directly against the gravity
    //     direction we are trying to estimate. What matters is that it moves
    //     the RIGHT WAY and stays bounded, and that the tilt is still right.
    {
        MapGravityEstimator est;
        const auto          res = runSequence(est, w, ImuIntegrationParams(), 30);

        const double trueNorm  = w.bias_acc.norm();
        const double recovered = (res.bias_acc.x * w.bias_acc.x + res.bias_acc.y * w.bias_acc.y +
                                  res.bias_acc.z * w.bias_acc.z) /
                                 (trueNorm * trueNorm);

        ASSERTMSG_(
            recovered > 0.3 && recovered <= 1.0,
            mrpt::format(
                "accel bias shrinkage out of the expected range: recovered fraction %.3f "
                "(estimated %s vs true %s)",
                recovered, res.bias_acc.asString().c_str(), w.bias_acc.asString().c_str()));

        const double dirErrDeg = angleBetweenDeg(res.gravity_in_map, w.gravity_in_map());
        ASSERT_LT_(dirErrDeg, 0.3);

        std::cout << "  default prior: recovered fraction " << recovered << ", dir err "
                  << dirErrDeg << " deg" << std::endl;
    }

    std::cout << "test_accel_bias_observability_with_rotation_excitation passed." << std::endl;
}

// Without rotation excitation the accel bias and a tilt are indistinguishable.
// The point of this test is that the estimator degrades GRACEFULLY (the prior
// dominates, nothing diverges), not that it magically recovers the truth.
void test_no_rotation_excitation_degrades_gracefully()
{
    SyntheticWorld w;
    w.map_tilt  = {mrpt::DEG2RAD(1.5), 0, 0};
    w.bias_acc  = {0.05, -0.04, 0.0};
    w.yaw_rate  = 0.0;  // no rotation at all
    w.body_tilt = {0, 0, 0};  // and level body

    MapGravityEstimator est;
    const auto          res = runSequence(est, w, ImuIntegrationParams(), 20);

    ASSERT_(std::isfinite(res.tilt));
    ASSERT_(std::isfinite(res.bias_acc.norm()));
    // The prior must keep the bias bounded rather than letting it run away:
    ASSERT_LT_(res.bias_acc.norm(), 3.0 * est.parameters.bias_acc_prior_sigma);
    // And the magnitude constraint must keep |g| sane:
    ASSERT_NEAR_(res.gravity_in_map.norm(), GRAVITY, 0.5);

    std::cout << "test_no_rotation_excitation_degrades_gracefully passed (tilt "
              << mrpt::RAD2DEG(res.tilt) << " deg, |b_a| " << res.bias_acc.norm() << ")."
              << std::endl;
}

// -----------------------------------------------------------------------
// Robustness.
// -----------------------------------------------------------------------
void test_huber_rejects_corrupted_velocity()
{
    SyntheticWorld w;
    w.map_tilt = {mrpt::DEG2RAD(2.0), mrpt::DEG2RAD(-1.5), 0};

    // Reference run, all intervals clean:
    MapGravityEstimator estClean;
    const auto          clean = runSequence(estClean, w, ImuIntegrationParams(), 20);

    // Same, but one interval has a grossly wrong odometry velocity:
    MapGravityEstimator est;
    for (int k = 0; k < 20; k++)
    {
        const double t0  = k * 0.5;
        auto         itv = makeInterval(w, t0, t0 + 0.5, 200.0, ImuIntegrationParams());
        if (k == 10)
        {
            itv.v_to = itv.v_to + TVector3D{5.0, -4.0, 3.0};  // a badly wrong LO velocity
        }
        ASSERT_(est.add_interval(itv));
    }
    ASSERT_(est.solve());
    const auto res = *est.latest_result();

    const double dirErrDeg = angleBetweenDeg(res.gravity_in_map, w.gravity_in_map());
    ASSERTMSG_(
        dirErrDeg < 0.5, mrpt::format(
                             "outlier dragged the estimate: %.4f deg (clean run: %.4f deg)",
                             dirErrDeg, angleBetweenDeg(clean.gravity_in_map, w.gravity_in_map())));

    std::cout << "test_huber_rejects_corrupted_velocity passed (dir err " << dirErrDeg << " deg)."
              << std::endl;
}

void test_survives_realistic_noise()
{
    SyntheticWorld w;
    w.map_tilt  = {mrpt::DEG2RAD(2.0), mrpt::DEG2RAD(-2.0), 0};
    w.bias_gyro = {0.002, 0.001, -0.001};

    const double rate   = 200.0;
    const double sigmaG = 3e-3;  // rad/s/sqrt(Hz)
    const double sigmaA = 3e-2;  // m/s^2/sqrt(Hz)

    ImuIntegrationParams imuParams;
    imuParams.cov_gyro.setDiagonal(sigmaG * sigmaG);
    imuParams.cov_acc.setDiagonal(sigmaA * sigmaA);

    std::mt19937 rng(1234);

    MapGravityEstimator est;
    const auto          res = runSequence(
                 est, w, imuParams, 30, 0.5, rate, &rng, sigmaA * std::sqrt(rate), sigmaG * std::sqrt(rate));

    const double dirErrDeg = angleBetweenDeg(res.gravity_in_map, w.gravity_in_map());
    ASSERTMSG_(
        dirErrDeg < 1.0,
        mrpt::format("noisy-run gravity direction error too large: %.4f deg", dirErrDeg));
    ASSERT_(res.converged);

    std::cout << "test_survives_realistic_noise passed (dir err " << dirErrDeg
              << " deg, reported sigma pitch/roll " << mrpt::RAD2DEG(res.pitch_sigma) << "/"
              << mrpt::RAD2DEG(res.roll_sigma) << " deg)." << std::endl;
}

// -----------------------------------------------------------------------
// Residual Jacobians.
//
// The two analytic Jacobians are the most error-prone code in the class, and
// the recovery tests above would still pass with a merely APPROXIMATE one (a
// wrong Jacobian costs iterations, not necessarily accuracy). So check them
// against numerical differentiation of the very same residual functions the
// solver calls.
// -----------------------------------------------------------------------
void checkJacobianColumn(
    const char* what, int col, const TVector3D& numeric, const Eigen::Vector3d& analytic)
{
    const Eigen::Vector3d numE(numeric.x, numeric.y, numeric.z);
    const double          err = (numE - analytic).norm();

    ASSERTMSG_(
        err < 1e-6, mrpt::format(
                        "%s: analytic Jacobian column %i disagrees with numerical "
                        "differentiation (err %g): analytic [%g %g %g] vs numeric [%g %g %g]",
                        what, col, err, analytic.x(), analytic.y(), analytic.z(), numeric.x,
                        numeric.y, numeric.z));
}

void test_residual_jacobians()
{
    SyntheticWorld w;
    w.map_tilt  = {mrpt::DEG2RAD(2.0), mrpt::DEG2RAD(-1.5), 0};
    w.bias_acc  = {0.03, -0.02, 0.04};
    w.bias_gyro = {0.005, -0.004, 0.003};
    w.yaw_rate  = 0.5;

    const auto itv = makeInterval(w, 1.0, 1.5, 200.0, ImuIntegrationParams());

    // Linearization point deliberately away from both the truth and the
    // interval's own bias linearization point, so that every first-order bias
    // term is actually exercised (they vanish at zero bias increment).
    const TVector3D g{0.35, -0.22, -9.65};
    const TVector3D ba{0.02, -0.03, 0.015};
    const TVector3D bg{0.004, -0.006, 0.002};

    constexpr double EPS = 1e-6;

    // ---- (a) gravity residual ------------------------------------------
    {
        mrpt::math::CMatrixDouble33 Jg;
        mrpt::math::CMatrixDouble33 Ja;
        mrpt::math::CMatrixDouble33 Jb;
        MapGravityEstimator::gravity_residual(itv, g, ba, bg, Jg, Ja, Jb);

        for (int i = 0; i < 3; i++)
        {
            TVector3D d{0, 0, 0};
            d[i] = EPS;

            const auto rgP = MapGravityEstimator::gravity_residual(itv, g + d, ba, bg);
            const auto rgM = MapGravityEstimator::gravity_residual(itv, g - d, ba, bg);
            checkJacobianColumn(
                "gravity residual d/d g", i, (rgP - rgM) * (1.0 / (2 * EPS)), Jg.asEigen().col(i));

            const auto raP = MapGravityEstimator::gravity_residual(itv, g, ba + d, bg);
            const auto raM = MapGravityEstimator::gravity_residual(itv, g, ba - d, bg);
            checkJacobianColumn(
                "gravity residual d/d b_a", i, (raP - raM) * (1.0 / (2 * EPS)),
                Ja.asEigen().col(i));

            const auto rbP = MapGravityEstimator::gravity_residual(itv, g, ba, bg + d);
            const auto rbM = MapGravityEstimator::gravity_residual(itv, g, ba, bg - d);
            checkJacobianColumn(
                "gravity residual d/d b_g", i, (rbP - rbM) * (1.0 / (2 * EPS)),
                Jb.asEigen().col(i));
        }
    }

    // ---- (b) relative-rotation residual --------------------------------
    {
        mrpt::math::CMatrixDouble33 Jb;
        MapGravityEstimator::rotation_residual(itv, bg, Jb);

        for (int i = 0; i < 3; i++)
        {
            TVector3D d{0, 0, 0};
            d[i] = EPS;

            const auto rP = MapGravityEstimator::rotation_residual(itv, bg + d);
            const auto rM = MapGravityEstimator::rotation_residual(itv, bg - d);
            checkJacobianColumn(
                "rotation residual d/d b_g", i, (rP - rM) * (1.0 / (2 * EPS)), Jb.asEigen().col(i));
        }
    }

    std::cout << "test_residual_jacobians passed." << std::endl;
}

// -----------------------------------------------------------------------
// Guards and bookkeeping.
// -----------------------------------------------------------------------

// A weak-but-usable estimate must be REPORTED with an honest large sigma, not
// filtered out. This pins the contract that replaced the former
// `max_tilt_sigma_deg` gate, whose default silently disabled the estimator on
// real hardware, where per-window tilt sigmas of several degrees are normal.
void test_weak_data_is_reported_not_gated()
{
    SyntheticWorld w;
    w.map_tilt = {mrpt::DEG2RAD(2.0), mrpt::DEG2RAD(-1.5), 0};

    MapGravityEstimator est;

    // Only the minimum number of intervals, and an external odometry whose
    // velocities are declared very uncertain:
    for (int k = 0; k < 4; k++)
    {
        const double t0  = k * 0.5;
        auto         itv = makeInterval(w, t0, t0 + 0.5, 200.0, ImuIntegrationParams());

        itv.cov_v_from.setDiagonal(mrpt::square(2.0));  // 2 m/s
        itv.cov_v_to.setDiagonal(mrpt::square(2.0));

        ASSERT_(est.add_interval(itv));
    }
    ASSERT_(est.solve());
    const auto res = *est.latest_result();

    ASSERTMSG_(res.converged, "a usable, if weak, estimate must be reported as converged");

    const double sigmaDeg = mrpt::RAD2DEG(std::max(res.pitch_sigma, res.roll_sigma));
    ASSERTMSG_(
        sigmaDeg > 3.0,
        mrpt::format(
            "expected an honestly large sigma from weak data, got %.3f deg; the test no longer "
            "exercises the ungated path",
            sigmaDeg));

    // ... and the estimate itself is still good: this is exactly what a
    // quality gate would have thrown away.
    const double dirErrDeg = angleBetweenDeg(res.gravity_in_map, w.gravity_in_map());
    ASSERT_LT_(dirErrDeg, 1.0);

    std::cout << "test_weak_data_is_reported_not_gated passed (sigma " << sigmaDeg
              << " deg, dir err " << dirErrDeg << " deg)." << std::endl;
}

void test_coverage_guard_rejects_truncated_window()
{
    SyntheticWorld w;

    MapGravityEstimator est;

    // An interval whose samples only cover the first half of it: exactly the
    // failure that silently corrupts a preintegrated constraint.
    auto itv = makeInterval(w, 0.0, 0.5, 200.0, ImuIntegrationParams());
    itv.t_to = itv.t_from + 1.0;  // claim twice the duration actually integrated

    ASSERT_(!est.add_interval(itv));
    ASSERT_(!est.last_rejection_reason().empty());
    ASSERT_EQUAL_(est.window_length(), 0U);

    std::cout << "test_coverage_guard_rejects_truncated_window passed (\""
              << est.last_rejection_reason() << "\")." << std::endl;
}

void test_input_guards()
{
    SyntheticWorld w;

    MapGravityEstimator est;

    // Empty delta:
    {
        MapGravityEstimator::Interval itv;
        itv.t_from = 0;
        itv.t_to   = 0.5;
        ASSERT_(!est.add_interval(itv));
    }
    // Non-positive duration:
    {
        auto itv = makeInterval(w, 0.0, 0.5, 200.0, ImuIntegrationParams());
        itv.t_to = itv.t_from;
        ASSERT_(!est.add_interval(itv));
    }
    // Non-finite velocity:
    {
        auto itv   = makeInterval(w, 0.0, 0.5, 200.0, ImuIntegrationParams());
        itv.v_to.x = std::numeric_limits<double>::quiet_NaN();
        ASSERT_(!est.add_interval(itv));
    }
    // Out-of-order:
    {
        ASSERT_(est.add_interval(makeInterval(w, 1.0, 1.5, 200.0, ImuIntegrationParams())));
        ASSERT_(!est.add_interval(makeInterval(w, 0.0, 0.5, 200.0, ImuIntegrationParams())));
    }

    // Solving with a single interval must not crash (it is simply weak):
    ASSERT_(est.solve());
    ASSERT_(est.latest_result().has_value());
    ASSERT_(!est.latest_result()->converged);  // too few intervals

    std::cout << "test_input_guards passed." << std::endl;
}

void test_window_slides_and_resets()
{
    SyntheticWorld w;

    MapGravityEstimator est;
    est.parameters.window_size = 5;

    for (int k = 0; k < 12; k++)
    {
        const double t0 = k * 0.5;
        ASSERT_(est.add_interval(makeInterval(w, t0, t0 + 0.5, 200.0, ImuIntegrationParams())));
    }
    ASSERT_EQUAL_(est.window_length(), 5U);

    est.reset();
    ASSERT_EQUAL_(est.window_length(), 0U);
    ASSERT_(!est.latest_result().has_value());

    std::cout << "test_window_slides_and_resets passed." << std::endl;
}

void test_correction_from_gravity_degenerate_inputs()
{
    // Zero vector: identity, no crash.
    ASSERT_LT_(so3_log(MapGravityEstimator::correction_from_gravity({0, 0, 0})).norm(), 1e-12);

    // Already level: identity.
    ASSERT_LT_(
        so3_log(MapGravityEstimator::correction_from_gravity({0, 0, -GRAVITY})).norm(), 1e-12);

    // Upside down: a valid 180 deg rotation, and it must still level.
    {
        const auto            C = MapGravityEstimator::correction_from_gravity({0, 0, GRAVITY});
        const Eigen::Vector3d leveled = C.asEigen() * Eigen::Vector3d(0, 0, GRAVITY);
        ASSERT_NEAR_(leveled.z(), -GRAVITY, 1e-9);
    }

    // The correction must carry NO yaw, so a pure-roll gravity tilt yields a
    // pure-roll correction.
    {
        const double    a = mrpt::DEG2RAD(5.0);
        const TVector3D g{0, GRAVITY * std::sin(a), -GRAVITY * std::cos(a)};
        const auto      C = MapGravityEstimator::correction_from_gravity(g);

        mrpt::poses::CPose3D p;
        p.setRotationMatrix(C);
        ASSERT_NEAR_(p.yaw(), 0.0, 1e-9);
        ASSERT_NEAR_(std::abs(p.roll()), a, 1e-9);
    }

    std::cout << "test_correction_from_gravity_degenerate_inputs passed." << std::endl;
}

}  // namespace

int main([[maybe_unused]] int argc, [[maybe_unused]] char** argv)
{
    try
    {
        test_recovers_map_tilt_while_always_moving();
        test_level_map_yields_no_correction();
        test_gyro_bias_observability();
        test_accel_bias_observability_with_rotation_excitation();
        test_no_rotation_excitation_degrades_gracefully();
        test_huber_rejects_corrupted_velocity();
        test_survives_realistic_noise();
        test_residual_jacobians();
        test_weak_data_is_reported_not_gated();
        test_coverage_guard_rejects_truncated_window();
        test_input_guards();
        test_window_slides_and_resets();
        test_correction_from_gravity_degenerate_inputs();
        std::cout << "All MapGravityEstimator tests passed." << std::endl;
        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Test failed: " << e.what() << std::endl;
        return 1;
    }
}
