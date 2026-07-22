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
 * @file   test-imu-preintegrator.cpp
 * @brief  Unit tests for full on-manifold IMU preintegration.
 * @author Jose Luis Blanco Claraco
 * @date   Jul 22, 2026
 *
 * The reference is a synthetic trajectory whose true pose/velocity is known in
 * closed form. The IMU readings are SYNTHESIZED from it, so the preintegrated
 * deltas are compared against the exact analytic answer, not against a second
 * implementation of the same equations.
 */

#include <mola_imu_preintegration/ImuPreintegrator.h>
#include <mrpt/core/exceptions.h>

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

using namespace mola::imu;
using mrpt::math::TVector3D;

namespace
{
constexpr double GRAVITY = 9.81;

/// A synthetic ground-truth trajectory: position, velocity and attitude as
/// analytic functions of time, plus the IMU readings they imply.
struct SyntheticTrajectory
{
    TVector3D omega{0.35, -0.20, 0.55};  //!< constant body angular rate [rad/s]
    TVector3D tilt0{0.08, -0.05, 0.0};  //!< initial attitude, as a rotation vector

    /// Position follows a Lissajous curve, so acceleration is nontrivial in
    /// all three axes.
    TVector3D amplitude{2.0, 1.5, 0.6};
    TVector3D frequency{0.7, 1.1, 1.7};  //!< [rad/s]

    SO3 R(double t) const
    {
        SO3 out;
        out.asEigen() = so3_exp(tilt0).asEigen() * so3_exp(omega * t).asEigen();
        return out;
    }

    TVector3D p(double t) const
    {
        return {
            amplitude.x * std::sin(frequency.x * t), amplitude.y * std::sin(frequency.y * t),
            amplitude.z * std::sin(frequency.z * t)};
    }

    TVector3D v(double t) const
    {
        return {
            amplitude.x * frequency.x * std::cos(frequency.x * t),
            amplitude.y * frequency.y * std::cos(frequency.y * t),
            amplitude.z * frequency.z * std::cos(frequency.z * t)};
    }

    /// World-frame coordinate acceleration.
    TVector3D a_world(double t) const
    {
        return {
            -amplitude.x * frequency.x * frequency.x * std::sin(frequency.x * t),
            -amplitude.y * frequency.y * frequency.y * std::sin(frequency.y * t),
            -amplitude.z * frequency.z * frequency.z * std::sin(frequency.z * t)};
    }

    /// What an ideal accelerometer reads: the proper acceleration (specific
    /// force) in the body frame, i.e. R^T (a_world - g).
    TVector3D accel_reading(double t) const
    {
        const TVector3D       aw = a_world(t);
        const Eigen::Vector3d d(aw.x, aw.y, aw.z + GRAVITY);
        const Eigen::Vector3d b = R(t).asEigen().transpose() * d;
        return {b.x(), b.y(), b.z()};
    }

    /// What an ideal gyroscope reads: the body-frame angular rate (constant).
    TVector3D gyro_reading(double /*t*/) const { return omega; }
};

/// Feeds the trajectory's readings into a preintegrator over [t0, t1].
void integrate(
    ImuPreintegrator& pim, const SyntheticTrajectory& traj, double t0, double t1, double dt,
    const TVector3D& biasAcc = {0, 0, 0}, const TVector3D& biasGyro = {0, 0, 0})
{
    const int n = static_cast<int>(std::llround((t1 - t0) / dt));
    for (int k = 0; k < n; k++)
    {
        // Sample at the START of each step, matching the integrator's own
        // convention (a reading is held constant over the FOLLOWING dt).
        // Sampling elsewhere would measure the mismatch between the two
        // conventions rather than the integrator's own accuracy.
        const double t = t0 + k * dt;
        pim.integrate_measurement(
            traj.accel_reading(t) + biasAcc, traj.gyro_reading(t) + biasGyro, dt);
    }
}

struct ExactDeltas
{
    SO3       dR;
    TVector3D dV;
    TVector3D dP;
};

/// The exact deltas implied by the trajectory definition.
ExactDeltas exactDeltas(const SyntheticTrajectory& traj, double t0, double t1)
{
    const double    dt = t1 - t0;
    const TVector3D g{0, 0, -GRAVITY};

    const Eigen::Matrix3d Ri = traj.R(t0).asEigen();

    ExactDeltas out;
    out.dR.asEigen() = Ri.transpose() * traj.R(t1).asEigen();

    const TVector3D vi = traj.v(t0);
    const TVector3D dv = traj.v(t1) - vi - g * dt;
    const TVector3D dp = traj.p(t1) - traj.p(t0) - vi * dt - g * (0.5 * dt * dt);

    const Eigen::Vector3d dve = Ri.transpose() * Eigen::Vector3d(dv.x, dv.y, dv.z);
    const Eigen::Vector3d dpe = Ri.transpose() * Eigen::Vector3d(dp.x, dp.y, dp.z);

    out.dV = {dve.x(), dve.y(), dve.z()};
    out.dP = {dpe.x(), dpe.y(), dpe.z()};
    return out;
}

double rotationErrorDeg(const SO3& a, const SO3& b)
{
    SO3 d;
    d.asEigen() = a.asEigen().transpose() * b.asEigen();
    return mrpt::RAD2DEG(so3_log(d).norm());
}

// -----------------------------------------------------------------------
// The gold-standard test: deltas must match the analytic answer.
// -----------------------------------------------------------------------
void test_matches_analytic_trajectory()
{
    const SyntheticTrajectory traj;
    const double              t0 = 0.0;
    const double              t1 = 2.5;

    ImuPreintegrator pim;
    pim.reset_integration();
    integrate(pim, traj, t0, t1, 1.0 / 400.0);

    const auto& st = pim.current_state();
    const auto  ex = exactDeltas(traj, t0, t1);

    ASSERT_NEAR_(st.deltaTij, t1 - t0, 1e-9);
    ASSERT_GT_(st.num_measurements, 900U);

    // Rotation is exact here regardless of dt: the gyro rate is constant, so
    // composing Exp(w*dt) reproduces Exp(w*T) exactly.
    const double rotErr = rotationErrorDeg(st.deltaRij, ex.dR);
    ASSERTMSG_(rotErr < 1e-3, mrpt::format("deltaR error too large: %g deg", rotErr));

    // Velocity/position carry the first-order discretization error inherent to
    // the standard Euler-forward preintegration form (see the class docs). At
    // 400 Hz over 2.5 s of aggressive motion that is a few cm/s; the
    // convergence test below is what pins down the error ORDER.
    const double velErr = (st.deltaVij - ex.dV).norm();
    const double posErr = (st.deltaPij - ex.dP).norm();
    ASSERTMSG_(velErr < 0.05, mrpt::format("deltaV error too large: %g m/s", velErr));
    ASSERTMSG_(posErr < 0.05, mrpt::format("deltaP error too large: %g m", posErr));

    std::cout << "test_matches_analytic_trajectory passed (rot err " << rotErr << " deg)."
              << std::endl;
}

// Discretization error must shrink as the sample rate grows: this catches an
// integration scheme that is biased rather than merely noisy.
void test_converges_with_sample_rate()
{
    const SyntheticTrajectory traj;
    const auto                ex = exactDeltas(traj, 0.0, 1.0);

    // Each rate is 4x the previous one. A first-order scheme must therefore
    // shrink the error by ~4x per step: asserting the RATIO (not just "it got
    // smaller") is what actually pins down the convergence order and would
    // catch a scheme that is accidentally zeroth-order in some term.
    double prevErr = std::numeric_limits<double>::max();
    for (const double rate : {50.0, 200.0, 800.0})
    {
        ImuPreintegrator pim;
        pim.reset_integration();
        integrate(pim, traj, 0.0, 1.0, 1.0 / rate);

        const double err = (pim.current_state().deltaVij - ex.dV).norm();
        if (prevErr != std::numeric_limits<double>::max())
        {
            const double ratio = prevErr / err;
            ASSERTMSG_(
                ratio > 3.0,
                mrpt::format(
                    "convergence slower than first order at rate=%g Hz: ratio=%g (%g -> %g)", rate,
                    ratio, prevErr, err));
        }
        prevErr = err;
    }
    ASSERT_LT_(prevErr, 0.02);

    std::cout << "test_converges_with_sample_rate passed (err at 800 Hz: " << prevErr << ")."
              << std::endl;
}

// -----------------------------------------------------------------------
// First-order bias update vs. full re-integration at the perturbed bias.
// -----------------------------------------------------------------------
void test_bias_correction_matches_reintegration()
{
    const SyntheticTrajectory traj;
    const double              t0 = 0.0;
    const double              t1 = 1.5;
    const double              dt = 1.0 / 400.0;

    // A realistic bias perturbation away from the linearization point:
    const TVector3D dba{0.02, -0.03, 0.01};  // m/s^2
    const TVector3D dbg{0.002, 0.001, -0.003};  // rad/s

    // (a) integrate at zero bias, then apply the first-order correction:
    ImuPreintegrator pimA;
    pimA.reset_integration({0, 0, 0}, {0, 0, 0});
    integrate(pimA, traj, t0, t1, dt);
    const auto corrected = pimA.current_state().bias_corrected(dba, dbg);

    // (b) re-integrate from scratch with the perturbed bias as linearization
    //     point (the readings are identical; only the subtracted bias differs):
    ImuPreintegrator pimB;
    pimB.reset_integration(dba, dbg);
    integrate(pimB, traj, t0, t1, dt);
    const auto& exact = pimB.current_state();

    const double uncorrectedRotErr =
        rotationErrorDeg(pimA.current_state().deltaRij, exact.deltaRij);
    const double correctedRotErr = rotationErrorDeg(corrected.deltaRij, exact.deltaRij);

    ASSERT_GT_(uncorrectedRotErr, 0.1);
    ASSERTMSG_(
        correctedRotErr < 0.02 * uncorrectedRotErr,
        mrpt::format(
            "first-order rotation correction too weak: uncorrected=%g deg corrected=%g deg",
            uncorrectedRotErr, correctedRotErr));

    const double uncorrectedVelErr = (pimA.current_state().deltaVij - exact.deltaVij).norm();
    const double correctedVelErr   = (corrected.deltaVij - exact.deltaVij).norm();

    ASSERT_GT_(uncorrectedVelErr, 0.01);
    ASSERTMSG_(
        correctedVelErr < 0.05 * uncorrectedVelErr,
        mrpt::format(
            "first-order velocity correction too weak: uncorrected=%g corrected=%g",
            uncorrectedVelErr, correctedVelErr));

    std::cout << "test_bias_correction_matches_reintegration passed (rot " << uncorrectedRotErr
              << " -> " << correctedRotErr << " deg)." << std::endl;
}

// -----------------------------------------------------------------------
// Bias Jacobians vs. numerical differentiation of a full re-integration.
// -----------------------------------------------------------------------
void test_bias_jacobians_match_numerical()
{
    const SyntheticTrajectory traj;
    const double              t0  = 0.0;
    const double              t1  = 1.0;
    const double              dt  = 1.0 / 400.0;
    const double              eps = 1e-6;

    ImuPreintegrator pim;
    pim.reset_integration();
    integrate(pim, traj, t0, t1, dt);
    const auto& st = pim.current_state();

    const auto reintegrate = [&](const TVector3D& ba, const TVector3D& bg)
    {
        ImuPreintegrator p;
        p.reset_integration(ba, bg);
        integrate(p, traj, t0, t1, dt);
        return p.current_state();
    };

    for (int i = 0; i < 3; i++)
    {
        TVector3D db{0, 0, 0};
        db[i] = eps;

        // --- accel bias: affects dV and dP only ---
        {
            const auto plus  = reintegrate(db, {0, 0, 0});
            const auto minus = reintegrate(db * -1.0, {0, 0, 0});

            const TVector3D numV = (plus.deltaVij - minus.deltaVij) * (1.0 / (2 * eps));
            const TVector3D numP = (plus.deltaPij - minus.deltaPij) * (1.0 / (2 * eps));

            for (int r = 0; r < 3; r++)
            {
                ASSERT_NEAR_(st.dV_dba(r, i), numV[r], 1e-6);
                ASSERT_NEAR_(st.dP_dba(r, i), numP[r], 1e-6);
            }
        }

        // --- gyro bias: affects dR, dV and dP ---
        {
            const auto plus  = reintegrate({0, 0, 0}, db);
            const auto minus = reintegrate({0, 0, 0}, db * -1.0);

            // dR_dbg is defined on the manifold: dR(b) = dR_bar * Exp(J * db),
            // so its numerical counterpart is Log(dR_bar^T dR(b)) / db.
            SO3 relPlus;
            SO3 relMinus;
            relPlus.asEigen()  = st.deltaRij.asEigen().transpose() * plus.deltaRij.asEigen();
            relMinus.asEigen() = st.deltaRij.asEigen().transpose() * minus.deltaRij.asEigen();

            const TVector3D numR = (so3_log(relPlus) - so3_log(relMinus)) * (1.0 / (2 * eps));
            const TVector3D numV = (plus.deltaVij - minus.deltaVij) * (1.0 / (2 * eps));
            const TVector3D numP = (plus.deltaPij - minus.deltaPij) * (1.0 / (2 * eps));

            for (int r = 0; r < 3; r++)
            {
                ASSERT_NEAR_(st.dR_dbg(r, i), numR[r], 1e-5);
                ASSERT_NEAR_(st.dV_dbg(r, i), numV[r], 1e-4);
                ASSERT_NEAR_(st.dP_dbg(r, i), numP[r], 1e-4);
            }
        }
    }

    std::cout << "test_bias_jacobians_match_numerical passed." << std::endl;
}

// -----------------------------------------------------------------------
// Covariance propagation vs. Monte Carlo.
// -----------------------------------------------------------------------
void test_covariance_matches_monte_carlo()
{
    const SyntheticTrajectory traj;
    const double              t0 = 0.0;
    const double              t1 = 0.5;
    const double              dt = 1.0 / 200.0;

    // Noise densities (sigma units: rad/s/sqrt(Hz) and m/s^2/sqrt(Hz)):
    const double sigmaG = 5e-3;
    const double sigmaA = 5e-2;

    ImuIntegrationParams params;
    params.cov_gyro.setDiagonal(sigmaG * sigmaG);
    params.cov_acc.setDiagonal(sigmaA * sigmaA);

    // Analytic propagation:
    ImuPreintegrator pim(params);
    pim.reset_integration();
    integrate(pim, traj, t0, t1, dt);
    const auto& st = pim.current_state();

    // Monte Carlo: the per-sample discrete noise stddev is density/sqrt(dt).
    const double                     sG = sigmaG / std::sqrt(dt);
    const double                     sA = sigmaA / std::sqrt(dt);
    std::mt19937                     rng(42);
    std::normal_distribution<double> nd(0.0, 1.0);

    const int       N         = 400;
    Eigen::MatrixXd empirical = Eigen::MatrixXd::Zero(9, 9);

    for (int trial = 0; trial < N; trial++)
    {
        ImuPreintegrator p(params);
        p.reset_integration();

        const int n = static_cast<int>(std::llround((t1 - t0) / dt));
        for (int k = 0; k < n; k++)
        {
            const double    t = t0 + (k + 0.5) * dt;
            const TVector3D na{nd(rng) * sA, nd(rng) * sA, nd(rng) * sA};
            const TVector3D ng{nd(rng) * sG, nd(rng) * sG, nd(rng) * sG};
            p.integrate_measurement(traj.accel_reading(t) + na, traj.gyro_reading(t) + ng, dt);
        }

        // Error w.r.t. the noise-free result, in the same [theta, v, p] order.
        SO3 dRerr;
        dRerr.asEigen() = st.deltaRij.asEigen().transpose() * p.current_state().deltaRij.asEigen();
        const TVector3D th = so3_log(dRerr);
        const TVector3D dv = p.current_state().deltaVij - st.deltaVij;
        const TVector3D dp = p.current_state().deltaPij - st.deltaPij;

        Eigen::VectorXd e(9);
        e << th.x, th.y, th.z, dv.x, dv.y, dv.z, dp.x, dp.y, dp.z;
        empirical += e * e.transpose();
    }
    empirical /= static_cast<double>(N);

    // Compare the diagonal standard deviations. With N=400 the sampling error
    // of a stddev is ~1/sqrt(2N) ~ 3.5%, so a 25% band is a meaningful test
    // (it catches wrong scaling, wrong dt powers or a missing term) while
    // staying robust to Monte-Carlo noise.
    for (int i = 0; i < 9; i++)
    {
        const double sAna = std::sqrt(st.cov(i, i));
        const double sEmp = std::sqrt(empirical(i, i));
        ASSERT_GT_(sAna, 0.0);
        ASSERTMSG_(
            std::abs(sEmp / sAna - 1.0) < 0.25,
            mrpt::format(
                "covariance mismatch at index %i: analytic sigma=%g empirical=%g", i, sAna, sEmp));
    }

    std::cout << "test_covariance_matches_monte_carlo passed." << std::endl;
}

// -----------------------------------------------------------------------
// Degenerate / guard cases.
// -----------------------------------------------------------------------
void test_rejects_non_positive_dt()
{
    ImuPreintegrator pim;
    pim.reset_integration();

    bool threwOnZero = false;
    try
    {
        pim.integrate_measurement({0, 0, GRAVITY}, {0, 0, 0}, 0.0);
    }
    catch (const std::exception&)
    {
        threwOnZero = true;
    }
    ASSERT_(threwOnZero);

    bool threwOnNegative = false;
    try
    {
        pim.integrate_measurement({0, 0, GRAVITY}, {0, 0, 0}, -0.01);
    }
    catch (const std::exception&)
    {
        threwOnNegative = true;
    }
    ASSERT_(threwOnNegative);

    std::cout << "test_rejects_non_positive_dt passed." << std::endl;
}

void test_reset_clears_state()
{
    ImuPreintegrator pim;
    pim.reset_integration();
    pim.integrate_measurement({0.1, 0.2, 9.7}, {0.01, 0.02, 0.03}, 0.01);
    ASSERT_GT_(pim.current_state().deltaTij, 0);

    pim.reset_integration();
    ASSERT_EQUAL_(pim.current_state().deltaTij, 0);
    ASSERT_EQUAL_(pim.current_state().num_measurements, 0U);
    ASSERT_LT_(pim.current_state().deltaVij.norm(), 1e-15);
    ASSERT_LT_(pim.current_state().cov.asEigen().cwiseAbs().maxCoeff(), 1e-15);

    std::cout << "test_reset_clears_state passed." << std::endl;
}

void test_stationary_level_sensor_reads_gravity_only()
{
    // At rest and level, the accelerometer reads (0,0,+g). The preintegrated
    // dV must then be exactly +g*dt in the body frame, since the delta
    // definition has gravity subtracted OUTSIDE the integration.
    ImuPreintegrator pim;
    pim.reset_integration();

    const double dt = 0.005;
    for (int i = 0; i < 200; i++)
    {
        pim.integrate_measurement({0, 0, GRAVITY}, {0, 0, 0}, dt);
    }

    const auto& st = pim.current_state();
    ASSERT_NEAR_(st.deltaTij, 1.0, 1e-12);
    ASSERT_LT_(rotationErrorDeg(st.deltaRij, SO3::Identity()), 1e-12);
    ASSERT_NEAR_(st.deltaVij.x, 0.0, 1e-12);
    ASSERT_NEAR_(st.deltaVij.y, 0.0, 1e-12);
    ASSERT_NEAR_(st.deltaVij.z, GRAVITY, 1e-12);

    std::cout << "test_stationary_level_sensor_reads_gravity_only passed." << std::endl;
}

}  // namespace

int main([[maybe_unused]] int argc, [[maybe_unused]] char** argv)
{
    try
    {
        test_matches_analytic_trajectory();
        test_converges_with_sample_rate();
        test_bias_correction_matches_reintegration();
        test_bias_jacobians_match_numerical();
        test_covariance_matches_monte_carlo();
        test_rejects_non_positive_dt();
        test_reset_clears_state();
        test_stationary_level_sensor_reads_gravity_only();
        std::cout << "All ImuPreintegrator tests passed." << std::endl;
        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Test failed: " << e.what() << std::endl;
        return 1;
    }
}
