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
 * @file   MapGravityEstimator.cpp
 * @brief  Estimates the gravity vector in the odometry/map frame.
 * @author Jose Luis Blanco Claraco
 * @date   Jul 22, 2026
 */

#include <mola_imu_preintegration/MapGravityEstimator.h>
#include <mrpt/core/bits_math.h>
#include <mrpt/core/exceptions.h>
#include <mrpt/poses/CPose3D.h>

#include <Eigen/Dense>
#include <cmath>
#include <iostream>

using namespace mola::imu;

namespace
{
Eigen::Vector3d       toEigen(const mrpt::math::TVector3D& v) { return {v.x, v.y, v.z}; }
mrpt::math::TVector3D fromEigen(const Eigen::Vector3d& v) { return {v.x(), v.y(), v.z()}; }

bool isFinite(const mrpt::math::TVector3D& v)
{
    return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

/// Extracts (pitch, roll) from a rotation matrix.
std::pair<double, double> pitchRollOf(const SO3& R)
{
    mrpt::poses::CPose3D p;
    p.setRotationMatrix(R);
    return {p.pitch(), p.roll()};
}

/// Ridge added to the normal equations, purely for numerical safety of the
/// LDLT: small enough to be irrelevant next to any real information.
constexpr double NUMERIC_RIDGE = 1e-12;

/// Below this vector norm a direction is considered undefined.
constexpr double DIRECTION_EPS = 1e-9;
}  // namespace

// ---------------------------------------------------------------------------
// Parameters
// ---------------------------------------------------------------------------
void MapGravityEstimator::Parameters::load_from(const mrpt::containers::yaml& cfg)
{
    MCP_LOAD_OPT(cfg, window_size);
    MCP_LOAD_OPT(cfg, gravity_magnitude);
    MCP_LOAD_OPT(cfg, gravity_magnitude_sigma);
    MCP_LOAD_OPT(cfg, min_interval_coverage);
    MCP_LOAD_OPT(cfg, min_interval_duration);
    MCP_LOAD_OPT(cfg, max_interval_duration);
    MCP_LOAD_OPT(cfg, odometry_relative_rotation_sigma_deg);
    MCP_LOAD_OPT(cfg, bias_acc_prior_sigma);
    MCP_LOAD_OPT(cfg, bias_gyro_prior_sigma);
    MCP_LOAD_OPT(cfg, huber_threshold);
    MCP_LOAD_OPT(cfg, max_iterations);
    MCP_LOAD_OPT(cfg, tolerance);
    MCP_LOAD_OPT(cfg, min_intervals_for_convergence);

    // Removed parameter: warn instead of ignoring it silently, since a config
    // that still sets it expects a filtering behavior that no longer exists.
    if (cfg.has("max_tilt_sigma_deg"))
    {
        std::cerr << "[MapGravityEstimator] WARNING: parameter 'max_tilt_sigma_deg' no longer "
                     "exists and is being IGNORED. Quality is now reported through the earned "
                     "pitch_sigma/roll_sigma of the result, not filtered out internally.\n";
    }

    ASSERT_GT_(window_size, 0U);
    ASSERT_GT_(gravity_magnitude, 0.0);
    ASSERT_GT_(gravity_magnitude_sigma, 0.0);
    ASSERT_GT_(min_interval_coverage, 0.0);
    ASSERT_GT_(min_interval_duration, 0.0);
    ASSERT_GT_(max_interval_duration, min_interval_duration);
    ASSERT_GT_(odometry_relative_rotation_sigma_deg, 0.0);
    ASSERT_GT_(bias_acc_prior_sigma, 0.0);
    ASSERT_GT_(bias_gyro_prior_sigma, 0.0);
    ASSERT_GT_(huber_threshold, 0.0);
}

void MapGravityEstimator::Parameters::save_to(mrpt::containers::yaml& cfg) const
{
    MCP_SAVE(cfg, window_size);
    MCP_SAVE(cfg, gravity_magnitude);
    MCP_SAVE(cfg, gravity_magnitude_sigma);
    MCP_SAVE(cfg, min_interval_coverage);
    MCP_SAVE(cfg, min_interval_duration);
    MCP_SAVE(cfg, max_interval_duration);
    MCP_SAVE(cfg, odometry_relative_rotation_sigma_deg);
    MCP_SAVE(cfg, bias_acc_prior_sigma);
    MCP_SAVE(cfg, bias_gyro_prior_sigma);
    MCP_SAVE(cfg, huber_threshold);
    MCP_SAVE(cfg, max_iterations);
    MCP_SAVE(cfg, tolerance);
    MCP_SAVE(cfg, min_intervals_for_convergence);
}

mrpt::math::CMatrixDouble33 MapGravityEstimator::Interval::defaultVelocityCov()
{
    mrpt::math::CMatrixDouble33 c;
    c.setDiagonal(mrpt::square(0.25));  // 0.25 m/s, deliberately loose
    return c;
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------
SO3 MapGravityEstimator::correction_from_gravity(const mrpt::math::TVector3D& gravityInMap)
{
    // The "up" direction implied by the measured gravity:
    const double n = gravityInMap.norm();
    if (n < DIRECTION_EPS)
    {
        return SO3::Identity();
    }
    const Eigen::Vector3d up = -toEigen(gravityInMap) / n;
    const Eigen::Vector3d z(0, 0, 1);

    // Minimal rotation taking `up` onto +Z: axis = up x z, angle between them.
    // Being the minimal (geodesic) rotation, it has no yaw component by
    // construction, so it re-levels the frame without spinning it.
    const Eigen::Vector3d axis = up.cross(z);
    const double          s    = axis.norm();
    const double          c    = up.dot(z);

    if (s < DIRECTION_EPS)
    {
        if (c > 0)
        {
            return SO3::Identity();  // already aligned
        }
        // Antiparallel: any perpendicular axis is valid; pick X.
        return so3_exp({M_PI, 0, 0});
    }

    const double angle = std::atan2(s, c);
    return so3_exp(fromEigen(Eigen::Vector3d(axis / s * angle)));
}

// ---------------------------------------------------------------------------
// Residuals
// ---------------------------------------------------------------------------
mrpt::math::TVector3D MapGravityEstimator::gravity_residual(
    const Interval& interval, const mrpt::math::TVector3D& gravityInMap,
    const LinearAcceleration& biasAcc, const AngularVelocity& biasGyro,
    const mrpt::optional_ref<mrpt::math::CMatrixDouble33>& J_gravity,
    const mrpt::optional_ref<mrpt::math::CMatrixDouble33>& J_bias_acc,
    const mrpt::optional_ref<mrpt::math::CMatrixDouble33>& J_bias_gyro)
{
    const auto& d = interval.delta;

    const Eigen::Matrix3d Ri = interval.R_from.asEigen();
    const double          dt = d.deltaTij;

    // Bias increments with respect to the interval's own linearization point:
    const Eigen::Vector3d dba = toEigen(biasAcc - d.bias_acc);
    const Eigen::Vector3d dbg = toEigen(biasGyro - d.bias_gyro);

    // First-order bias update of the preintegrated velocity delta:
    const Eigen::Vector3d dV =
        toEigen(d.deltaVij) + d.dV_dba.asEigen() * dba + d.dV_dbg.asEigen() * dbg;

    const Eigen::Vector3d r =
        toEigen(interval.v_to) - toEigen(interval.v_from) - toEigen(gravityInMap) * dt - Ri * dV;

    if (J_gravity.has_value())
    {
        J_gravity->get().asEigen() = -Eigen::Matrix3d::Identity() * dt;
    }
    if (J_bias_acc.has_value())
    {
        J_bias_acc->get().asEigen() = -Ri * d.dV_dba.asEigen();
    }
    if (J_bias_gyro.has_value())
    {
        J_bias_gyro->get().asEigen() = -Ri * d.dV_dbg.asEigen();
    }

    return fromEigen(r);
}

mrpt::math::TVector3D MapGravityEstimator::rotation_residual(
    const Interval& interval, const AngularVelocity& biasGyro,
    const mrpt::optional_ref<mrpt::math::CMatrixDouble33>& J_bias_gyro)
{
    const auto& d = interval.delta;

    const Eigen::Matrix3d Ri = interval.R_from.asEigen();
    const Eigen::Matrix3d Rj = interval.R_to.asEigen();

    const Eigen::Vector3d dbg = toEigen(biasGyro - d.bias_gyro);

    // First-order bias update of the preintegrated rotation, on the manifold:
    const Eigen::Vector3d       Jdbg   = d.dR_dbg.asEigen() * dbg;
    const mrpt::math::TVector3D JdbgV  = fromEigen(Jdbg);
    const Eigen::Matrix3d       dRcorr = d.deltaRij.asEigen() * so3_exp(JdbgV).asEigen();

    SO3 M;
    M.asEigen() = dRcorr.transpose() * Ri.transpose() * Rj;

    const mrpt::math::TVector3D e = so3_log(M);

    if (J_bias_gyro.has_value())
    {
        const Eigen::Matrix3d JrInv = inverse_right_jacobian_so3(e).asEigen();

        J_bias_gyro->get().asEigen() = -JrInv * so3_exp(e * -1.0).asEigen() *
                                       right_jacobian_so3(JdbgV).asEigen() * d.dR_dbg.asEigen();
    }

    return e;
}

void MapGravityEstimator::set_bias_prior(
    const LinearAcceleration& biasAcc, const AngularVelocity& biasGyro)
{
    bias_acc_prior_  = biasAcc;
    bias_gyro_prior_ = biasGyro;
    if (intervals_.empty())
    {
        bias_acc_  = biasAcc;
        bias_gyro_ = biasGyro;
    }
}

void MapGravityEstimator::reset()
{
    intervals_.clear();
    latest_result_.reset();
    last_rejection_reason_.clear();
    gravity_initialized_ = false;
    bias_acc_            = bias_acc_prior_;
    bias_gyro_           = bias_gyro_prior_;
}

bool MapGravityEstimator::reject(const std::string& why)
{
    last_rejection_reason_ = why;
    return false;
}

// ---------------------------------------------------------------------------
// Window management
// ---------------------------------------------------------------------------
bool MapGravityEstimator::add_interval(const Interval& interval)
{
    const auto& d = interval.delta;

    // ---- Input guards. Rejecting is safe: the unknowns are global, so a
    //      skipped interval leaves nothing unconstrained. -------------------
    if (d.num_measurements == 0)
    {
        return reject("no IMU samples in the interval");
    }

    const double wallDt = interval.t_to - interval.t_from;
    if (!(wallDt > 0))
    {
        return reject(mrpt::format("non-positive interval duration (%g s)", wallDt));
    }
    if (wallDt < parameters.min_interval_duration)
    {
        return reject(mrpt::format(
            "interval too short (%g s < %g s)", wallDt, parameters.min_interval_duration));
    }
    if (wallDt > parameters.max_interval_duration)
    {
        return reject(mrpt::format(
            "interval too long (%g s > %g s)", wallDt, parameters.max_interval_duration));
    }
    if (!(d.deltaTij > 0))
    {
        return reject("non-positive integrated time");
    }

    // The coverage guard: a partially covered interval yields a delta that
    // does NOT span the gap, silently corrupting the constraint.
    const double coverage = d.deltaTij / wallDt;
    if (coverage < parameters.min_interval_coverage)
    {
        return reject(mrpt::format(
            "insufficient IMU coverage of the interval: %.1f%% < %.1f%% (integrated %g s of "
            "%g s, %zu samples)",
            100.0 * coverage, 100.0 * parameters.min_interval_coverage, d.deltaTij, wallDt,
            d.num_measurements));
    }

    if (!isFinite(interval.v_from) || !isFinite(interval.v_to) || !isFinite(d.deltaVij) ||
        !d.deltaRij.asEigen().allFinite() || !interval.R_from.asEigen().allFinite() ||
        !interval.R_to.asEigen().allFinite())
    {
        return reject("non-finite values in the interval");
    }

    // Consecutive intervals are expected to chain exactly (t_from == previous
    // t_to); a small tolerance keeps floating-point stamp arithmetic from
    // tripping the guard, without letting genuinely out-of-order data through.
    constexpr double STAMP_TOLERANCE = 1e-6;  // [s]
    if (!intervals_.empty() && interval.t_from < intervals_.back().t_to - STAMP_TOLERANCE)
    {
        return reject(mrpt::format(
            "out-of-order interval (t_from=%.6f < newest t_to=%.6f)", interval.t_from,
            intervals_.back().t_to));
    }

    intervals_.push_back(interval);

    while (intervals_.size() > parameters.window_size)
    {
        intervals_.pop_front();
    }

    if (!gravity_initialized_)
    {
        gravity_             = {0, 0, -parameters.gravity_magnitude};
        gravity_initialized_ = true;
    }

    last_rejection_reason_.clear();
    return true;
}

// ---------------------------------------------------------------------------
// Solver: Gauss-Newton over x = [ g (3), bias_gyro (3), bias_acc (3) ]
// ---------------------------------------------------------------------------
bool MapGravityEstimator::solve()
{
    if (intervals_.empty())
    {
        return false;
    }

    constexpr Eigen::Index IDX_G  = 0;
    constexpr Eigen::Index IDX_BG = 3;
    constexpr Eigen::Index IDX_BA = 6;
    constexpr Eigen::Index N      = 9;

    const double wBa  = 1.0 / parameters.bias_acc_prior_sigma;
    const double wBg  = 1.0 / parameters.bias_gyro_prior_sigma;
    const double wMag = 1.0 / parameters.gravity_magnitude_sigma;
    const double covOdomRot =
        mrpt::square(mrpt::DEG2RAD(parameters.odometry_relative_rotation_sigma_deg));
    const double huberK = parameters.huber_threshold;

    Eigen::Matrix<double, N, N> H              = Eigen::Matrix<double, N, N>::Zero();
    std::size_t                 iterationsDone = 0;
    double                      incrementNorm  = 0;

    for (std::size_t iter = 0; iter < parameters.max_iterations; iter++)
    {
        H.setZero();
        Eigen::Matrix<double, N, 1> b = Eigen::Matrix<double, N, 1>::Zero();

        for (const auto& itv : intervals_)
        {
            const auto& d = itv.delta;

            const Eigen::Matrix3d Ri = itv.R_from.asEigen();

            // ---- (a) gravity residual ------------------------------------
            //   r = (v_j - v_i) - g*dt - R_i * dV_ij(b)
            {
                mrpt::math::CMatrixDouble33 Jg;
                mrpt::math::CMatrixDouble33 JaM;
                mrpt::math::CMatrixDouble33 JbM;

                const Eigen::Vector3d r =
                    toEigen(gravity_residual(itv, gravity_, bias_acc_, bias_gyro_, Jg, JaM, JbM));

                Eigen::Matrix3d Jg_ = Jg.asEigen();  // d/d gravity
                Eigen::Matrix3d Ja  = JaM.asEigen();
                Eigen::Matrix3d Jb  = JbM.asEigen();

                // Weight = preintegration velocity covariance (rotated into
                // the map frame) PLUS both external velocity covariances.
                // This is the "earned" sigma: never hand-tuned.
                const Eigen::Matrix3d covDV = d.cov.asEigen().block<3, 3>(3, 3);
                const Eigen::Matrix3d cov =
                    Ri * covDV * Ri.transpose() + itv.cov_v_from.asEigen() + itv.cov_v_to.asEigen();

                const Eigen::LLT<Eigen::Matrix3d> llt(
                    cov + Eigen::Matrix3d::Identity() * NUMERIC_RIDGE);
                const auto L = llt.matrixL();

                Eigen::Vector3d rw = L.solve(r);
                Jg_                = L.solve(Jg_);
                Ja                 = L.solve(Ja);
                Jb                 = L.solve(Jb);

                // Huber: a bad odometry velocity or a jolt cannot drag the
                // estimate.
                const double rn = rw.norm();
                if (rn > huberK)
                {
                    const double sqrtW = std::sqrt(huberK / rn);
                    rw *= sqrtW;
                    Jg_ *= sqrtW;
                    Ja *= sqrtW;
                    Jb *= sqrtW;
                }

                H.block<3, 3>(IDX_G, IDX_G) += Jg_.transpose() * Jg_;
                H.block<3, 3>(IDX_G, IDX_BA) += Jg_.transpose() * Ja;
                H.block<3, 3>(IDX_BA, IDX_G) += Ja.transpose() * Jg_;
                H.block<3, 3>(IDX_G, IDX_BG) += Jg_.transpose() * Jb;
                H.block<3, 3>(IDX_BG, IDX_G) += Jb.transpose() * Jg_;

                H.block<3, 3>(IDX_BA, IDX_BA) += Ja.transpose() * Ja;
                H.block<3, 3>(IDX_BG, IDX_BG) += Jb.transpose() * Jb;
                H.block<3, 3>(IDX_BA, IDX_BG) += Ja.transpose() * Jb;
                H.block<3, 3>(IDX_BG, IDX_BA) += Jb.transpose() * Ja;

                b.segment<3>(IDX_G) += Jg_.transpose() * rw;
                b.segment<3>(IDX_BA) += Ja.transpose() * rw;
                b.segment<3>(IDX_BG) += Jb.transpose() * rw;
            }

            // ---- (b) relative-rotation residual (gyro bias) --------------
            //   r = Log( dR_ij(b_g)^T * R_i^T * R_j )
            {
                mrpt::math::CMatrixDouble33 JbM;

                const mrpt::math::TVector3D e = rotation_residual(itv, bias_gyro_, JbM);

                Eigen::Matrix3d Jb = JbM.asEigen();

                // Weight: preintegrated rotation covariance PLUS the assumed
                // uncertainty of the external odometry's relative rotation.
                const Eigen::Matrix3d cov =
                    d.cov.asEigen().block<3, 3>(0, 0) + Eigen::Matrix3d::Identity() * covOdomRot;

                const Eigen::LLT<Eigen::Matrix3d> llt(
                    cov + Eigen::Matrix3d::Identity() * NUMERIC_RIDGE);
                const auto L = llt.matrixL();

                const Eigen::Vector3d rw = L.solve(toEigen(e));
                Jb                       = L.solve(Jb);

                H.block<3, 3>(IDX_BG, IDX_BG) += Jb.transpose() * Jb;
                b.segment<3>(IDX_BG) += Jb.transpose() * rw;
            }
        }

        // ---- Gravity magnitude soft constraint ---------------------------
        {
            const double gn = gravity_.norm();
            if (gn > DIRECTION_EPS)
            {
                const double             r = wMag * (gn - parameters.gravity_magnitude);
                const Eigen::RowVector3d J = wMag * (toEigen(gravity_) / gn).transpose();

                H.block<3, 3>(IDX_G, IDX_G) += J.transpose() * J;
                b.segment<3>(IDX_G) += J.transpose() * r;
            }
        }

        // ---- Zero-anchored bias priors -----------------------------------
        {
            const Eigen::Vector3d rg = wBg * toEigen(bias_gyro_ - bias_gyro_prior_);
            const Eigen::Vector3d ra = wBa * toEigen(bias_acc_ - bias_acc_prior_);

            H.block<3, 3>(IDX_BG, IDX_BG) += Eigen::Matrix3d::Identity() * (wBg * wBg);
            H.block<3, 3>(IDX_BA, IDX_BA) += Eigen::Matrix3d::Identity() * (wBa * wBa);

            b.segment<3>(IDX_BG) += wBg * rg;
            b.segment<3>(IDX_BA) += wBa * ra;
        }

        // ---- Solve and update --------------------------------------------
        H.diagonal().array() += NUMERIC_RIDGE;

        const Eigen::LDLT<Eigen::Matrix<double, N, N>> ldlt(H);
        if (ldlt.info() != Eigen::Success)
        {
            return false;
        }
        const Eigen::Matrix<double, N, 1> dx = ldlt.solve(-b);
        if (!dx.allFinite())
        {
            return false;
        }

        gravity_   = gravity_ + fromEigen(dx.segment<3>(IDX_G));
        bias_gyro_ = bias_gyro_ + fromEigen(dx.segment<3>(IDX_BG));
        bias_acc_  = bias_acc_ + fromEigen(dx.segment<3>(IDX_BA));

        iterationsDone = iter + 1;
        incrementNorm  = dx.norm();

        if (incrementNorm < parameters.tolerance)
        {
            break;
        }
    }

    // ---- Build the result -------------------------------------------------
    Result res;
    res.timestamp      = intervals_.back().t_to;
    res.gravity_in_map = gravity_;
    res.correction     = correction_from_gravity(gravity_);
    res.bias_acc       = bias_acc_;
    res.bias_gyro      = bias_gyro_;
    res.num_intervals  = intervals_.size();
    res.iterations     = iterationsDone;
    res.increment_norm = incrementNorm;

    std::tie(res.pitch_correction, res.roll_correction) = pitchRollOf(res.correction);
    res.tilt                                            = so3_log(res.correction).norm();

    // "Earned" sigmas: the marginal covariance of the gravity block, pushed
    // through the numerical Jacobian of (pitch, roll) of the correction.
    bool covOk = false;
    {
        const Eigen::LDLT<Eigen::Matrix<double, N, N>> ldlt(H);
        if (ldlt.info() == Eigen::Success)
        {
            Eigen::Matrix<double, N, 3> sel = Eigen::Matrix<double, N, 3>::Zero();
            sel.block<3, 3>(IDX_G, 0)       = Eigen::Matrix3d::Identity();

            const Eigen::Matrix<double, N, 3> covCols = ldlt.solve(sel);
            const Eigen::Matrix3d             covG    = covCols.block<3, 3>(IDX_G, 0);

            if (covG.allFinite())
            {
                Eigen::Matrix<double, 2, 3> J   = Eigen::Matrix<double, 2, 3>::Zero();
                const double                eps = 1e-6;
                for (int i = 0; i < 3; i++)
                {
                    mrpt::math::TVector3D dg{0, 0, 0};
                    dg[i] = eps;

                    const auto [pp, rp] = pitchRollOf(correction_from_gravity(gravity_ + dg));
                    const auto [pm, rm] = pitchRollOf(correction_from_gravity(gravity_ - dg));

                    J(0, i) = (pp - pm) / (2 * eps);
                    J(1, i) = (rp - rm) / (2 * eps);
                }

                const Eigen::Matrix2d covPR = J * covG * J.transpose();
                if (covPR(0, 0) >= 0 && covPR(1, 1) >= 0)
                {
                    res.pitch_sigma = std::sqrt(covPR(0, 0));
                    res.roll_sigma  = std::sqrt(covPR(1, 1));
                    covOk           = true;
                }
            }
        }
    }

    // "Converged" means USABLE: enough data, and finite earned sigmas. There
    // is deliberately no accuracy threshold here; a large sigma is an honest
    // weak estimate that the consumer can still weight correctly, and hiding
    // it behind a boolean is what silently disables the estimator on real
    // hardware.
    res.converged = covOk && res.num_intervals >= parameters.min_intervals_for_convergence;

    latest_result_ = res;
    return true;
}
