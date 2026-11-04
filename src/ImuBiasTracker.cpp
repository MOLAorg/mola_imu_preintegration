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

#include <mola_imu_preintegration/ImuBiasTracker.h>
#include <mrpt/core/bits_math.h>
#include <mrpt/core/exceptions.h>
#include <mrpt/system/datetime.h>

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>

using namespace mola::imu;

namespace
{
using Mat33 = mrpt::math::CMatrixDouble33;

Mat33 diag(double sigma)
{
    Mat33 m;
    m.setDiagonal(sigma * sigma);
    return m;
}

Eigen::Vector3d toEigen(const mrpt::math::TVector3D& v) { return {v.x, v.y, v.z}; }

mrpt::math::TVector3D fromEigen(const Eigen::Vector3d& v) { return {v.x(), v.y(), v.z()}; }

/// Kalman update of a 3D state with a full 3D measurement of it.
void kalman_update(
    mrpt::math::TVector3D& x, Mat33& P, const mrpt::math::TVector3D& z, const Mat33& R)
{
    const Eigen::Matrix3d S  = P.asEigen() + R.asEigen();
    const Eigen::Matrix3d K  = P.asEigen() * S.inverse();
    const Eigen::Vector3d xn = toEigen(x) + K * (toEigen(z) - toEigen(x));
    Eigen::Matrix3d       Pn = (Eigen::Matrix3d::Identity() - K) * P.asEigen();
    Pn                       = 0.5 * (Pn + Pn.transpose()).eval();  // keep it symmetric
    x                        = fromEigen(xn);
    P                        = Mat33(Pn);
}

/// Kalman update of a 3D state with a scalar measurement z = h^T x, variance r.
void kalman_update_scalar(
    mrpt::math::TVector3D& x, Mat33& P, const Eigen::Vector3d& h, double z, double r)
{
    const Eigen::Vector3d Ph = P.asEigen() * h;
    const double          s  = h.dot(Ph) + r;
    const Eigen::Vector3d K  = Ph / s;
    const Eigen::Vector3d xn = toEigen(x) + K * (z - h.dot(toEigen(x)));
    Eigen::Matrix3d       Pn = P.asEigen() - K * Ph.transpose();
    Pn                       = 0.5 * (Pn + Pn.transpose()).eval();
    x                        = fromEigen(xn);
    P                        = Mat33(Pn);
}

}  // namespace

void ImuBiasTracker::Parameters::load_from(const mrpt::containers::yaml& cfg)
{
    MCP_LOAD_OPT(cfg, standstill_window);
    MCP_LOAD_OPT(cfg, standstill_min_samples);
    MCP_LOAD_OPT(cfg, max_sample_gap);
    MCP_LOAD_OPT(cfg, standstill_max_gyro_rms);
    MCP_LOAD_OPT(cfg, standstill_max_gyro_mean_offset);
    MCP_LOAD_OPT(cfg, standstill_max_acc_rms);
    MCP_LOAD_OPT(cfg, standstill_max_acc_norm_error);
    MCP_LOAD_OPT(cfg, gravity);
    MCP_LOAD_OPT(cfg, require_static_hint);
    MCP_LOAD_OPT(cfg, gyro_bias_random_walk);
    MCP_LOAD_OPT(cfg, acc_bias_random_walk);
    MCP_LOAD_OPT(cfg, gyro_bias_prior_sigma);
    MCP_LOAD_OPT(cfg, acc_bias_prior_sigma);
    MCP_LOAD_OPT(cfg, gyro_measurement_min_sigma);
    MCP_LOAD_OPT(cfg, acc_measurement_min_sigma);
}

void ImuBiasTracker::Parameters::save_to(mrpt::containers::yaml& cfg) const
{
    MCP_SAVE(cfg, standstill_window);
    MCP_SAVE(cfg, standstill_min_samples);
    MCP_SAVE(cfg, max_sample_gap);
    MCP_SAVE(cfg, standstill_max_gyro_rms);
    MCP_SAVE(cfg, standstill_max_gyro_mean_offset);
    MCP_SAVE(cfg, standstill_max_acc_rms);
    MCP_SAVE(cfg, standstill_max_acc_norm_error);
    MCP_SAVE(cfg, gravity);
    MCP_SAVE(cfg, require_static_hint);
    MCP_SAVE(cfg, gyro_bias_random_walk);
    MCP_SAVE(cfg, acc_bias_random_walk);
    MCP_SAVE(cfg, gyro_bias_prior_sigma);
    MCP_SAVE(cfg, acc_bias_prior_sigma);
    MCP_SAVE(cfg, gyro_measurement_min_sigma);
    MCP_SAVE(cfg, acc_measurement_min_sigma);
}

ImuBiasTracker::PerImu& ImuBiasTracker::imu_state(const std::string& sensorLabel, TimeStamp t)
{
    auto it = imus_.find(sensorLabel);
    if (it == imus_.end())
    {
        PerImu s;
        s.est.cov_gyro = diag(parameters.gyro_bias_prior_sigma);
        s.est.cov_acc  = diag(parameters.acc_bias_prior_sigma);
        s.est.stamp    = t;
        it             = imus_.emplace(sensorLabel, std::move(s)).first;
    }
    return it->second;
}

void ImuBiasTracker::predict(Estimate& e, TimeStamp t) const
{
    const double dt = t - e.stamp;
    if (dt <= 0)
    {
        return;
    }
    for (int i = 0; i < 3; i++)
    {
        e.cov_gyro(i, i) += mrpt::square(parameters.gyro_bias_random_walk) * dt;
        e.cov_acc(i, i) += mrpt::square(parameters.acc_bias_random_walk) * dt;
    }
    e.stamp = t;
}

bool ImuBiasTracker::add(const mrpt::obs::CObservationIMU& imu)
{
    using namespace mrpt::obs;
    if (!imu.has(IMU_WX) || !imu.has(IMU_WY) || !imu.has(IMU_WZ))
    {
        return false;
    }
    const TimeStamp t = mrpt::Clock::toDouble(imu.timestamp);
    auto&           s = imu_state(imu.sensorLabel, t);

    Sample smp;
    smp.t = t;
    smp.w = {imu.get(IMU_WX), imu.get(IMU_WY), imu.get(IMU_WZ)};
    if (imu.has(IMU_X_ACC) && imu.has(IMU_Y_ACC) && imu.has(IMU_Z_ACC))
    {
        smp.a = mrpt::math::TVector3D(imu.get(IMU_X_ACC), imu.get(IMU_Y_ACC), imu.get(IMU_Z_ACC));
    }
    if (!std::isfinite(smp.w.x) || !std::isfinite(smp.w.y) || !std::isfinite(smp.w.z))
    {
        return false;
    }

    // Out of order or a gap: start a new window.
    if (!s.window.empty() &&
        (t <= s.window.back().t || t - s.window.back().t > parameters.max_sample_gap))
    {
        s.window.clear();
    }
    s.window.push_back(smp);

    // Hints older than a few windows are no longer needed:
    const TimeStamp oldest = t - 4 * parameters.standstill_window;
    while (!hints_.empty() && hints_.begin()->first < oldest)
    {
        hints_.erase(hints_.begin());
    }

    if (s.window.back().t - s.window.front().t < parameters.standstill_window)
    {
        return false;
    }
    const bool updated = process_window(s);
    s.window.clear();
    return updated;
}

bool ImuBiasTracker::hints_allow(TimeStamp from, TimeStamp to) const
{
    bool anyStatic = false;
    for (auto it = hints_.lower_bound(from); it != hints_.end() && it->first <= to; ++it)
    {
        if (!it->second)
        {
            return false;  // the vehicle was moving at some point
        }
        anyStatic = true;
    }
    return anyStatic || !parameters.require_static_hint;
}

bool ImuBiasTracker::process_window(PerImu& s)
{
    const auto& p = parameters;
    const auto  n = s.window.size();
    if (n < p.standstill_min_samples || n < 2)
    {
        return false;
    }
    const TimeStamp t0 = s.window.front().t;
    const TimeStamp t1 = s.window.back().t;
    if (!hints_allow(t0, t1))
    {
        return false;
    }

    // Gyroscope: mean and per-axis variance.
    Eigen::Vector3d wMean = Eigen::Vector3d::Zero();
    for (const auto& smp : s.window)
    {
        wMean += toEigen(smp.w);
    }
    wMean /= static_cast<double>(n);
    Eigen::Vector3d wVar = Eigen::Vector3d::Zero();
    for (const auto& smp : s.window)
    {
        wVar += (toEigen(smp.w) - wMean).cwiseAbs2();
    }
    wVar /= static_cast<double>(n - 1);
    if (std::sqrt(wVar.sum()) > p.standstill_max_gyro_rms)
    {
        return false;  // vibrating or rotating: moving
    }
    if ((wMean - toEigen(s.est.bias_gyro)).norm() > p.standstill_max_gyro_mean_offset)
    {
        return false;  // a steady rotation, not a bias
    }

    // Accelerometer (optional): mean, spread and norm.
    std::size_t     na    = 0;
    Eigen::Vector3d aMean = Eigen::Vector3d::Zero();
    for (const auto& smp : s.window)
    {
        if (smp.a)
        {
            aMean += toEigen(*smp.a);
            na++;
        }
    }
    std::optional<double> accAlongVar;
    Eigen::Vector3d       u = Eigen::Vector3d::UnitZ();
    if (na >= 2)
    {
        aMean /= static_cast<double>(na);
        double aVar     = 0;
        double alongVar = 0;
        u               = aMean.normalized();
        for (const auto& smp : s.window)
        {
            if (smp.a)
            {
                const Eigen::Vector3d d = toEigen(*smp.a) - aMean;
                aVar += d.squaredNorm();
                alongVar += mrpt::square(d.dot(u));
            }
        }
        aVar /= static_cast<double>(na - 1);
        alongVar /= static_cast<double>(na - 1);
        if (std::sqrt(aVar) > p.standstill_max_acc_rms)
        {
            return false;
        }
        if (std::abs(aMean.norm() - p.gravity) > p.standstill_max_acc_norm_error)
        {
            return false;  // not just gravity
        }
        accAlongVar = alongVar / static_cast<double>(na);
    }

    // A standstill: update the estimate at the window end.
    auto& e = s.est;
    predict(e, t1);

    Mat33 Rg;
    for (int i = 0; i < 3; i++)
    {
        Rg(i, i) =
            std::max(wVar[i] / static_cast<double>(n), mrpt::square(p.gyro_measurement_min_sigma));
    }
    kalman_update(e.bias_gyro, e.cov_gyro, fromEigen(wMean), Rg);

    if (accAlongVar)
    {
        // |mean(a)| = |R^T g_up + b_acc| ~= g + u^T b_acc, so only u^T b_acc is observed:
        const double z = aMean.norm() - p.gravity;
        const double r = std::max(*accAlongVar, mrpt::square(p.acc_measurement_min_sigma));
        kalman_update_scalar(e.bias_acc, e.cov_acc, u, z, r);
    }

    e.last_update = t1;
    e.standstill_updates++;
    return true;
}

void ImuBiasTracker::add_motion_hint(TimeStamp t, bool isStatic) { hints_[t] = isStatic; }

void ImuBiasTracker::add_measurement(
    const std::string& sensorLabel, TimeStamp t,
    const std::optional<std::pair<AngularVelocity, mrpt::math::CMatrixDouble33>>&    gyro,
    const std::optional<std::pair<LinearAcceleration, mrpt::math::CMatrixDouble33>>& acc)
{
    if (!gyro && !acc)
    {
        return;
    }
    auto& e = imu_state(sensorLabel, t).est;
    if (t < e.stamp)
    {
        // An older measurement than the current state: fused as if it were now,
        // with the covariance it had (a conservative approximation).
        t = e.stamp;
    }
    predict(e, t);
    if (gyro)
    {
        kalman_update(e.bias_gyro, e.cov_gyro, gyro->first, gyro->second);
    }
    if (acc)
    {
        kalman_update(e.bias_acc, e.cov_acc, acc->first, acc->second);
    }
    e.last_update = std::max(e.last_update, t);
    e.external_updates++;
}

std::optional<ImuBiasTracker::Estimate> ImuBiasTracker::estimate(
    const std::string& sensorLabel, const std::optional<TimeStamp>& t) const
{
    const auto it = imus_.find(sensorLabel);
    if (it == imus_.end())
    {
        return std::nullopt;
    }
    Estimate e = it->second.est;
    if (t)
    {
        predict(e, *t);
    }
    return e;
}

std::vector<std::string> ImuBiasTracker::sensor_labels() const
{
    std::vector<std::string> out;
    out.reserve(imus_.size());
    for (const auto& [label, s] : imus_)
    {
        out.push_back(label);
    }
    return out;
}

void ImuBiasTracker::reset()
{
    imus_.clear();
    hints_.clear();
}
