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

bool all_finite(const mrpt::math::TVector3D& v)
{
    return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

/// Larger backward jumps restart an IMU (e.g. a replayed dataset).
constexpr double CLOCK_JUMP_BACK = 1.0;  // [s]
/// Hints kept at most, should readings stop.
constexpr std::size_t MAX_HINTS = 10000;

}  // namespace

void ImuBiasTracker::Parameters::load_from(const mrpt::containers::yaml& cfg)
{
    MCP_LOAD_OPT(cfg, standstill_window);
    MCP_LOAD_OPT(cfg, standstill_min_samples);
    MCP_LOAD_OPT(cfg, max_sample_gap);
    MCP_LOAD_OPT(cfg, standstill_max_gyro_rms);
    MCP_LOAD_OPT(cfg, standstill_max_gyro_mean_offset);
    MCP_LOAD_OPT(cfg, standstill_gyro_offset_sigmas);
    MCP_LOAD_OPT(cfg, standstill_max_acc_rms);
    MCP_LOAD_OPT(cfg, standstill_max_acc_norm_error);
    MCP_LOAD_OPT(cfg, gravity);
    MCP_LOAD_OPT(cfg, require_static_hint);
    MCP_LOAD_OPT(cfg, motion_hint_hold);
    MCP_LOAD_OPT(cfg, max_hint_latency);
    MCP_LOAD_OPT(cfg, bias_gyro_random_walk);
    MCP_LOAD_OPT(cfg, bias_acc_random_walk);
    MCP_LOAD_OPT(cfg, bias_gyro_prior_sigma);
    MCP_LOAD_OPT(cfg, bias_acc_prior_sigma);
    MCP_LOAD_OPT(cfg, standstill_gyro_min_sigma);
    MCP_LOAD_OPT(cfg, standstill_acc_min_sigma);
    MCP_LOAD_OPT(cfg, bias_gyro_min_sigma);
    MCP_LOAD_OPT(cfg, bias_acc_min_sigma);
}

void ImuBiasTracker::Parameters::save_to(mrpt::containers::yaml& cfg) const
{
    MCP_SAVE(cfg, standstill_window);
    MCP_SAVE(cfg, standstill_min_samples);
    MCP_SAVE(cfg, max_sample_gap);
    MCP_SAVE(cfg, standstill_max_gyro_rms);
    MCP_SAVE(cfg, standstill_max_gyro_mean_offset);
    MCP_SAVE(cfg, standstill_gyro_offset_sigmas);
    MCP_SAVE(cfg, standstill_max_acc_rms);
    MCP_SAVE(cfg, standstill_max_acc_norm_error);
    MCP_SAVE(cfg, gravity);
    MCP_SAVE(cfg, require_static_hint);
    MCP_SAVE(cfg, motion_hint_hold);
    MCP_SAVE(cfg, max_hint_latency);
    MCP_SAVE(cfg, bias_gyro_random_walk);
    MCP_SAVE(cfg, bias_acc_random_walk);
    MCP_SAVE(cfg, bias_gyro_prior_sigma);
    MCP_SAVE(cfg, bias_acc_prior_sigma);
    MCP_SAVE(cfg, standstill_gyro_min_sigma);
    MCP_SAVE(cfg, standstill_acc_min_sigma);
    MCP_SAVE(cfg, bias_gyro_min_sigma);
    MCP_SAVE(cfg, bias_acc_min_sigma);
}

void ImuBiasTracker::Window::push(const Sample& s)
{
    samples.push_back(s);
    sum_w += s.w;
    sum_ww += mrpt::math::TVector3D(s.w.x * s.w.x, s.w.y * s.w.y, s.w.z * s.w.z);
    if (s.a)
    {
        na++;
        sum_a += *s.a;
        const Eigen::Vector3d a = toEigen(*s.a);
        sum_aa.asEigen() += a * a.transpose();
    }
}

void ImuBiasTracker::Window::pop()
{
    const Sample& s = samples.front();
    sum_w -= s.w;
    sum_ww -= mrpt::math::TVector3D(s.w.x * s.w.x, s.w.y * s.w.y, s.w.z * s.w.z);
    if (s.a)
    {
        na--;
        sum_a -= *s.a;
        const Eigen::Vector3d a = toEigen(*s.a);
        sum_aa.asEigen() -= a * a.transpose();
    }
    samples.pop_front();
}

ImuBiasTracker::PerImu& ImuBiasTracker::imu_state(const std::string& sensorLabel, TimeStamp t)
{
    auto it = imus_.find(sensorLabel);
    if (it == imus_.end())
    {
        PerImu s;
        s.est.cov_gyro = diag(parameters.bias_gyro_prior_sigma);
        s.est.cov_acc  = diag(parameters.bias_acc_prior_sigma);
        s.est.stamp    = t;
        s.newest       = t;
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
        e.cov_gyro(i, i) += mrpt::square(parameters.bias_gyro_random_walk) * dt;
        e.cov_acc(i, i) += mrpt::square(parameters.bias_acc_random_walk) * dt;
    }
    e.stamp = t;
}

void ImuBiasTracker::apply_floors(Estimate& e) const
{
    for (int i = 0; i < 3; i++)
    {
        e.cov_gyro(i, i) = std::max(e.cov_gyro(i, i), mrpt::square(parameters.bias_gyro_min_sigma));
        e.cov_acc(i, i)  = std::max(e.cov_acc(i, i), mrpt::square(parameters.bias_acc_min_sigma));
    }
}

bool ImuBiasTracker::add(const mrpt::obs::CObservationIMU& imu)
{
    using namespace mrpt::obs;
    if (!imu.has(IMU_WX) || !imu.has(IMU_WY) || !imu.has(IMU_WZ))
    {
        return false;
    }
    Sample smp;
    smp.t = mrpt::Clock::toDouble(imu.timestamp);
    smp.w = {imu.get(IMU_WX), imu.get(IMU_WY), imu.get(IMU_WZ)};
    if (!all_finite(smp.w) || !std::isfinite(smp.t))
    {
        return false;
    }
    if (imu.has(IMU_X_ACC) && imu.has(IMU_Y_ACC) && imu.has(IMU_Z_ACC))
    {
        const mrpt::math::TVector3D a(imu.get(IMU_X_ACC), imu.get(IMU_Y_ACC), imu.get(IMU_Z_ACC));
        if (all_finite(a))
        {
            smp.a = a;
        }
    }

    const bool isNew = imus_.count(imu.sensorLabel) == 0;
    auto&      s     = imu_state(imu.sensorLabel, smp.t);
    if (!isNew && smp.t <= s.newest)
    {
        if (s.newest - smp.t <= CLOCK_JUMP_BACK)
        {
            return false;  // duplicated or out of order: dropped
        }
        // The clock went back: restart this IMU, keeping its estimate.
        s.input.clear();
        s.window.clear();
        s.est.stamp = smp.t;
        if (s.est.has_measurements())
        {
            s.est.last_update = smp.t;
        }
        hints_.erase(hints_.upper_bound(smp.t), hints_.end());
        newest_hint_ = hints_.empty() ? 0 : hints_.rbegin()->first;
    }
    s.newest = smp.t;
    s.input.push_back(smp);
    const bool updated = drain(s);
    prune_hints();
    return updated;
}

bool ImuBiasTracker::drain(PerImu& s)
{
    const bool hintsInUse = hints_used_ || parameters.require_static_hint;
    bool       updated    = false;
    while (!s.input.empty())
    {
        const TimeStamp t = s.input.front().t;
        if (hintsInUse && t > newest_hint_ && s.newest - t < parameters.max_hint_latency)
        {
            break;  // wait for hints covering this reading
        }
        const Sample smp = s.input.front();
        s.input.pop_front();
        if (process_sample(s, smp))
        {
            updated = true;
        }
    }
    return updated;
}

bool ImuBiasTracker::process_sample(PerImu& s, const Sample& smp)
{
    auto& w = s.window;
    if (!w.samples.empty() && smp.t - w.samples.back().t > parameters.max_sample_gap)
    {
        w.clear();
    }
    w.push(smp);
    if (w.samples.back().t - w.samples.front().t < parameters.standstill_window)
    {
        return false;
    }
    if (try_update(s))
    {
        w.clear();  // accepted windows do not overlap
        return true;
    }
    // Slide: keep the newest readings, spanning less than a window.
    while (w.samples.size() > 1 &&
           w.samples.back().t - w.samples.front().t >= parameters.standstill_window)
    {
        w.pop();
    }
    return false;
}

bool ImuBiasTracker::hints_allow(TimeStamp from, TimeStamp to) const
{
    bool anyStatic = false;
    auto it        = hints_.upper_bound(from);
    // The hint before the window holds into it, for a while:
    if (it != hints_.begin())
    {
        const auto prev = std::prev(it);
        if (from - prev->first <= parameters.motion_hint_hold)
        {
            it = prev;
        }
    }
    for (; it != hints_.end() && it->first <= to; ++it)
    {
        if (!it->second)
        {
            return false;  // the vehicle was moving at some point
        }
        anyStatic = true;
    }
    return anyStatic || !parameters.require_static_hint;
}

bool ImuBiasTracker::try_update(PerImu& s)
{
    const auto& p   = parameters;
    const auto& win = s.window;
    const auto  n   = win.samples.size();
    if (n < p.standstill_min_samples || n < 2)
    {
        return false;
    }
    const TimeStamp t0 = win.samples.front().t;
    const TimeStamp t1 = win.samples.back().t;
    if (!hints_allow(t0, t1))
    {
        return false;
    }

    // Gyroscope: mean and per-axis variance.
    const double          nd    = static_cast<double>(n);
    const Eigen::Vector3d wMean = toEigen(win.sum_w) / nd;
    const Eigen::Vector3d wVar =
        ((toEigen(win.sum_ww) - nd * wMean.cwiseAbs2()) / (nd - 1)).cwiseMax(0.0);
    if (std::sqrt(wVar.sum()) > p.standstill_max_gyro_rms)
    {
        return false;  // vibrating or rotating: moving
    }
    auto&        e         = s.est;
    const double sigmaBias = std::sqrt(e.cov_gyro.asEigen().trace());
    const double maxOffset =
        std::max(p.standstill_max_gyro_mean_offset, p.standstill_gyro_offset_sigmas * sigmaBias);
    const double gyroOffset = (wMean - toEigen(e.bias_gyro)).norm();
    if (gyroOffset > maxOffset)
    {
        return false;  // a steady rotation, not a bias
    }

    // Accelerometer (optional): mean, spread and norm.
    std::optional<double> accAlongVar;
    Eigen::Vector3d       aMean = Eigen::Vector3d::Zero();
    Eigen::Vector3d       u     = Eigen::Vector3d::UnitZ();
    if (win.na >= 2)
    {
        const double na = static_cast<double>(win.na);
        aMean           = toEigen(win.sum_a) / na;
        const Eigen::Matrix3d aCov =
            (win.sum_aa.asEigen() - na * aMean * aMean.transpose()) / (na - 1);
        u = aMean.normalized();
        if (std::sqrt(std::max(aCov.trace(), 0.0)) > p.standstill_max_acc_rms)
        {
            return false;
        }
        if (std::abs(aMean.norm() - p.gravity) > p.standstill_max_acc_norm_error)
        {
            return false;  // not just gravity
        }
        accAlongVar = std::max(u.dot(aCov * u), 0.0) / na;
    }

    // A standstill: update the estimate at the window end.
    predict(e, t1);

    Mat33 Rg;
    for (int i = 0; i < 3; i++)
    {
        Rg(i, i) = std::max(wVar[i] / nd, mrpt::square(p.standstill_gyro_min_sigma));
    }
    kalman_update(e.bias_gyro, e.cov_gyro, fromEigen(wMean), Rg);

    if (accAlongVar)
    {
        // |mean(a)| = |R^T g_up + b_acc| ~= g + u^T b_acc, so only u^T b_acc is observed:
        const double z = aMean.norm() - p.gravity;
        const double r = std::max(*accAlongVar, mrpt::square(p.standstill_acc_min_sigma));
        kalman_update_scalar(e.bias_acc, e.cov_acc, u, z, r);
    }
    apply_floors(e);

    e.last_update = std::max(e.last_update, t1);
    e.standstill_updates++;
    return true;
}

std::vector<std::string> ImuBiasTracker::add_motion_hint(TimeStamp t, bool isStatic)
{
    std::vector<std::string> updated;
    if (!std::isfinite(t))
    {
        return updated;
    }
    // Two hints at the same time: "moving" wins.
    auto [it, inserted] = hints_.emplace(t, isStatic);
    if (!inserted)
    {
        it->second = it->second && isStatic;
    }
    hints_used_  = true;
    newest_hint_ = std::max(newest_hint_, t);
    for (auto& [label, s] : imus_)
    {
        if (drain(s))
        {
            updated.push_back(label);
        }
    }
    prune_hints();
    return updated;
}

void ImuBiasTracker::prune_hints()
{
    // Readings still to be checked start at the oldest window or waiting reading:
    std::optional<TimeStamp> oldest;
    for (const auto& [label, s] : imus_)
    {
        TimeStamp t = s.newest;
        if (!s.input.empty())
        {
            t = std::min(t, s.input.front().t);
        }
        if (!s.window.samples.empty())
        {
            t = std::min(t, s.window.samples.front().t);
        }
        oldest = oldest ? std::min(*oldest, t) : t;
    }
    if (oldest)
    {
        // Keep the newest hint before it, which holds into it:
        const TimeStamp limit = *oldest - parameters.motion_hint_hold;
        while (hints_.size() > 1 && std::next(hints_.begin())->first < limit)
        {
            hints_.erase(hints_.begin());
        }
    }
    while (hints_.size() > MAX_HINTS)
    {
        hints_.erase(hints_.begin());
    }
}

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
    // An older measurement than the estimate: the bias may have moved since.
    const double age = std::max(0.0, e.stamp - t);
    predict(e, t);
    if (gyro)
    {
        Mat33 R = gyro->second;
        for (int i = 0; i < 3; i++)
        {
            R(i, i) += mrpt::square(parameters.bias_gyro_random_walk) * age;
        }
        kalman_update(e.bias_gyro, e.cov_gyro, gyro->first, R);
    }
    if (acc)
    {
        Mat33 R = acc->second;
        for (int i = 0; i < 3; i++)
        {
            R(i, i) += mrpt::square(parameters.bias_acc_random_walk) * age;
        }
        kalman_update(e.bias_acc, e.cov_acc, acc->first, R);
    }
    apply_floors(e);
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
    predict(e, t ? *t : it->second.newest);
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
    hints_used_  = false;
    newest_hint_ = 0;
}
