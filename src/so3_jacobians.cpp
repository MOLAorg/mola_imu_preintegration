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
 * @file   so3_jacobians.cpp
 * @brief  SO(3) right Jacobian and small helpers used by IMU preintegration.
 * @author Jose Luis Blanco Claraco
 * @date   Jul 22, 2026
 */

#include <mola_imu_preintegration/so3_jacobians.h>
#include <mrpt/poses/Lie/SO.h>

#include <Eigen/Dense>
#include <cmath>

using namespace mola::imu;

namespace
{
/// Below this rotation angle [rad] the closed forms suffer cancellation and a
/// Taylor expansion is used instead. At this threshold the truncated terms are
/// of order 1e-17, i.e. already below double precision.
constexpr double SMALL_ANGLE_THRESHOLD = 1e-4;
}  // namespace

mrpt::math::CMatrixDouble33 mola::imu::skew_symmetric(const mrpt::math::TVector3D& v)
{
    mrpt::math::CMatrixDouble33 K;
    K(0, 0) = 0;
    K(0, 1) = -v.z;
    K(0, 2) = v.y;
    K(1, 0) = v.z;
    K(1, 1) = 0;
    K(1, 2) = -v.x;
    K(2, 0) = -v.y;
    K(2, 1) = v.x;
    K(2, 2) = 0;
    return K;
}

SO3 mola::imu::so3_exp(const mrpt::math::TVector3D& phi)
{
    return mrpt::poses::Lie::SO<3>::exp(mrpt::math::CVectorFixedDouble<3>(phi));
}

mrpt::math::TVector3D mola::imu::so3_log(const SO3& R)
{
    const auto v = mrpt::poses::Lie::SO<3>::log(R);
    return {v[0], v[1], v[2]};
}

mrpt::math::CMatrixDouble33 mola::imu::right_jacobian_so3(const mrpt::math::TVector3D& phi)
{
    const double t2 = phi.sqrNorm();
    const double t  = std::sqrt(t2);

    // Coefficients of  Jr = I - c1 * K + c2 * K^2
    double c1 = 0;
    double c2 = 0;
    if (t < SMALL_ANGLE_THRESHOLD)
    {
        // (1-cos t)/t^2  and  (t - sin t)/t^3
        c1 = 0.5 - t2 / 24.0;
        c2 = 1.0 / 6.0 - t2 / 120.0;
    }
    else
    {
        c1 = (1.0 - std::cos(t)) / t2;
        c2 = (t - std::sin(t)) / (t2 * t);
    }

    const auto K  = skew_symmetric(phi);
    const auto Ke = K.asEigen();

    mrpt::math::CMatrixDouble33 Jr;
    Jr.asEigen() = Eigen::Matrix3d::Identity() - c1 * Ke + c2 * (Ke * Ke);
    return Jr;
}

mrpt::math::CMatrixDouble33 mola::imu::inverse_right_jacobian_so3(const mrpt::math::TVector3D& phi)
{
    const double t2 = phi.sqrNorm();
    const double t  = std::sqrt(t2);

    // Coefficient of K^2 in  Jr^-1 = I + K/2 + c * K^2
    double c = 0;
    if (t < SMALL_ANGLE_THRESHOLD)
    {
        c = 1.0 / 12.0 + t2 / 720.0;
    }
    else
    {
        c = 1.0 / t2 - (1.0 + std::cos(t)) / (2.0 * t * std::sin(t));
    }

    const auto K  = skew_symmetric(phi);
    const auto Ke = K.asEigen();

    mrpt::math::CMatrixDouble33 Jinv;
    Jinv.asEigen() = Eigen::Matrix3d::Identity() + 0.5 * Ke + c * (Ke * Ke);
    return Jinv;
}
