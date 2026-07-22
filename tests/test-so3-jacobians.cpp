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
 * @file   test-so3-jacobians.cpp
 * @brief  Unit tests for the SO(3) right Jacobian helpers.
 * @author Jose Luis Blanco Claraco
 * @date   Jul 22, 2026
 *
 * The right Jacobian is validated against its DEFINING identity,
 *   Exp(phi + delta) ~= Exp(phi) * Exp( Jr(phi) * delta ),
 * evaluated numerically, rather than against a second closed form: a typo
 * duplicated in both places would otherwise go unnoticed.
 */

#include <mola_imu_preintegration/so3_jacobians.h>
#include <mrpt/core/exceptions.h>
#include <mrpt/math/geometry.h>

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <vector>

using namespace mola::imu;
using mrpt::math::TVector3D;

namespace
{
/// Numerical Jacobian of  delta -> Log( Exp(phi)^T Exp(phi + delta) )  at delta=0,
/// which by definition equals Jr(phi).
Eigen::Matrix3d numericalRightJacobian(const TVector3D& phi, double eps = 1e-6)
{
    const auto            R  = so3_exp(phi);
    const Eigen::Matrix3d Rt = R.asEigen().transpose();

    Eigen::Matrix3d J = Eigen::Matrix3d::Zero();
    for (int i = 0; i < 3; i++)
    {
        TVector3D dp{0, 0, 0};
        dp[i] = eps;

        SO3 Dp;
        SO3 Dm;
        Dp.asEigen() = Rt * so3_exp(phi + dp).asEigen();
        Dm.asEigen() = Rt * so3_exp(phi - dp).asEigen();

        const TVector3D lp = so3_log(Dp);
        const TVector3D lm = so3_log(Dm);

        for (int r = 0; r < 3; r++)
        {
            J(r, i) = (lp[r] - lm[r]) / (2 * eps);
        }
    }
    return J;
}

const std::vector<TVector3D> testVectors()
{
    return {
        {0, 0, 0},    {1e-9, 0, 0},     {1e-5, -2e-5, 3e-6}, {0.1, 0, 0},     {0, 0.3, 0},
        {0, 0, -0.7}, {0.3, -0.4, 0.5}, {1.5, 0.2, -1.1},    {2.9, 0.1, 0.1},
    };
}

void test_skew_symmetric_is_cross_product()
{
    const TVector3D a{0.3, -1.2, 0.7};
    const TVector3D b{-0.5, 0.9, 2.1};

    const Eigen::Vector3d viaSkew = skew_symmetric(a).asEigen() * Eigen::Vector3d(b.x, b.y, b.z);

    TVector3D viaCross;
    mrpt::math::crossProduct3D(a, b, viaCross);

    ASSERT_NEAR_(viaSkew.x(), viaCross.x, 1e-12);
    ASSERT_NEAR_(viaSkew.y(), viaCross.y, 1e-12);
    ASSERT_NEAR_(viaSkew.z(), viaCross.z, 1e-12);

    std::cout << "test_skew_symmetric_is_cross_product passed." << std::endl;
}

void test_right_jacobian_matches_definition()
{
    for (const auto& phi : testVectors())
    {
        const Eigen::Matrix3d Jnum = numericalRightJacobian(phi);
        const Eigen::Matrix3d Jana = right_jacobian_so3(phi).asEigen();

        const double err = (Jnum - Jana).cwiseAbs().maxCoeff();
        ASSERTMSG_(
            err < 1e-6,
            mrpt::format("Jr mismatch for phi=%s : max abs error=%g", phi.asString().c_str(), err));
    }
    std::cout << "test_right_jacobian_matches_definition passed." << std::endl;
}

void test_inverse_right_jacobian_is_inverse()
{
    for (const auto& phi : testVectors())
    {
        const Eigen::Matrix3d J    = right_jacobian_so3(phi).asEigen();
        const Eigen::Matrix3d Jinv = inverse_right_jacobian_so3(phi).asEigen();

        const double err = (J * Jinv - Eigen::Matrix3d::Identity()).cwiseAbs().maxCoeff();
        ASSERTMSG_(
            err < 1e-9,
            mrpt::format(
                "Jr*Jr^-1 != I for phi=%s : max abs error=%g", phi.asString().c_str(), err));
    }
    std::cout << "test_inverse_right_jacobian_is_inverse passed." << std::endl;
}

void test_small_angle_branch_is_continuous()
{
    // Straddle the Taylor/closed-form threshold (1e-4) and check the two
    // branches agree, i.e. no discontinuity was introduced.
    const TVector3D dir{0.6, -0.48, 0.64};

    for (const double t : {0.99e-4, 1.01e-4})
    {
        const TVector3D phi = dir * t;

        const Eigen::Matrix3d Jnum = numericalRightJacobian(phi, 1e-7);
        const Eigen::Matrix3d Jana = right_jacobian_so3(phi).asEigen();

        const double err = (Jnum - Jana).cwiseAbs().maxCoeff();
        ASSERTMSG_(err < 1e-7, mrpt::format("Jr branch discontinuity at t=%g : error=%g", t, err));
    }
    std::cout << "test_small_angle_branch_is_continuous passed." << std::endl;
}

void test_exp_log_round_trip()
{
    for (const auto& phi : testVectors())
    {
        const TVector3D back = so3_log(so3_exp(phi));
        ASSERT_NEAR_(back.x, phi.x, 1e-10);
        ASSERT_NEAR_(back.y, phi.y, 1e-10);
        ASSERT_NEAR_(back.z, phi.z, 1e-10);
    }
    std::cout << "test_exp_log_round_trip passed." << std::endl;
}

}  // namespace

int main([[maybe_unused]] int argc, [[maybe_unused]] char** argv)
{
    try
    {
        test_skew_symmetric_is_cross_product();
        test_right_jacobian_matches_definition();
        test_inverse_right_jacobian_is_inverse();
        test_small_angle_branch_is_continuous();
        test_exp_log_round_trip();
        std::cout << "All SO(3) Jacobian tests passed." << std::endl;
        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Test failed: " << e.what() << std::endl;
        return 1;
    }
}
