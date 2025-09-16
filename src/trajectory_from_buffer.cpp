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
 * @file   trajectory_from_buffer.cpp
 * @brief  Reconstruct a trajectory from a LocalVelocityBuffer
 * @author Jose Luis Blanco Claraco
 * @date   Sep 5, 2025
 */

#include <mola_imu_preintegration/trajectory_from_buffer.h>

using namespace mola::imu;

std::string TrajectoryPoint::asString() const
{
    using mrpt::format;

    auto vecToStr = [](const mrpt::math::TVector3D& x)
    { return mrpt::format("[%.3f %.3f %.3f]", x[0], x[1], x[2]); };

    auto optVecToStr = [&](const std::optional<mrpt::math::TVector3D>& x)
    { return x ? vecToStr(*x) : std::string{"<none>"}; };

    auto matToStr = [](const mrpt::math::CMatrixDouble33& M)
    {
        mrpt::poses::CPose3D pp;
        pp.setRotationMatrix(M);
        return mrpt::format(
            "(ypr)=(%.02f,%.02f,%.02f) [deg]", mrpt::RAD2DEG(pp.yaw()), mrpt::RAD2DEG(pp.pitch()),
            mrpt::RAD2DEG(pp.roll()));
    };

    auto optMatToStr = [&](const std::optional<mrpt::math::CMatrixDouble33>& M)
    { return M ? matToStr(*M) : std::string{"<none>"}; };

    std::ostringstream oss;
    oss << "TrajectoryPoint{"
        << "\n  pose = " << pose.asString()  //
        << "\n  R_ga = " << optMatToStr(R_ga)  //
        << "\n  v    = " << optVecToStr(v)  //
        << "\n  w_b  = " << optVecToStr(w_b)  //
        << "\n  a_b  = " << optVecToStr(a_b)  //
        << "\n  ac_b = " << optVecToStr(ac_b)  //
        << "\n  alpha= " << vecToStr(alpha_b)  //
        << "\n  j_b  = " << vecToStr(j_b)  //
        << "\n}";
    return oss.str();
}
