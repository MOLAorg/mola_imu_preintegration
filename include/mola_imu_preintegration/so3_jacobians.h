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
 * @file   so3_jacobians.h
 * @brief  SO(3) right Jacobian and small helpers used by IMU preintegration.
 * @author Jose Luis Blanco Claraco
 * @date   Jul 22, 2026
 */
#pragma once

#include <mola_imu_preintegration/types.h>
#include <mrpt/math/CMatrixFixed.h>
#include <mrpt/math/TPoint3D.h>

namespace mola::imu
{
/** Skew-symmetric ("hat") matrix of a 3-vector, such that
 *  `skew_symmetric(a) * b == a.cross(b)`.
 */
mrpt::math::CMatrixDouble33 skew_symmetric(const mrpt::math::TVector3D& v);

/** SO(3) exponential map, `phi -> Exp(phi)`. Thin wrapper over
 *  `mrpt::poses::Lie::SO<3>::exp()` taking a TVector3D.
 */
SO3 so3_exp(const mrpt::math::TVector3D& phi);

/** SO(3) logarithm map, `R -> Log(R)`. Thin wrapper over
 *  `mrpt::poses::Lie::SO<3>::log()` returning a TVector3D.
 */
mrpt::math::TVector3D so3_log(const SO3& R);

/** Right Jacobian of SO(3), `Jr(phi)`, defined by the first-order identity
 *
 *    Exp(phi + delta) ~= Exp(phi) * Exp( Jr(phi) * delta )
 *
 *  Closed form:
 *
 *    Jr = I - ((1-cos t)/t^2) K + ((t - sin t)/t^3) K^2,   K = skew(phi), t = |phi|
 *
 *  A Taylor expansion is used for small `t` to avoid the catastrophic
 *  cancellation of `t - sin(t)`.
 *
 *  \note This is NOT what `mrpt::poses::Lie::SO<3>::jacob_dexpe_de()` returns:
 *        that one is the 9x3 Jacobian of the (column-major flattened) rotation
 *        matrix w.r.t. the tangent vector, a different object.
 */
mrpt::math::CMatrixDouble33 right_jacobian_so3(const mrpt::math::TVector3D& phi);

/** Inverse of right_jacobian_so3(), in closed form:
 *
 *    Jr^-1 = I + K/2 + ( 1/t^2 - (1+cos t)/(2 t sin t) ) K^2
 *
 *  with a Taylor expansion for small `t`.
 */
mrpt::math::CMatrixDouble33 inverse_right_jacobian_so3(const mrpt::math::TVector3D& phi);

}  // namespace mola::imu
