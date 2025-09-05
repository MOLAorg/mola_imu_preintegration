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
 * @file   RotationIntegrator.cpp
 * @brief  Integrator of IMU accelerations and angular velocity readings.
 * @author Jose Luis Blanco Claraco
 * @date   Sep 20, 2021
 */

#include <mola_imu_preintegration/RotationIntegrator.h>
#include <mrpt/poses/Lie/SO.h>

using namespace mola;

void RotationIntegrator::initialize(const mrpt::containers::yaml& cfg)
{
    reset_integration();

    // Load params:
    params_.load_from(cfg);
}

void RotationIntegrator::reset_integration()
{
    // reset:
    state_ = IntegrationState();
}

void RotationIntegrator::integrate_measurement(const mrpt::math::TVector3D& w, double dt)
{
    const auto incrR = mola::incremental_rotation(w, params_, dt);

    // Update integration state:
    state_.deltaTij_ += dt;
    state_.deltaRij_ = state_.deltaRij_ * incrR;

    // TODO: Update Jacobian
}

mrpt::math::CMatrixDouble33 mola::incremental_rotation(
    const mrpt::math::TVector3D& w, const RotationIntegrationParams& params, double dt,
    const mrpt::optional_ref<mrpt::math::CMatrixDouble33>& D_incrR_integratedOmega)
{
    using mrpt::math::TVector3D;

    // Bias:
    TVector3D correctedW = w - params.gyroBias;

    // Translate to vehicle frame:
    if (params.sensorPose.has_value()) correctedW = params.sensorPose->rotateVector(correctedW);

    // Integrate:
    const TVector3D w_dt = correctedW * dt;

    if (D_incrR_integratedOmega.has_value())
    {
        // TODO: Jacobian: mrpt::poses::Lie::SO<3>::jacob_dexpe_de()
        THROW_EXCEPTION("Jacobian not implemented yet");
    }

    return mrpt::poses::Lie::SO<3>::exp(mrpt::math::CVectorFixedDouble<3>(w_dt));
}
