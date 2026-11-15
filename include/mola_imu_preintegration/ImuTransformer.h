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
 * @file   ImuTransformer.h
 * @brief  Transforms IMU readings from an arbitrary sensor pose to the vehicle reference base_link
 * @author Jose Luis Blanco Claraco
 * @date   Sep 5, 2025
 */

#pragma once

#include <mrpt/math/TPoint3D.h>
#include <mrpt/obs/obs_frwds.h>

namespace mola::imu
{
/** Transforms IMU readings from an arbitrary sensor pose to the vehicle reference "base_link".
 *
 * This basically rotates the angular velocity vector according to `raw_imu.sensorPose` and
 * transforms the acceleration vector according to the lever between vehicle "base_link" and the
 * sensor pose.
 * The returned IMU reading object has `sensorPose` set to the SE(3) identity to reflect its new
 * frame.
 *
 * Note that bias is *not* subtracted here, it must be done by the caller.
 * Also, gravity is not subtracted here since this class has no knowledge of global orientation.
 *
 * Since this class has a state (the stamp of the last observation), one object must be instantiated
 * for each IMU sensor.
 *
 * \ingroup mola_imu_preintegration_grp
 */
class ImuTransformer
{
   public:
    ImuTransformer() = default;

    struct Parameters
    {
        /** Cut-off frequency [Hz] of the low-pass filter applied to the angular
         * acceleration, which is estimated by finite differences of the gyroscope
         * readings and used in the tangential lever-arm term α×t.
         *
         * The filter is a first-order exponential moving average whose
         * per-sample factor follows the actual time step of each reading,
         * `alpha = 1 - exp(-2·π·f_c·dt)`, so its response does not depend on the
         * IMU rate. A value <= 0 disables the filter (raw finite differences).
         *
         * The default (20 Hz) keeps most of the true angular acceleration of
         * hand-held or legged motion, which a cut-off around 1 Hz removes.
         */
        double ang_acc_lpf_cutoff_hz = 20.0;

        /** Cut-off frequency [Hz] of the low-pass filter applied to the copy of
         * the angular velocity used in the centripetal lever-arm term ω×(ω×t).
         * Same filter type and conventions as `ang_acc_lpf_cutoff_hz`.
         *
         * The angular-velocity channels written to the output observation
         * (IMU_WX/WY/WZ) are always the unfiltered, rotated readings.
         */
        double ang_vel_lpf_cutoff_hz = 20.0;

        /** Legacy per-sample smoothing factor of the angular-acceleration
         * filter, in (0, 1]. If > 0, it is used for every sample instead of
         * `ang_acc_lpf_cutoff_hz`, as in earlier versions (values above 1 are
         * clamped to 1); its cut-off then depends on the IMU rate:
         * `f_c = -fs/(2·π)·ln(1 - alpha)`.
         * Default: 0, meaning unset (use the cut-off frequency).
         */
        double ang_acc_lpf_alpha = 0;

        /** Legacy per-sample smoothing factor of the angular-velocity filter,
         * with the same conventions as `ang_acc_lpf_alpha`.
         * Default: 0 (unset: use `ang_vel_lpf_cutoff_hz`).
         */
        double ang_vel_lpf_alpha = 0;
    };

    Parameters parameters;

    mrpt::obs::CObservationIMU process(const mrpt::obs::CObservationIMU& raw_imu);

   private:
    double                last_stamp_   = 0;
    bool                  first_sample_ = true;  ///< guards against potential first-call spike
    mrpt::math::TVector3D last_ang_vel_body_{0, 0, 0};
    mrpt::math::TVector3D filtered_ang_acc_{0, 0, 0};  ///< EMA state for angular acceleration
    mrpt::math::TVector3D filtered_ang_vel_{0, 0, 0};  ///< EMA state for ω used in centripetal
};

}  // namespace mola::imu
