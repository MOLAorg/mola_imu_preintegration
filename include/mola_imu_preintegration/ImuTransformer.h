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
        /** Smoothing factor for the angular-acceleration EMA low-pass filter.
         *
         * Range (0, 1]:
         *   - Values close to 1: almost no filtering (raw finite difference).
         *   - Values close to 0: heavy smoothing (very slow response).
         *
         * The exact -3 dB cut-off frequency of the EMA filter is:
         *
         * \verbatim
         *   f_c = -fs / (2·π) · ln(1 - alpha)
         * \endverbatim
         *
         * Approximate cut-off frequencies for common IMU rates:
         *
         * \verbatim
         *   alpha |  fs = 100 Hz  |  fs = 400 Hz  |  fs = 600 Hz
         *   ------+--------------+---------------+--------------
         *   0.02  |    0.32 Hz   |    1.29 Hz    |    1.93 Hz  (default)
         *   0.05  |    0.82 Hz   |    3.27 Hz    |    4.90 Hz
         *   0.10  |    1.68 Hz   |    6.71 Hz    |   10.06 Hz
         *   0.20  |    3.55 Hz   |   14.21 Hz    |   21.31 Hz
         *   0.30  |    5.68 Hz   |   22.71 Hz    |   34.06 Hz
         * \endverbatim
         *
         * The default (0.02) targets a cut-off of ~0.32 Hz at 100 Hz, ~1.29 Hz at
         * 400 Hz, and ~1.93 Hz at 600 Hz.  This rejects amplified gyro noise from
         * the finite-difference step while remaining responsive to genuine
         * angular-acceleration transients at typical robot motion bandwidths.
         * Raise alpha if the lever-arm correction feels sluggish; lower it if
         * accelerometer readings are still noisy after transformation.
         */
        double ang_acc_lpf_alpha = 0.02;

        /** Smoothing factor for the angular-velocity EMA low-pass filter used
         * in the centripetal lever-arm correction term ω×(ω×t).
         *
         * Range (0, 1]:
         *   - Values close to 1: almost no filtering (raw rotated gyro reading).
         *   - Values close to 0: heavy smoothing (very slow response).
         *
         * Important: this filter is applied only to the ω copy used internally
         * in the lever-arm cross-product.  The angular-velocity channels written
         * to the output observation (IMU_WX/WY/WZ) are always the unfiltered,
         * rotated readings so that downstream consumers (e.g. preintegration)
         * receive the full-bandwidth signal.
         *
         * The exact -3 dB cut-off frequency of the EMA filter is:
         *
         * \verbatim
         *   f_c = -fs / (2·π) · ln(1 - alpha)
         * \endverbatim
         *
         * Approximate cut-off frequencies for common IMU rates:
         *
         * \verbatim
         *   alpha |  fs = 100 Hz  |  fs = 400 Hz  |  fs = 600 Hz
         *   ------+--------------+---------------+--------------
         *   0.02  |    0.32 Hz   |    1.29 Hz    |    1.93 Hz
         *   0.05  |    0.82 Hz   |    3.27 Hz    |    4.90 Hz
         *   0.10  |    1.68 Hz   |    6.71 Hz    |   10.06 Hz
         *   0.20  |    3.55 Hz   |   14.21 Hz    |   21.31 Hz
         *   0.30  |    5.68 Hz   |   22.71 Hz    |   34.06 Hz   (default)
         * \endverbatim
         *
         * A higher default than ang_acc_lpf_alpha is appropriate here because
         * ω is a direct gyro measurement (not a finite difference), so it
         * carries far less amplified noise.  The default (0.30) gives ~5.7 Hz
         * at 100 Hz and ~22.7 Hz at 400 Hz, which tracks real angular-velocity
         * transients faithfully while still attenuating high-frequency noise
         * in the quadratic centripetal term.
         */
        double ang_vel_lpf_alpha = 0.3;
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
