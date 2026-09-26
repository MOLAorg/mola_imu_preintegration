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
 * @file   ImuAverager.h
 * @brief  Averages high-rate IMU readings down to a lower, fixed rate.
 * @author Jose Luis Blanco Claraco
 */
#pragma once

#include <mrpt/core/Clock.h>
#include <mrpt/math/TPoint3D.h>
#include <mrpt/obs/CObservationIMU.h>

#include <cstddef>
#include <optional>

namespace mola::imu
{
/** Turns a high-rate IMU stream into one reading per `period` seconds, each the
 *  average of all readings received in that period.
 *
 *  Keeping one raw reading per period instead would alias any vibration faster
 *  than the output rate (e.g. motors or propellers) into the kept readings,
 *  which are then fused as if they were the platform motion.
 *
 *  Averaged: angular velocity and linear acceleration, each over the readings
 *  that have it with finite values. The orientation quaternion, if any, is
 *  absolute, so it is taken from the newest reading. The output timestamp is
 *  the mean of the batch timestamps, which is where an average of a smoothly
 *  varying signal is best referred to.
 *
 *  If readings must also be moved to the vehicle frame with ImuTransformer, do
 *  that first, on the raw stream: its lever-arm terms are not linear in the
 *  angular velocity, and it differentiates consecutive raw readings.
 */
class ImuAverager
{
   public:
    ImuAverager() = default;

    /** Adds one reading. Returns the average of the current batch once it spans
     *  at least `period` seconds (from its first to its newest reading), and
     *  starts a new batch. With `period <= 0`, every reading is returned as is.
     */
    [[nodiscard]] std::optional<mrpt::obs::CObservationIMU> add(
        const mrpt::obs::CObservationIMU& imu, double period);

    /** Discards the current batch. */
    void reset();

   private:
    std::optional<mrpt::Clock::time_point> batchStart_;
    double                                 sumStampOffset_ = 0;  // [s] from batchStart_
    std::size_t                            count_          = 0;

    mrpt::math::TVector3D sumAcc_ = {0, 0, 0};
    std::size_t           nAcc_   = 0;
    mrpt::math::TVector3D sumW_   = {0, 0, 0};
    std::size_t           nW_     = 0;
};

}  // namespace mola::imu
