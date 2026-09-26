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
 * @file   ImuAverager.cpp
 * @brief  Averages high-rate IMU readings down to a lower, fixed rate.
 * @author Jose Luis Blanco Claraco
 */

#include <mola_imu_preintegration/ImuAverager.h>
#include <mrpt/system/datetime.h>

#include <cmath>

using namespace mola::imu;

namespace
{
bool readVector(
    const mrpt::obs::CObservationIMU& imu, mrpt::obs::TIMUDataIndex ix, mrpt::obs::TIMUDataIndex iy,
    mrpt::obs::TIMUDataIndex iz, mrpt::math::TVector3D& out)
{
    if (!imu.has(ix) || !imu.has(iy) || !imu.has(iz))
    {
        return false;
    }
    out = {imu.get(ix), imu.get(iy), imu.get(iz)};
    return std::isfinite(out.x) && std::isfinite(out.y) && std::isfinite(out.z);
}
}  // namespace

std::optional<mrpt::obs::CObservationIMU> ImuAverager::add(
    const mrpt::obs::CObservationIMU& imu, double period)
{
    if (period <= 0)
    {
        return imu;
    }

    if (!batchStart_)
    {
        batchStart_ = imu.timestamp;
    }

    const double offset = mrpt::system::timeDifference(*batchStart_, imu.timestamp);
    sumStampOffset_ += offset;
    count_++;

    mrpt::math::TVector3D v;
    if (readVector(imu, mrpt::obs::IMU_X_ACC, mrpt::obs::IMU_Y_ACC, mrpt::obs::IMU_Z_ACC, v))
    {
        sumAcc_ += v;
        nAcc_++;
    }
    if (readVector(imu, mrpt::obs::IMU_WX, mrpt::obs::IMU_WY, mrpt::obs::IMU_WZ, v))
    {
        sumW_ += v;
        nW_++;
    }

    if (offset < period)
    {
        return std::nullopt;
    }

    // The batch is complete: the newest reading carries everything that is not
    // averaged (sensor pose, label, orientation), then the averages go on top.
    mrpt::obs::CObservationIMU out = imu;
    out.timestamp =
        mrpt::system::timestampAdd(*batchStart_, sumStampOffset_ / static_cast<double>(count_));

    if (nAcc_ > 0)
    {
        const auto a = sumAcc_ * (1.0 / static_cast<double>(nAcc_));
        out.set(mrpt::obs::IMU_X_ACC, a.x);
        out.set(mrpt::obs::IMU_Y_ACC, a.y);
        out.set(mrpt::obs::IMU_Z_ACC, a.z);
    }
    if (nW_ > 0)
    {
        const auto w = sumW_ * (1.0 / static_cast<double>(nW_));
        out.set(mrpt::obs::IMU_WX, w.x);
        out.set(mrpt::obs::IMU_WY, w.y);
        out.set(mrpt::obs::IMU_WZ, w.z);
    }

    reset();
    return out;
}

void ImuAverager::reset() { *this = ImuAverager(); }
