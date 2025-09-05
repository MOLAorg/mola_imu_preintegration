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
 * @file   ImuInitialCalibrator.h
 * @brief  Provides a rough initial calibration and attitude for IMUs without
 *         proper bias calibration
 * @author Jose Luis Blanco Claraco
 * @date   Sep 5, 2025
 */

#pragma once

#include <mrpt/core/exceptions.h>
#include <mrpt/obs/obs_frwds.h>

#include <cstdlib>
#include <map>
#include <memory>

namespace mola
{
/** Provides a rough initial calibration and attitude for IMUs without proper
 * bias calibration.
 *
 * \ingroup mola_imu_preintegration_grp
 */
class ImuInitialCalibrator
{
   public:
    ImuInitialCalibrator(
        const std::size_t required_samples,  // NOLINT
        const double      max_samples_age)
        : required_samples_(required_samples), max_samples_age_(max_samples_age)
    {
        ASSERT_(required_samples > 0);
    }

    void add(const std::shared_ptr<const mrpt::obs::CObservationIMU>& obs);

    [[nodiscard]] bool                       isReady() const;
    [[nodiscard]] std::tuple<double, double> getPitchRoll() const;

   private:
    std::size_t required_samples_;
    double      max_samples_age_;

    std::map<double, std::shared_ptr<const mrpt::obs::CObservationIMU>> samples_;
};

}  // namespace mola
