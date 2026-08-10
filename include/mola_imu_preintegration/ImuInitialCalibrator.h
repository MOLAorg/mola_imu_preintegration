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

#include <mola_imu_preintegration/ImuTransformer.h>
#include <mrpt/obs/CObservationIMU.h>

#include <cstdlib>
#include <map>
#include <optional>
#include <string>

namespace mola::imu
{
/** Provides a rough initial calibration and attitude for IMUs without a proper
 *  bias calibration.
 *
 * \ingroup mola_imu_preintegration_grp
 */
class ImuInitialCalibrator
{
   public:
    ImuInitialCalibrator() = default;

    struct Parameters
    {
        Parameters() = default;

        std::size_t required_samples = 50;  //!< Minimum required samples
        double      max_samples_age  = .75;  //!< Maximum samples age [seconds]
        double      gravity          = 9.81;  //<! Gravity magnitude [m/s²]

        /** If provided by the IMU, prefer gravity-aligned orientation from the sensor instead of
         * accelerometer data. Should be normally preferable, so default is true.
         */
        bool use_imu_orientation = true;

        /** Length of the averaging window, in **seconds** [s]. This is the knob that actually
         * determines the accuracy of the estimated attitude, since the error is dominated by
         * platform motion during the window and not by sensor noise. Being a duration, it is
         * independent of the IMU rate: `required_samples` alone means a different averaging
         * time on every sensor, and cannot even be satisfied by a low-rate IMU if it does not
         * fit within `max_samples_age`.
         *
         * When >0, `required_samples` degrades to a minimum-sample sanity floor, the buffer is
         * pruned to (slightly more than) this window instead of `max_samples_age`, and the
         * calibration averages only the samples inside it.
         *
         * 0 (default) disables it and restores the legacy sample-count-only readiness rule.
         */
        double window_seconds = .0;

        /** Maximum RMS angular dispersion [rad] of the buffered accelerometer directions for the
         * window to be accepted as actually measuring gravity. A window taken while the platform
         * is being jostled freezes a reading that is not gravity, so isReady() defers and the
         * caller is expected to retry later with a fresher window.
         *
         * 0 (default) disables the gate.
         */
        double max_direction_dispersion = .0;

        /** How long [s] the dispersion gate may defer readiness before it gives up and accepts
         * the most recent window anyway. Measured from the first sample ever inserted, in sensor
         * time. Without it, a platform that never becomes quiet would never initialize.
         *
         * 0 (default) means "no timeout", i.e. wait indefinitely for a quiet window.
         */
        double dispersion_timeout = .0;
    };

    Parameters parameters;

    /// Initial calibration results from getCalibration()
    struct Results
    {
        Results() = default;

        mrpt::math::TVector3D bias_gyro{0, 0, 0};  //!< Gyroscope bias
        mrpt::math::TVector3D bias_acc_b{0, 0, 0};  //!< Accelerometer bias, in base_link frame

        /// Gyroscope additive noise standard deviation
        mrpt::math::TVector3D noise_stddev_gyro{0, 0, 0};
        /// Accelerometer additive noise standard deviation, in the base_link/body frame
        mrpt::math::TVector3D noise_stddev_acc{0, 0, 0};

        double pitch = 0;  //!< Estimated pitch angle, assuming being at rest during calibration
        double roll  = 0;  //!< Estimated roll angle, assuming being at rest during calibration

        /// Prints to a string a human-readable representation of all values
        std::string asString() const;
    };

    /// Detailed outcome of the readiness check, for logging and testing.
    struct Readiness
    {
        Readiness() = default;

        bool ready = false;  //!< Whether getCalibration() would return a value

        std::size_t samples   = 0;  //!< Samples inside the averaging window
        double      span      = 0;  //!< Time span of the averaging window [s]
        double      elapsed   = 0;  //!< Time since the first sample ever inserted [s]
        bool        timed_out = false;  //!< Accepted only because the gate timed out
        std::string reason;  //!< Why it is not ready yet (empty if ready)

        /// RMS angular dispersion of the accelerometer directions [rad], if computable
        std::optional<double> dispersion;
    };

    /// Inserts an IMU observation into the queue
    void add(const mrpt::obs::CObservationIMU::ConstPtr& obs);

    /// Returns true if there are already samples enough in the buffer to call getCalibration()
    [[nodiscard]] bool isReady() const;

    /// Same as isReady(), plus the quantities the decision was made on
    [[nodiscard]] Readiness readiness() const;

    /// RMS angular spread [rad] of the buffered accelerometer directions about their mean,
    /// over the averaging window. A quasi-static platform gives ~0; anything larger means the
    /// accelerometer is measuring platform motion on top of gravity. Returns nothing if there
    /// are fewer than 2 usable samples.
    [[nodiscard]] std::optional<double> directionDispersion() const;

    /// If enough samples are given, it computes the initial rough IMU calibration
    [[nodiscard]] std::optional<Results> getCalibration() const;

   private:
    std::map<std::string /*sensorLabel*/, ImuTransformer> imu_transformers_;

    /// Samples here have been already transformed to be on the base_link frame
    /// (accel, gyro AND the absolute-orientation quaternion; placeholder identity
    /// orientations from non-attitude IMUs are dropped on insertion):
    std::map<double, const mrpt::obs::CObservationIMU> samples_;

    /// Timestamp of the first sample ever inserted, to measure the gate timeout in sensor time
    std::optional<double> first_sample_time_;

    /// How far back samples are kept [s]. In time-window mode this is the window itself (plus a
    /// small margin, so the buffer can actually span it), and `max_samples_age` is unused.
    [[nodiscard]] double bufferHorizon() const;

    using const_iterator = std::map<double, const mrpt::obs::CObservationIMU>::const_iterator;

    /// First sample of the averaging window: in time-window mode, the oldest one no older than
    /// `window_seconds`; otherwise the whole buffer.
    [[nodiscard]] const_iterator windowBegin() const;
};

// Remove when all distros have the version with "use_imu_orientation"
#define MOLA_IMU_PREINT_HAS_USE_IMU_ORIENT_PARAM

// Remove when all distros have the version with the time window and dispersion gate
#define MOLA_IMU_PREINT_HAS_INIT_TIME_WINDOW

}  // namespace mola::imu
