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
 * @file   test-imu-averager.cpp
 * @brief  Unit test for ImuAverager: high-rate IMU decimation must average the
 *         readings of each period, not keep one of them, so vibration does not
 *         alias into the output.
 * @author Jose Luis Blanco Claraco
 */

#include <mola_imu_preintegration/ImuAverager.h>
#include <mrpt/core/exceptions.h>
#include <mrpt/system/datetime.h>

#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

using mola::imu::ImuAverager;

namespace
{
constexpr double T0 = 1000.0;  // [s]

mrpt::obs::CObservationIMU makeImu(double t, double w, double acc)
{
    mrpt::obs::CObservationIMU obs;
    obs.timestamp = mrpt::Clock::fromDouble(t);
    obs.set(mrpt::obs::IMU_WX, 0.0);
    obs.set(mrpt::obs::IMU_WY, 0.0);
    obs.set(mrpt::obs::IMU_WZ, w);
    obs.set(mrpt::obs::IMU_X_ACC, acc);
    obs.set(mrpt::obs::IMU_Y_ACC, 0.0);
    obs.set(mrpt::obs::IMU_Z_ACC, 9.81);
    return obs;
}

// A period of zero passes every reading through untouched.
void test_passthrough()
{
    ImuAverager avg;
    for (int i = 0; i < 5; i++)
    {
        const auto in  = makeImu(T0 + 0.001 * i, 0.1 * i, 0.2 * i);
        const auto out = avg.add(in, 0.0);
        ASSERT_(out.has_value());
        ASSERT_(out->timestamp == in.timestamp);
        ASSERT_EQUAL_(out->get(mrpt::obs::IMU_WZ), in.get(mrpt::obs::IMU_WZ));
    }
    std::cout << "test_passthrough: OK\n";
}

// Vibration alternating at the input rate, of amplitude A around the true
// values: each output must be the true value, not a +-A raw sample.
void test_vibration_is_averaged_out()
{
    constexpr double RATE   = 1000.0;  // [Hz]
    constexpr double PERIOD = 0.1;  // [s]
    constexpr double W_TRUE = 0.5;  // [rad/s]
    constexpr double A_TRUE = 1.5;  // [m/s^2]
    constexpr double VIB    = 3.0;  // vibration amplitude

    ImuAverager avg;
    size_t      nOut = 0;
    for (int i = 0; i < 1000; i++)
    {
        const double sign = (i % 2 == 0) ? 1.0 : -1.0;
        const double t    = T0 + i / RATE;
        const auto   out  = avg.add(makeImu(t, W_TRUE + sign * VIB, A_TRUE + sign * VIB), PERIOD);
        if (!out)
        {
            continue;
        }
        nOut++;
        // A batch has ~100 readings, so an unpaired one leaves at most VIB/100:
        ASSERT_LT_(std::abs(out->get(mrpt::obs::IMU_WZ) - W_TRUE), 0.05);
        ASSERT_LT_(std::abs(out->get(mrpt::obs::IMU_X_ACC) - A_TRUE), 0.05);
        ASSERT_LT_(std::abs(out->get(mrpt::obs::IMU_Z_ACC) - 9.81), 1e-9);
    }
    // One output per period, give or take the batch boundary:
    ASSERT_GE_(nOut, 9U);
    ASSERT_LE_(nOut, 10U);
    std::cout << "test_vibration_is_averaged_out: OK (" << nOut << " outputs)\n";
}

// The output is stamped at the mean of its readings' stamps, and each batch
// starts right after the previous one.
void test_output_stamp_is_batch_mean()
{
    ImuAverager         avg;
    std::vector<double> batch;
    std::vector<double> outStamps;
    std::vector<double> expectedStamps;
    for (int i = 0; i <= 40; i++)
    {
        const double t = T0 + 0.01 * i;
        batch.push_back(t);
        const auto out = avg.add(makeImu(t, 0, 0), 0.1);
        if (!out)
        {
            continue;
        }
        double sum = 0;
        for (const double s : batch)
        {
            sum += s;
        }
        expectedStamps.push_back(sum / static_cast<double>(batch.size()));
        outStamps.push_back(mrpt::Clock::toDouble(out->timestamp));
        batch.clear();
    }
    ASSERT_EQUAL_(outStamps.size(), expectedStamps.size());
    ASSERT_GE_(outStamps.size(), 3U);
    for (size_t k = 0; k < outStamps.size(); k++)
    {
        ASSERT_LT_(std::abs(outStamps[k] - expectedStamps[k]), 1e-5);
    }
    std::cout << "test_output_stamp_is_batch_mean: OK\n";
}

// Readings without a channel, or with non-finite values in it, do not enter
// that channel's average; the quaternion comes from the newest reading.
void test_missing_and_nonfinite_channels()
{
    ImuAverager avg;

    auto a      = makeImu(T0, 1.0, 1.0);
    auto b      = makeImu(T0 + 0.05, std::numeric_limits<double>::quiet_NaN(), 3.0);
    auto c      = mrpt::obs::CObservationIMU();
    c.timestamp = mrpt::Clock::fromDouble(T0 + 0.11);
    c.set(mrpt::obs::IMU_WX, 0.0);
    c.set(mrpt::obs::IMU_WY, 0.0);
    c.set(mrpt::obs::IMU_WZ, 3.0);  // no accelerometer in this one
    c.set(mrpt::obs::IMU_ORI_QUAT_W, 1.0);
    c.set(mrpt::obs::IMU_ORI_QUAT_X, 0.0);
    c.set(mrpt::obs::IMU_ORI_QUAT_Y, 0.0);
    c.set(mrpt::obs::IMU_ORI_QUAT_Z, 0.0);

    ASSERT_(!avg.add(a, 0.1).has_value());
    ASSERT_(!avg.add(b, 0.1).has_value());
    const auto out = avg.add(c, 0.1);
    ASSERT_(out.has_value());

    ASSERT_LT_(std::abs(out->get(mrpt::obs::IMU_WZ) - 2.0), 1e-9);  // mean of 1 and 3
    ASSERT_LT_(std::abs(out->get(mrpt::obs::IMU_X_ACC) - 2.0), 1e-9);  // mean of 1 and 3
    ASSERT_(out->has(mrpt::obs::IMU_ORI_QUAT_W));
    ASSERT_EQUAL_(out->get(mrpt::obs::IMU_ORI_QUAT_W), 1.0);
    std::cout << "test_missing_and_nonfinite_channels: OK\n";
}

void test_reset_discards_batch()
{
    ImuAverager avg;
    ASSERT_(!avg.add(makeImu(T0, 100.0, 0), 0.1).has_value());
    avg.reset();
    ASSERT_(!avg.add(makeImu(T0 + 1.0, 1.0, 0), 0.1).has_value());
    const auto out = avg.add(makeImu(T0 + 1.11, 1.0, 0), 0.1);
    ASSERT_(out.has_value());
    ASSERT_LT_(std::abs(out->get(mrpt::obs::IMU_WZ) - 1.0), 1e-9);
    std::cout << "test_reset_discards_batch: OK\n";
}
}  // namespace

int main()
{
    try
    {
        test_passthrough();
        test_vibration_is_averaged_out();
        test_output_stamp_is_batch_mean();
        test_missing_and_nonfinite_channels();
        test_reset_discards_batch();
        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
