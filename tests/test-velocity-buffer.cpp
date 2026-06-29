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
 * @file   test-mp2p_velocity_buffer.cpp
 * @brief  Unit tests for
 * @author Jose Luis Blanco Claraco
 * @date   Aug 16, 2025
 */

#include <mola_imu_preintegration/LocalVelocityBuffer.h>
#include <mrpt/containers/yaml.h>
#include <mrpt/core/exceptions.h>
#include <mrpt/core/round.h>
#include <mrpt/poses/Lie/SO.h>
#include <mrpt/system/os.h>

#include <iostream>

using namespace mola::imu;

namespace
{
void unit_test_basic_api()
{
    mola::imu::LocalVelocityBuffer buffer;

    const double reference_time = 1755345252;

    // Add some velocities
    buffer.add_linear_velocity(reference_time, {1.0, 2.0, 3.0});
    buffer.add_angular_velocity(reference_time, {0.1, 0.2, 0.3});

    // Check the contents
    ASSERT_EQUAL_(buffer.get_linear_velocities().size(), 1);
    ASSERT_EQUAL_(buffer.get_angular_velocities().size(), 1);

    // Check the reference time
    buffer.set_reference_zero_time(reference_time);
    ASSERT_EQUAL_(buffer.get_reference_zero_time(), reference_time);

    // Clear the buffer and check it's empty
    buffer.clear();
    ASSERT_EQUAL_(buffer.get_linear_velocities().size(), 0);
    ASSERT_EQUAL_(buffer.get_angular_velocities().size(), 0);

    std::cout << "✅ LocalVelocityBuffer unit_test_basic_api passed!" << std::endl;
}

void unit_test_basic_yaml()
{
    mola::imu::LocalVelocityBuffer buffer;

    const double reference_time = 1755345252;

    // Add some velocities
    buffer.add_linear_velocity(reference_time, {1.0, 2.0, 3.0});
    buffer.add_angular_velocity(reference_time, {0.1, 0.2, 0.3});
    buffer.add_linear_velocity(reference_time + 0.1, {4.0, 5.0, 6.0});
    buffer.add_angular_velocity(reference_time + 0.1, {0.4, 0.5, 0.6});

    buffer.set_reference_zero_time(reference_time);

    // Serialize to YAML
    const auto yaml    = buffer.toYAML();
    const auto yamlStr = [&]()
    {
        std::stringstream ss;
        ss << yaml;
        return ss.str();
    }();

    std::cout << "Serialized YAML:\n" << yamlStr << "\n\n";

    // Deserialize from YAML
    mola::imu::LocalVelocityBuffer buffer2;
    buffer2.fromYAML(mrpt::containers::yaml::FromText(yamlStr));

    std::cout << "Done parsing.\n";

    // Check that the contents are the same
    ASSERT_EQUAL_(buffer2.get_reference_zero_time(), buffer.get_reference_zero_time());
    ASSERT_EQUAL_(buffer2.get_linear_velocities().size(), buffer.get_linear_velocities().size());
    ASSERT_EQUAL_(buffer2.get_angular_velocities().size(), buffer.get_angular_velocities().size());

    for (const auto& [time, vel] : buffer.get_linear_velocities())
    {
        const auto it = buffer2.get_linear_velocities().find(time);
        ASSERT_(it != buffer2.get_linear_velocities().end());
        ASSERT_(it->second == vel);  // Assuming operator== is defined for LinearVelocity
    }

    for (const auto& [time, w] : buffer.get_angular_velocities())
    {
        const auto it = buffer2.get_angular_velocities().find(time);
        ASSERT_(it != buffer2.get_angular_velocities().end());
        ASSERT_(it->second == w);  // Assuming operator== is defined for AngularVelocity
    }
    std::cout << "✅ LocalVelocityBuffer unit_test_basic_yaml passed!" << std::endl;
}

void unit_test_yaml_roundtrip()
{
    LocalVelocityBuffer buf;

    // Fill with sample data
    buf.parameters.max_time_window = 2.0;
    buf.set_reference_zero_time(123.456);

    const LinearVelocity     v{1.0, 2.0, 3.0};
    const AngularVelocity    w{0.1, 0.2, 0.3};
    const LinearAcceleration a{9.8, 0.0, -9.8};
    const SO3                R = mrpt::poses::CPose3D(0, 0, 0, 0.7, 0.0, 0.0).getRotationMatrix();

    buf.add_linear_velocity(123.400, v);
    buf.add_angular_velocity(123.410, w);
    buf.add_linear_acceleration(123.420, a);
    buf.add_orientation(123.430, R);

    // Serialize
    const auto yml = buf.toYAML();

    yml.printAsYAML();

    // Deserialize into a new buffer
    LocalVelocityBuffer buf2;
    buf2.fromYAML(yml);

    // Checks
    ASSERT_EQUAL_(buf2.parameters.max_time_window, buf.parameters.max_time_window);
    ASSERT_EQUAL_(buf2.get_reference_zero_time(), buf.get_reference_zero_time());

    // linear velocities
    ASSERT_EQUAL_(buf2.get_linear_velocities().size(), buf.get_linear_velocities().size());
    for (const auto& [t, v1] : buf.get_linear_velocities())
    {
        const auto it = buf2.get_linear_velocities().find(t);
        ASSERT_(it != buf2.get_linear_velocities().end());
        ASSERT_(v1 == it->second);
    }

    // angular velocities
    ASSERT_EQUAL_(buf2.get_angular_velocities().size(), buf.get_angular_velocities().size());
    for (const auto& [t, w1] : buf.get_angular_velocities())
    {
        const auto it = buf2.get_angular_velocities().find(t);
        ASSERT_(it != buf2.get_angular_velocities().end());
        ASSERT_(w1 == it->second);
    }

    // linear accelerations
    ASSERT_EQUAL_(buf2.get_linear_accelerations().size(), buf.get_linear_accelerations().size());
    for (const auto& [t, a1] : buf.get_linear_accelerations())
    {
        const auto it = buf2.get_linear_accelerations().find(t);
        ASSERT_(it != buf2.get_linear_accelerations().end());
        ASSERT_(a1 == it->second);
    }

    // orientations
    ASSERT_EQUAL_(buf2.get_orientations().size(), buf.get_orientations().size());
    for (const auto& [t, R1] : buf.get_orientations())
    {
        const auto it = buf2.get_orientations().find(t);
        ASSERT_(it != buf2.get_orientations().end());
        mrpt::math::CMatrixDouble33 error_R = R1.inverse() * it->second;
        if (mrpt::poses::Lie::SO<3>::log(error_R).norm() > 1e-4)
        {
            std::cerr << "R1:\n" << R1 << "\nit->second:\n" << it->second << "\n";
            THROW_EXCEPTION("Error in rotation matrix");
        }
    }

    std::cout << "✅ LocalVelocityBuffer unit_test_yaml_roundtrip passed!" << std::endl;
}

void unit_test_yaml_epoch_precision()
{
    // Absolute UNIX-epoch timestamps (~1.7e9) with sub-millisecond fractions
    // exceed the precision of the legacy "%.09lf" key format. Verify the
    // lossless encoding recovers the exact double-valued keys on round-trip.
    LocalVelocityBuffer buf;
    buf.set_reference_zero_time(0.0);
    buf.parameters.max_time_window = 1.0e12;  // never prune in this test

    const std::vector<TimeStamp> stamps = {
        1755345252.123456789, 1755345252.987654321, 1755345253.000000001};

    for (size_t i = 0; i < stamps.size(); ++i)
    {
        buf.add_linear_velocity(stamps[i], {static_cast<double>(i), 0.0, 0.0});
    }

    LocalVelocityBuffer buf2;
    {
        std::stringstream ss;
        ss << buf.toYAML();
        buf2.fromYAML(mrpt::containers::yaml::FromText(ss.str()));
    }

    ASSERT_EQUAL_(buf2.get_linear_velocities().size(), stamps.size());
    for (const auto t : stamps)
    {
        const auto it = buf2.get_linear_velocities().find(t);
        if (it == buf2.get_linear_velocities().end())
        {
            THROW_EXCEPTION_FMT("Timestamp %.17g not recovered exactly after YAML round-trip", t);
        }
    }

    std::cout << "✅ LocalVelocityBuffer unit_test_yaml_epoch_precision passed!" << std::endl;
}

void unit_test_window_since()
{
    mola::imu::LocalVelocityBuffer buffer;

    const double t0 = 1000.0;

    // Five samples at t0, t0+0.1, t0+0.2, t0+0.3, t0+0.4 in each channel.
    for (int i = 0; i < 5; ++i)
    {
        const double t = t0 + 0.1 * i;
        buffer.add_linear_velocity(t, {static_cast<double>(i), 0.0, 0.0});
        buffer.add_linear_acceleration(t, {0.0, static_cast<double>(i), 0.0});
        buffer.add_angular_velocity(t, {0.0, 0.0, static_cast<double>(i)});
    }

    // 1) Open-ended window: strictly greater than `from`.
    {
        const auto w = buffer.window_since(t0 + 0.15);
        // Expects samples at t0+0.2, t0+0.3, t0+0.4 (3 of each).
        ASSERT_EQUAL_(w.v_b.size(), 3);
        ASSERT_EQUAL_(w.a_b.size(), 3);
        ASSERT_EQUAL_(w.w_b.size(), 3);

        const auto first_v = w.v_b.begin();
        ASSERT_NEAR_(first_v->first, t0 + 0.2, 1e-9);
        ASSERT_NEAR_(first_v->second.x, 2.0, 1e-9);
    }

    // 2) Bounded window (from, to]: strictly greater than `from`, less than
    //    or equal to `to`.
    {
        const auto w = buffer.window_since(t0, t0 + 0.2);
        // (1000, 1000.2] expects t0+0.1 and t0+0.2 (2 of each).
        ASSERT_EQUAL_(w.v_b.size(), 2);
        ASSERT_EQUAL_(w.a_b.size(), 2);
        ASSERT_EQUAL_(w.w_b.size(), 2);

        const auto last_a = w.a_b.rbegin();
        ASSERT_NEAR_(last_a->first, t0 + 0.2, 1e-9);
        ASSERT_NEAR_(last_a->second.y, 2.0, 1e-9);
    }

    // 3) Boundary semantics: `from` itself is excluded; `to` itself is
    //    included.
    {
        const auto w = buffer.window_since(t0, t0);
        ASSERT_EQUAL_(w.v_b.size(), 0);

        const auto w2 = buffer.window_since(t0 - 1.0, t0);
        ASSERT_EQUAL_(w2.v_b.size(), 1);
        ASSERT_NEAR_(w2.v_b.begin()->first, t0, 1e-9);
    }

    // 4) Empty result when `from` is past the last sample.
    {
        const auto w = buffer.window_since(t0 + 100.0);
        ASSERT_EQUAL_(w.v_b.size(), 0);
        ASSERT_EQUAL_(w.a_b.size(), 0);
        ASSERT_EQUAL_(w.w_b.size(), 0);
    }

    std::cout << "OK LocalVelocityBuffer unit_test_window_since passed!" << std::endl;
}

void unit_test_pruning_uses_latest_timestamp()
{
    // An out-of-order (older) insertion must not evict newer entries.
    // Window: 0.5 s.  Insert t=1000.0, t=1000.4, then t=999.9 (older).
    // After the late insert, the true latest is still 1000.4, so the
    // cutoff is 999.9 -- t=1000.0 must survive.
    LocalVelocityBuffer buf;
    buf.parameters.max_time_window = 0.5;

    buf.add_linear_velocity(1000.0, {1.0, 0.0, 0.0});
    buf.add_linear_velocity(1000.4, {2.0, 0.0, 0.0});
    // Out-of-order sample: slightly older than the first entry.
    buf.add_linear_velocity(999.9, {0.0, 0.0, 0.0});

    // latest = 1000.4, cutoff = 999.9 --> all three survive (1000.4-999.9=0.5,
    // exactly at the boundary; entries with age == max_time_window are kept).
    ASSERT_EQUAL_(buf.get_linear_velocities().size(), 3);

    // Now add a sample far enough ahead that t=999.9 should be pruned.
    buf.add_linear_velocity(1000.6, {3.0, 0.0, 0.0});
    // latest = 1000.6, cutoff = 1000.1 --> 999.9 and 1000.0 are evicted,
    // 1000.4 and 1000.6 survive.
    ASSERT_EQUAL_(buf.get_linear_velocities().size(), 2);
    ASSERT_(buf.get_linear_velocities().count(1000.4) == 1);
    ASSERT_(buf.get_linear_velocities().count(1000.6) == 1);

    std::cout << "OK LocalVelocityBuffer unit_test_pruning_uses_latest_timestamp passed!"
              << std::endl;
}

}  // namespace

int main([[maybe_unused]] int argc, [[maybe_unused]] char** argv)
{
    try
    {
        unit_test_basic_api();
        unit_test_basic_yaml();
        unit_test_yaml_roundtrip();
        unit_test_yaml_epoch_precision();
        unit_test_window_since();
        unit_test_pruning_uses_latest_timestamp();
    }
    catch (std::exception& e)
    {
        std::cerr << mrpt::exception_to_str(e) << "\n";
        return 1;
    }
}
