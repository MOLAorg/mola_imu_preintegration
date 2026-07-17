^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Changelog for package mola_imu_preintegration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1.17.1 (2026-07-17)
-------------------
* Merge pull request `#10 <https://github.com/MOLAorg/mola_imu_preintegration/issues/10>`_ from MOLAorg/fix/imu-transformer-near-duplicate-timestamp-dt
  fix: guard against near-zero dt in ImuTransformer angular-accel finite difference
* style: clang-format fix
* fix: guard against near-zero dt in ImuTransformer angular-accel finite difference
  Some IMU drivers (e.g. the Hesai built-in IMU) emit consecutive messages with
  near-duplicate timestamps a few microseconds apart. Dividing the finite-
  difference angular-acceleration term by such a near-zero dt amplified an
  otherwise negligible angular-velocity delta into a huge spurious lever-arm
  correction, corrupting the transformed acceleration for that sample.
  Reproduced end-to-end: ImuInitialCalibrator's pitch/roll came out ~33 deg off
  on a dataset where the true tilt (per the bag's own static transforms) is
  ~2 deg, with about a third of samples in the calibration window affected.
  Widen the existing fallback-rate guard to treat any dt at or below 1 ms the
  same as the dt<=0 case. Add a regression test.
* Contributors: Jose Luis Blanco-Claraco

1.17.0 (2026-07-04)
-------------------
* Merge pull request `#9 <https://github.com/MOLAorg/mola_imu_preintegration/issues/9>`_ from MOLAorg/fix/trajectory-empty-on-insufficient-data
* Merge pull request `#8 <https://github.com/MOLAorg/mola_imu_preintegration/issues/8>`_ from MOLAorg/fix/imu-initial-calibrator-rotated-mount-orientation
* fix: apply IMU extrinsics to orientation in ImuInitialCalibrator
* Merge pull request `#7 <https://github.com/MOLAorg/mola_imu_preintegration/issues/7>`_ from MOLAorg/improve/imu-review-fixes
* ci: update and fix CI jobs
* fix: prune relative to true latest timestamp, not just-inserted time
* fix: clear() preserves parameters; documents what is reset
* docs: add agents.md with key insights for the package
* test: add unit tests for ImuTransformer
* feat: lossless timestamp encoding in LocalVelocityBuffer YAML
* refactor: remove unused tolerance_search_stamp parameter
* docs: fix incorrect doxygen comments
* fix: guard against empty IMU samples in trajectory_from_buffer
* docs: add ROS 2 Lyrical badge row, update Rolling to Ubuntu 26.04 (resolute)
* ci: fix Jazzy Jalisco EOL date in CI workflow comment (May 2024 - May 2029)
* Contributors: Jose Luis Blanco-Claraco

1.16.1 (2026-05-11)
-------------------
* Merge pull request `#6 <https://github.com/MOLAorg/mola_imu_preintegration/issues/6>`_ from MOLAorg/simplify-ci
  CI: simplify clang-format helpers and use ros: docker images for stable builds
* CI: simplify clang-format helpers and use ros: docker images for stable builds
  - Replace complex Python-based clang_git_format with a simple scripts/formatter.sh supporting --check
  - Standardize formatter to scripts/formatter.sh (renamed from clang-formatter.sh)
  - Simplify check-clang-format.yml to just apt-install clang-format-14 and run the script
  - Use ros:humble / ros:jazzy pre-built images for stable CI builds (faster, no setup-ros needed)
  Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
* Merge pull request `#5 <https://github.com/MOLAorg/mola_imu_preintegration/issues/5>`_ from MOLAorg/bump-cmake-version
  bump min req cmake version to 3.22
* CI: sensible job names
* bump min req cmake version to 3.22
* Merge pull request `#4 <https://github.com/MOLAorg/mola_imu_preintegration/issues/4>`_ from MOLAorg/feat/window_since2
  feat(mola_imu_preintegration): MOLA_IMU_PREINTEGRATION_HAS_WINDOW_SINCE macro + tests
* feat(mola_imu_preintegration): MOLA_IMU_PREINTEGRATION_HAS_WINDOW_SINCE macro + tests
  Expose MOLA_IMU_PREINTEGRATION_HAS_WINDOW_SINCE at the end of the
  LocalVelocityBuffer header so downstream packages in separate repos
  (mola_lidar_odometry's online gravity rebake) can guard usage with
  __has_include + this macro and stay buildable against older
  mola_imu_preintegration checkouts.
  Add unit tests for window_since covering open-ended and bounded
  windows, the (from, to] boundary semantics, and the empty-result
  case.
  Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
* Merge pull request `#3 <https://github.com/MOLAorg/mola_imu_preintegration/issues/3>`_ from MOLAorg/feat/window_since
  feat: Add API entry window_since()
* feat: Add API entry window_since()
* Contributors: Jose Luis Blanco-Claraco

1.16.0 (2026-04-22)
-------------------
* Merge pull request `#2 <https://github.com/MOLAorg/mola_imu_preintegration/issues/2>`_ from MOLAorg/feat/implement-low-pass-filters
  Implement low-pass filtering for acc/omega (more stable transformed IMU vectors)
* Add formal CLA
* Contributors: Jose Luis Blanco-Claraco

1.15.0 (2026-02-23)
-------------------
* Merge pull request `#1 <https://github.com/MOLAorg/mola_imu_preintegration/issues/1>`_ from MOLAorg/feat/transform-without-w
  ImuTransformer: work on acceleration even without omega
* BUG FIX: Wrong copy-paste sign error in centripetal acceleration
* ImuTransformer: work on acceleration even without omega
* Contributors: Jose Luis Blanco-Claraco

1.14.1 (2025-12-26)
-------------------
* Add option 'use_imu_orientation' to disable using IMU orientation for initialization
* CI: Use standard checkout action
* Contributors: Jose Luis Blanco-Claraco

1.14.0 (2025-11-11)
-------------------
* BUGFIX: Estimated acc. bias and initial roll had a bug in their formulas.
* Add unit tests for initial calibrator
* Update project website URL in package.xml
* Contributors: Jose Luis Blanco-Claraco

1.13.2 (2025-10-20)
-------------------
* Make use of ConstPtr in API
* Contributors: Jose Luis Blanco-Claraco

1.13.1 (2025-10-15)
-------------------
* ImuCalibrator: Robust against invalid orientation data in Imu msgs
* Fix random failures in unit test
* CI: add running unit tests
* Contributors: Jose Luis Blanco-Claraco

1.13.0 (2025-10-11)
-------------------
* FIX: Picked wrong reference stamp for integrating trajectories
* Add asString() debugging methods to sample structs
* Contributors: Jose Luis Blanco-Claraco

1.12.0 (2025-10-09)
-------------------
* FIX: May leave trajectory poses without populating raw IMU data
* ImuInitialCalibrator: Use IMU orientation, if available
* LocalVelocityBuffer: reduce default max_time_window to 0.5 s
* ImuIntegrationParams: Add save_to() method, and finish missing load_from() fields
* Fix comment typos
* Readme: Add missing ROS build farm badges for bin packages
* Contributors: Jose Luis Blanco-Claraco

1.11.0 (2025-09-24)
-------------------
* New unit tests for IMU integration
* Move everything into namespace mola::imu to avoid ns pollution
* Move LocalVelocityBuffer class here from mp2p_icp repository
* Contributors: Jose Luis Blanco-Claraco

1.10.0 (2025-09-07)
-------------------
* Add class ImuTransformer
* Update copyright notice
* Import IMU initialization class ImuInitialCalibrator, refactored from the mola LO package
* Contributors: Jose Luis Blanco-Claraco

1.9.0 (2025-06-06)
------------------

1.8.1 (2025-05-25)
------------------
* fixes for clang-tidy
* Contributors: Jose Luis Blanco-Claraco

1.8.0 (2025-03-15)
------------------

1.7.0 (2025-02-22)
------------------

1.6.1 (2025-01-10)
------------------
* Fix package.xml URLs
* Contributors: Jose Luis Blanco-Claraco

1.6.0 (2025-01-03)
------------------

1.5.0 (2024-12-26)
------------------

1.4.1 (2024-12-20)
------------------

1.4.0 (2024-12-18)
------------------

1.3.0 (2024-12-11)
------------------

1.2.1 (2024-09-29)
------------------

1.2.0 (2024-09-16)
------------------

1.1.3 (2024-08-28)
------------------
* Depend on new mrpt_lib packages (deprecate mrpt2)
* Contributors: Jose Luis Blanco-Claraco

1.1.2 (2024-08-26)
------------------

1.1.1 (2024-08-23)
------------------

1.1.0 (2024-08-18)
------------------
* Update clang-format style; add reformat bash script
* Merge pull request `#62 <https://github.com/MOLAorg/mola/issues/62>`_ from MOLAorg/docs-fixes
  Docs fixes
* Fix ament_xmllint warnings in package.xml
* Contributors: Jose Luis Blanco-Claraco

1.0.8 (2024-07-29)
------------------
* ament_lint_cmake: clean warnings
* Contributors: Jose Luis Blanco-Claraco

1.0.7 (2024-07-24)
------------------

1.0.6 (2024-06-21)
------------------

1.0.5 (2024-05-28)
------------------

1.0.4 (2024-05-14)
------------------
* bump cmake_minimum_required to 3.5
* Contributors: Jose Luis Blanco-Claraco

1.0.3 (2024-04-22)
------------------
* Update docs
* Fix package.xml website URL
* Contributors: Jose Luis Blanco-Claraco

1.0.2 (2024-04-04)
------------------

1.0.1 (2024-03-28)
------------------

1.0.0 (2024-03-19)
------------------
* API changes for new package mola_navstate_fuse
* Contributors: Jose Luis Blanco-Claraco

0.2.2 (2023-09-08)
------------------
* Update copyright year
* Correct references to license
* Ported to ROS2 colcon build system
* Delete WIP files.
* first unit tests
* progress, unit tests
* Contributors: Jose Luis Blanco-Claraco

0.2.1 (2021-09-18)
------------------
* Initial commit
* Contributors: Jose Luis Blanco-Claraco
