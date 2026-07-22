# agents.md — `mola_imu_preintegration`

Minimal key insights for working in this repo. Keep this file short and in sync with the code.

## What this package is
A lightweight MOLA package (ROS 2, depends on MRPT 2.x) for IMU data manipulation, rough
calibration, and gyro-based trajectory reconstruction for LiDAR deskew. **It is not a full
visual-inertial preintegrator** (no ΔV/ΔP, no covariance, no bias Jacobians yet).

## Layout
- `include/mola_imu_preintegration/` — public headers (namespace `mola::imu`).
- `src/` — implementations. `register.cpp` is an (empty) RTTI initializer.
- `tests/` — `mola_add_test(...)` per file; run via colcon/ctest.

## Core components
- `ImuIntegrator` + `incremental_rotation()` — gyro-only rotation preintegration,
  `ΔR ← ΔR · Exp((ω − b_g)·dt)`. Bias subtracted in IMU frame, then rotated to body via `sensor_pose`.
- `ImuTransformer` — moves raw readings to `base_link`: rotates ω/accel and applies the rigid-body
  lever arm `a_body = R·a_imu − α×t − ω×(ω×t)`. Stateful (one instance per IMU); uses EMA low-pass
  filters for α (finite-diff) and the centripetal ω. Output `IMU_W*` channels stay full-bandwidth.
- `ImuInitialCalibrator` — at-rest bias + pitch/roll from averaged gravity; prefers IMU orientation
  quaternion when present. Accel bias is only observable along gravity unless orientation is used.
- `LocalVelocityBuffer` — time-keyed ring of v/a/ω/orientation; YAML (de)serialize; window queries.
- `trajectory_from_buffer()` — dead-reckons an SE(3) trajectory around `t=0` from one orientation
  anchor + one velocity anchor + IMU ω/accel. Specific force → coordinate accel via
  `ac_b = a_b + Rᵀ_ga·g_world`. `use_higher_order` adds jerk/α Taylor terms.

## Conventions / gotchas
- Frames: angular/linear quantities in body (`base_link`) frame; orientations are global,
  gravity-aligned. Gravity vector default `(0,0,-9.81)`; accel-at-rest reads `(0,0,+g)`.
- Quaternion YAML order is `[x,y,z,w]` (ROS/REP-103).
- `cov_gyro`/`cov_acc` params exist but are **not yet consumed** by any algorithm.
- `trajectory_from_buffer` uses the sample closest to the **latest** timestamp as the anchor and
  integrates forward+backward; if the buffer lacks the minimum anchors (one orientation, one
  velocity, one gyro and one accel sample) it returns an **empty trajectory** instead of throwing,
  so callers skip integration for that scan.
- `ImuTransformer::process()` clamps `dt` to a default rate whenever it falls at or below
  `MIN_SANE_DT` (1 ms), not just when `<= 0`. Some IMU drivers emit near-duplicate timestamps
  (microseconds apart); dividing the finite-difference angular acceleration by such a near-zero
  `dt` amplifies a negligible angular-velocity delta into a huge spurious lever-arm spike.

## Full on-manifold preintegration + gravity-in-map estimation (2026-07)

Three additions, all GTSAM-free (MRPT Lie algebra only):

- `so3_jacobians.h` -- SO(3) right Jacobian `Jr` and its inverse, closed form
  with Taylor branches below 1e-4 rad. NOTE: `Jr` is NOT
  `mrpt::poses::Lie::SO<3>::jacob_dexpe_de()`, which returns the 9x3 Jacobian of
  the flattened rotation MATRIX w.r.t. the tangent vector: a different object.
- `ImuPreintegrator` -- real Forster preintegration: `dR`/`dV`/`dP`, the 9x9
  covariance (`[theta,v,p]` order), and the five first-order bias Jacobians, so
  a bias update needs no re-integration. `ImuIntegrator` stays as the cheap
  ROTATION-ONLY path (its name and docstring still oversell it, see
  `~/plans/800_imu.md`). Discretization is standard Euler-forward (reading held
  over the FOLLOWING dt, rotation taken at the step start), chosen over a
  midpoint rule because it keeps state, Jacobians and covariance mutually
  consistent; error is first order in dt, far below sensor noise at real rates.
  Inputs MUST already be body-frame (use `ImuTransformer`); an assertion
  enforces `sensor_pose` is unset or identity.
- `MapGravityEstimator` -- estimates the GRAVITY VECTOR IN THE MAP FRAME plus
  both IMU biases, over a sliding window of intervals, from preintegrated IMU
  aided by an external odometry's attitudes and velocities. Rearranging the
  delta definition, `g = (v_j - v_i - R_i dV_ij)/dt`, makes every interval a
  direct gravity measurement in which body acceleration CANCELS: unlike
  averaging accelerometers, it stays valid while moving arbitrarily. Unknowns
  are 9 scalars regardless of window length, so there are no attitude variables
  and no gauge to pin. `correction()` is the yaw-free rotation that levels the
  frame; `correction * R_odometry` is the gravity-consistent attitude.

Load-bearing details, do not "simplify" them away:
- The zero-anchored BIAS PRIORS and the |g| soft constraint are what keep the
  problem conditioned. Accel bias is only separable from a tilt given ROTATION
  EXCITATION; on a straight, level, constant-attitude run the two are
  indistinguishable and the prior is all that bounds the answer.
- The INTERVAL-COVERAGE GUARD (`min_interval_coverage`, default 0.95) rejects
  an interval whose IMU samples do not span it. A partial delta attributed to a
  whole interval corrupts the constraint SILENTLY; this exact bug cost a long
  debugging session in `mola_mapper`.
- Rejecting an interval is always safe here (the unknowns are global, so
  nothing is left unconstrained), unlike a full-state formulation which needs
  fallback priors.
- `ImuIntegrationParams` defaults `cov_gyro`/`cov_acc` to IDENTITY, i.e.
  sigma=1, ~1000x worse than any real MEMS IMU. Always set real noise
  densities, or every "earned" sigma comes out enormous and the estimator
  never reports convergence.

See `~/plans/801_lio_imu_preintegration_gravity.md` for the design and status.

## Build & test
Standard MOLA/colcon build (see root MOLA repo). Tests are plain executables that return non-zero
on failure and use MRPT `ASSERT_*` macros.

## Style
clang-format-14 (`scripts/formatter.sh`); American spelling; no one-line `if`; one variable per
declaration; anonymous namespaces over `static`; avoid em dashes.

## Active improvement plan
See `~/plans/800_imu.md` for the in-progress review (math gaps, missing `ImuTransformer` tests,
covariance/Jacobian TODOs). Keep that doc's checkboxes updated as work lands.
