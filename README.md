# mola-imu-preintegration
Integrator of IMU angular velocity readings.

This repository provides:
* `IMUIntegrator` and `RotationIntegrator`: C++ classes to integrate IMU accelerations and angular velocities.
* `IMUIntegrationPublisher`: A MOLA module to process CObservationIMU readings and emit relative pose updates with respect to the last keyframe.

## Build and install
Refer to the [root MOLA repository](https://github.com/MOLAorg/mola).

## Docs and examples
See this package page [in the documentation](https://docs.mola-slam.org/latest/modules.html).

## License
This package is released under the GNU GPL v3 license.
