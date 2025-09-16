# Unit tests:
mola_add_test(
  TARGET  test-imu-rotation-integrator
  SOURCES test-imu-rotation-integrator.cpp
  LINK_LIBRARIES
    mola::mola_imu_preintegration
)

mola_add_test(
  TARGET  test-velocity-buffer
  SOURCES test-velocity-buffer.cpp
  LINK_LIBRARIES
    mola::mola_imu_preintegration
)
