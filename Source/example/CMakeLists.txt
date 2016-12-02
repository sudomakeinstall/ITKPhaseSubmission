# This is the root ITK CMakeLists file.
cmake_minimum_required(VERSION 2.4)
if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif()

# This project is designed to be built outside the Insight source tree.
project(dv-simulated_phase)

# Find ITK.
find_package(ITK 4.7 REQUIRED)
include(${ITK_USE_FILE})

add_executable(dv-simulated_phase dv-simulated_phase.cxx )

target_link_libraries(dv-simulated_phase ${ITK_LIBRARIES})
