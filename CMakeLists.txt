cmake_minimum_required(VERSION 3.14)
project(EmberCSG LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable BMI2 for fast integer operations (required by typed-geometry fixed_int)
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-mbmi2" HAS_BMI2)
check_cxx_compiler_flag("-mlzcnt" HAS_LZCNT)
if(HAS_BMI2)
    add_compile_options(-mbmi2)
endif()
if(HAS_LZCNT)
    add_compile_options(-mlzcnt)
endif()

# Dependencies from mesh-kernel
set(MESH_KERNEL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../mesh-kernel")

if(NOT TARGET clean-core)
    add_subdirectory("${MESH_KERNEL_DIR}/extern/clean-core" "${CMAKE_BINARY_DIR}/clean-core")
endif()
if(NOT TARGET typed-geometry)
    add_subdirectory("${MESH_KERNEL_DIR}/extern/typed-geometry" "${CMAKE_BINARY_DIR}/typed-geometry")
endif()

# EMBER CSG library
file(GLOB_RECURSE EMBER_HEADERS "include/*.hh")
file(GLOB_RECURSE EMBER_SOURCES "src/*.cc")

add_library(ember-csg STATIC ${EMBER_SOURCES} ${EMBER_HEADERS})
target_include_directories(ember-csg PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${MESH_KERNEL_DIR}/src"
)
target_link_libraries(ember-csg PUBLIC typed-geometry clean-core)
target_compile_features(ember-csg PUBLIC cxx_std_20)

# Test executables
add_executable(ember-test test/test_cubes.cc)
target_link_libraries(ember-test PRIVATE ember-csg)

add_executable(ember-roundtrip test/test_roundtrip.cc)
target_link_libraries(ember-roundtrip PRIVATE ember-csg)



add_executable(ember-perf test/test_perf.cc)
target_link_libraries(ember-perf PRIVATE ember-csg)

add_executable(ember-manifold test/test_manifold.cc)
target_link_libraries(ember-manifold PRIVATE ember-csg)

add_executable(ember-bvh-validate test/test_bvh_validate.cc)
target_link_libraries(ember-bvh-validate PRIVATE ember-csg)
