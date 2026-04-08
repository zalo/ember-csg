cmake_minimum_required(VERSION 3.15...3.27)
project(ember-python LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Scikit-build-core sets these values
if (DEFINED SKBUILD)
  set(PYTHON_EXECUTABLE "${SKBUILD_PYTHON_EXECUTABLE}")
  set(PYTHON_INCLUDE_DIR "${SKBUILD_PYTHON_INCLUDE_DIR}")
  set(PYTHON_LIBRARY "${SKBUILD_PYTHON_LIBRARY}")
  set(CMAKE_INSTALL_PREFIX "${SKBUILD_PLATLIB_DIR}")
endif()

# Find Python
if (CMAKE_VERSION VERSION_LESS 3.18)
  set(DEV_MODULE Development)
else()
  set(DEV_MODULE Development.Module)
endif()

find_package(Python 3.8 COMPONENTS Interpreter ${DEV_MODULE} REQUIRED
               OPTIONAL_COMPONENTS Development.SABIModule)

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

# Nanobind — use the copy in VoCD/extern or a local extern
if (DEFINED SKBUILD)
  add_subdirectory(extern/nanobind)
else()
  set(nanobind_DIR "${CMAKE_CURRENT_SOURCE_DIR}/extern/nanobind/cmake")
  find_package(nanobind CONFIG REQUIRED)
endif()
include_directories(${nanobind_INCLUDE_DIRS})

# Enable BMI2 / LZCNT (required by typed-geometry fixed_int)
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-mbmi2" HAS_BMI2)
check_cxx_compiler_flag("-mlzcnt" HAS_LZCNT)
if(HAS_BMI2)
    add_compile_options(-mbmi2)
endif()
if(HAS_LZCNT)
    add_compile_options(-mlzcnt)
endif()
# Enable SIMD for tinybvh performance
add_compile_options(-msse4.2 -mavx -mavx2 -mfma)

# Dependencies from mesh-kernel
set(EMBER_CSG_DIR "${CMAKE_CURRENT_SOURCE_DIR}/..")
set(MESH_KERNEL_DIR "${EMBER_CSG_DIR}/../mesh-kernel")
set(TINYBVH_DIR "${EMBER_CSG_DIR}/../tinybvh")

if(NOT TARGET clean-core)
    add_subdirectory("${MESH_KERNEL_DIR}/extern/clean-core" "${CMAKE_BINARY_DIR}/clean-core")
endif()
if(NOT TARGET typed-geometry)
    add_subdirectory("${MESH_KERNEL_DIR}/extern/typed-geometry" "${CMAKE_BINARY_DIR}/typed-geometry")
endif()

# EMBER CSG static library
file(GLOB_RECURSE EMBER_SOURCES "${EMBER_CSG_DIR}/src/*.cc")
file(GLOB_RECURSE EMBER_HEADERS "${EMBER_CSG_DIR}/include/*.hh")

add_library(ember-csg STATIC ${EMBER_SOURCES} ${EMBER_HEADERS})
target_include_directories(ember-csg PUBLIC
    "${EMBER_CSG_DIR}/include"
    "${TINYBVH_DIR}"
    "${MESH_KERNEL_DIR}/src"
)
target_link_libraries(ember-csg PUBLIC typed-geometry clean-core)
target_compile_features(ember-csg PUBLIC cxx_std_20)

# Python module
nanobind_add_module(ember_ext STABLE_ABI LTO ember_ext.cpp)
target_link_libraries(ember_ext PRIVATE ember-csg)

install(TARGETS ember_ext LIBRARY DESTINATION ember)
