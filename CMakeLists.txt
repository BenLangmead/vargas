cmake_minimum_required(VERSION 3.10)

project(vargas CXX C)

# Remove fixed output paths - let CMake use build directory
# set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
# set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

include_directories(doctest/doctest include cxxopts/src)

# --- External dependencies: vcflib + htslib (discovered via pkg-config) --------
# Both ship .pc files (system paths on Debian; under conda set
# PKG_CONFIG_PATH=$CONDA_PREFIX/lib/pkgconfig before invoking cmake).
find_package(PkgConfig REQUIRED)
pkg_check_modules(VCFLIB REQUIRED vcflib)
pkg_check_modules(HTS REQUIRED htslib)

# vcflib's own headers use bare includes (e.g. <Fasta.h>, <tabix.hpp>) while
# vargas includes them with the "vcflib/" prefix, so we need BOTH the directory
# vcflib.pc points at and its parent. The parent MUST come first: some installs
# ship two tabix.hpp (a working one in the top-level include dir and a broken one
# inside vcflib/), and the top-level parent dir is the one that resolves cleanly.
foreach(_vcfdir ${VCFLIB_INCLUDE_DIRS})
    include_directories(${_vcfdir}/.. ${_vcfdir})
endforeach()
include_directories(${HTS_INCLUDE_DIRS})
link_directories(${VCFLIB_LIBRARY_DIRS} ${HTS_LIBRARY_DIRS})

# tabixpp is a separate library on Debian but is bundled into libvcflib in newer
# (1.0.x) releases, which ship no standalone libtabixpp. Link it only if found.
find_library(TABIXPP_LIB tabixpp HINTS ${VCFLIB_LIBRARY_DIRS} ${HTS_LIBRARY_DIRS})

# Embed the resolved library dirs as an RPATH so the binary runs without the
# caller having to set LD_LIBRARY_PATH.
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_INSTALL_RPATH "${VCFLIB_LIBRARY_DIRS};${HTS_LIBRARY_DIRS}")

set(MAIN_SOURCES
        src/main.cpp
        src/graph.cpp
        src/utils.cpp
        src/varfile.cpp
        src/sim.cpp
        src/fasta.cpp
        src/sam.cpp
        src/align_main.cpp
        src/scoring.cpp
        src/graphman.cpp)

set(HEADERS
        include/alignment.h
        include/dyn_bitset.h
        include/fasta.h
        include/graph.h
        include/main.h
        include/sam.h
        include/sim.h
        include/utils.h
        include/varfile.h
        include/align_main.h
        include/scoring.h
        include/simd.h)

option(BUILD_AVX512BW_INTEL "Use Intel compiler to build for AVX512BW" OFF)
option(BUILD_AVX512BW_GCC "Use GCC compiler to build for AVX512BW" OFF)
option(BUILD_AVX2_INTEL "Use Intel compiler to build for AVX2" OFF)
option(BUILD_AVX2_GCC "Use GCC compiler to build for AVX2" OFF)

if(BUILD_AVX512BW_INTEL)
    message("   Building for AVX512BW (Intel)")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd597 -xCORE-AVX512 -DVA_SIMD_USE_AVX512")
elseif(BUILD_AVX512BW_GCC)
    message("   Building for AVX512BW (GCC)")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx512bw -DVA_SIMD_USE_AVX512")
elseif(BUILD_AVX2_INTEL)
    message("   Building for AVX2 (Intel)")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd597 -xCORE-AVX2 -DVA_SIMD_USE_AVX2")
elseif(BUILD_AVX2_GCC)
    message("   Building for AVX2 (GCC)")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -DVA_SIMD_USE_AVX2")
else()
    message("   Building for SSE4.1 (GCC)")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.1 -DVA_SIMD_USE_SSE")
endif()

add_executable(vargas ${MAIN_SOURCES})
target_link_libraries(vargas ${VCFLIB_LIBRARIES} ${HTS_LIBRARIES})
if(TABIXPP_LIB)
    target_link_libraries(vargas ${TABIXPP_LIB})
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -pthread")
# -fno-inline seems to be necessary to allow the release-mode test cases to pass
# TODO: investigate exactly where in the code that issue is
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG -fno-inline")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -pedantic -Wall -Wextra -Wno-unused-function -fno-inline")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS_DEBUG}")
