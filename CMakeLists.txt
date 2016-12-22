# Copyright (C) 2015   Michael Lahnert, Markus Baur, Jonas Harsch, Stephan Lunowa
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

project (NUMSIM)

cmake_minimum_required (VERSION 2.8)

# force to set a build type
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING
    "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# set compiler
set (CMAKE_CC_COMPILER "/usr/bin/mpicc")
set (CMAKE_CXX_COMPILER "/usr/bin/mpic++")

# set some good compiler flags
set (GOOD_FLAGS "-Wall -Wextra -pedantic -O3 -march=native -std=c++11")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GOOD_FLAGS}")

# debug visu?
option (DEBUG_VISU
  "Enable debugging visualization" ON)

if (DEBUG_VISU)
  # find sdl 2
  INCLUDE(FindPkgConfig)
  PKG_SEARCH_MODULE(SDL2 REQUIRED sdl2)
  INCLUDE_DIRECTORIES(${INCLUDE_DIRECTORIES} ${SDL2_INCLUDE_DIRS})
endif (DEBUG_VISU)

# write vtk?
option (WRITE_VTK
  "Writing out vtk files for use with paraview" ON)

# find mpi
find_package (MPI REQUIRED)
INCLUDE_DIRECTORIES(${INCLUDE_DIRECTORIES} ${MPI_INCLUDE_PATH})

# Using glob is discouraged in cmake.
set (HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/src/comm.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/compute.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/geometry.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/grid.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/iterator.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/parameter.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/solver.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/typedef.hpp
  )

set (SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/comm.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/compute.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/geometry.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/grid.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/iterator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/parameter.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/solver.cpp
  )

# add header and source of visu if it shall be built
# additionally set a flag
if (DEBUG_VISU)
  set (HEADERS
    ${HEADERS}
    ${CMAKE_CURRENT_SOURCE_DIR}/src/visu.hpp
    )

  set (SOURCES
    ${SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/src/visu.cpp
    )

  add_definitions(-DDEBUG_VISU)
endif (DEBUG_VISU)

# add header and source of vtk if it shall be built
# additionally set a flag
if (WRITE_VTK)
  set (HEADERS
    ${HEADERS}
    ${CMAKE_CURRENT_SOURCE_DIR}/src/vtk.hpp
    )

  set (SOURCES
    ${SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/src/vtk.cpp
    )

  add_definitions(-DWRITE_VTK)
endif (WRITE_VTK)

# build executable and link libs against it
add_executable (numsim ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)
add_executable (creator ${CMAKE_CURRENT_SOURCE_DIR}/src/creator.cpp)
add_executable (test_mpi ${CMAKE_CURRENT_SOURCE_DIR}/src/test_mpi.cpp)
add_executable (test_iter ${CMAKE_CURRENT_SOURCE_DIR}/src/test_iterator.cpp)
add_executable (test_sor ${CMAKE_CURRENT_SOURCE_DIR}/src/test_sor.cpp)
add_executable (test_rbsor ${CMAKE_CURRENT_SOURCE_DIR}/src/test_RBsor.cpp)
add_executable (test_dc ${CMAKE_CURRENT_SOURCE_DIR}/src/test_dc.cpp)
add_library (_numsim ${SOURCES})

target_link_libraries (numsim _numsim ${SDL2_LIBRARIES} ${MPI_LIBRARIES})
target_link_libraries (test_mpi _numsim ${SDL2_LIBRARIES} ${MPI_LIBRARIES})
target_link_libraries (test_iter _numsim ${SDL2_LIBRARIES} ${MPI_LIBRARIES})
target_link_libraries (test_sor _numsim ${SDL2_LIBRARIES} ${MPI_LIBRARIES})
target_link_libraries (test_rbsor _numsim ${SDL2_LIBRARIES} ${MPI_LIBRARIES})
target_link_libraries (test_dc _numsim ${SDL2_LIBRARIES} ${MPI_LIBRARIES})

if (MPI_COMPILE_FLAGS)
  set_target_properties(numsim PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if (MPI_LINK_FLAGS)
  set_target_properties (numsim PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
