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
add_definitions(-Wall -Wextra -pedantic -O3 -march=native -std=c++11)

# debug visu?
option (DEBUG_VISU
  "Enable debugging visualization" ON)

# write vtk?
option (WRITE_VTK
  "Writing out vtk files for use with paraview" ON)

# output for exercise sheet 4?
option (BLATT4
  "Produce output for exercise sheet 4 (Blatt 4). Disables WRITE_VTK and DEBUG_VISU." OFF)

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

# set flag for exercise sheet 4
if (BLATT4)
  add_definitions(-DBLATT4)
else (BLATT4)
# add header and source of visu if it shall be built
# additionally set a flag
if (DEBUG_VISU)
  # find sdl 2
  INCLUDE(FindPkgConfig)
  PKG_SEARCH_MODULE(SDL2 REQUIRED sdl2)
  INCLUDE_DIRECTORIES(${INCLUDE_DIRECTORIES} ${SDL2_INCLUDE_DIRS})

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
endif (BLATT4)

# build executable and link libs against it
add_library (_numsim ${SOURCES})

add_executable (numsim ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)
target_link_libraries (numsim _numsim ${SDL2_LIBRARIES} ${MPI_LIBRARIES})

add_executable (compare_solvers ${CMAKE_CURRENT_SOURCE_DIR}/src/compare_solvers.cpp)
target_link_libraries (compare_solvers _numsim ${MPI_LIBRARIES})

add_executable (creator ${CMAKE_CURRENT_SOURCE_DIR}/src/creator.cpp)

if (CMAKE_BUILD_TYPE MATCHES Debug)
  message ("Building Debug.")
  add_executable (test_mpi ${CMAKE_CURRENT_SOURCE_DIR}/src/test_mpi.cpp)
  target_link_libraries (test_mpi _numsim ${SDL2_LIBRARIES} ${MPI_LIBRARIES})

  add_executable (test_iter ${CMAKE_CURRENT_SOURCE_DIR}/src/test_iterator.cpp)
  target_link_libraries (test_iter _numsim ${SDL2_LIBRARIES} ${MPI_LIBRARIES})

  add_executable (test_sor ${CMAKE_CURRENT_SOURCE_DIR}/src/test_sor.cpp)
  target_link_libraries (test_sor _numsim ${SDL2_LIBRARIES} ${MPI_LIBRARIES})

  add_executable (test_rbsor ${CMAKE_CURRENT_SOURCE_DIR}/src/test_RBsor.cpp)
  target_link_libraries (test_rbsor _numsim ${SDL2_LIBRARIES} ${MPI_LIBRARIES})

  add_executable (test_mg ${CMAKE_CURRENT_SOURCE_DIR}/src/test_MG.cpp)
  target_link_libraries (test_mg _numsim ${SDL2_LIBRARIES} ${MPI_LIBRARIES})

  add_executable (test_dc ${CMAKE_CURRENT_SOURCE_DIR}/src/test_dc.cpp)
  target_link_libraries (test_dc _numsim ${SDL2_LIBRARIES} ${MPI_LIBRARIES})
endif (CMAKE_BUILD_TYPE MATCHES Debug)
