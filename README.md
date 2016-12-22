Numerical solver for the Navier-Stokes equations
================================================

Requirements
------------

* C++ compiler, which at least supports C++ standard 2011
* MPI implementation (openMPI / MPICH)
* SDL2 library headers (libsdl2-dev)

How To Compile
--------------

Run the following instructions inside the main folder:

    cmake .
    make

If You wish to disable real-time visualisation, run instead:

    cmake . -DDEBUG_VISU=off
    make

Usage
-----

The program *creator* produces geometry and parameter files.
The program *numsim* solves the Navier-Stokes equations on a given geometry
with the given parameters. It is parallelized using MPI, such that
`mpirun numsim` will start the program with several threads.

For creating geometry and parameter files please look at `creator -h`.
For running the solver please look at `numsim -h`.

### Example

To use the fluid simulation for a Karman vortex street of real length 6x1 with
180x30 grid cells (written to kar30.geom and kar30.param), one runs:

    ./creator -o kar30 -pre 4 -length 6x1 -size 180x30
    ./numsim -I kar30

The vtk-output is then inside the directory `VTK/`.

