# N Body distributed simulation

This project consists of two main parts: a simulation program (`calcul.cpp`) and a visualization program (`visualisation.cpp`). The simulation program calculates the positions of bodies in a gravitational field and writes the results to a CSV file. The visualization program reads the CSV file and displays the positions of the bodies.

## Speedup

The following table shows the speedup achieved by parallelizing the simulation program using MPI. The simulation was run on a cluster of 18 nodes, each with Intel I9-10900k 10 cores @3.70GHz processors for a total of 180 cores.

| Number of bodies | Number of steps | Number of cores | Time   |
|------------------|-----------------|-----------------|--------|
| 10000            | 300             | 1               | 834,45s |
| 10000            | 300             | 10              | 86,61s |
| 10000            | 300             | 20              | 54,96s |
| 10000            | 300             | 180             | 24,1s  |

## Prerequisites

- C++ compiler (supporting C++11)
- MPI (Message Passing Interface) for parallel computing
- OpenGL, GLUT libraries for visualization
- CMake (version 3.26 or higher)

## Compilation and Execution

### Simulation Program

The simulation program is parallelized using MPI and can be compiled and run using the `mpicxx` compiler and `mpiexec` command.

To compile the program, navigate to the directory containing `calcul.cpp` and run the following command:

```bash
mpicxx -std=c++11 -o nbody_simulation calcul.cpp
```

This will create an executable named `nbody_simulation`.

To run the program, use the `mpiexec` command. The `-hostfile` option specifies the file containing the host list. The program takes two arguments: the total number of bodies and the number of steps to simulate. For example, to simulate 5000 bodies for 500 steps, run:

```bash
mpiexec -hostfile <file path> ./nbody_simulation 5000 500
```

This will create a CSV file named `nbody_simulation.csv` in the `data` directory.

### Visualization Program

The visualization program uses OpenGL and GLUT for graphics. It can be compiled using CMake.

First, navigate to the project root directory and create a new directory for the build:

```bash
mkdir build
cd build
```

Then, run CMake to generate the Makefile:

```bash
cmake ..
```

Next, compile the program using the generated Makefile:

```bash
make
```

This will create an executable named `Projet_Systemes_Distribues_visualisation`.

To run the visualization program for 5000 bodies and 500 steps, simply execute the following command:

```bash
./Projet_Systemes_Distribues_visualisation 5000 500
```

This will open a window displaying the positions of the bodies as calculated by the simulation program.