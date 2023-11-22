//mpicxx -o nbody_simulation calcul.cpp
//time mpiexec -n 1 nbody_simulation
#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <mpi.h>

const double G = 6.67430e-11; // Gravitational constant

struct Body {
    double mass;
    double x, y, z; // Coordinates
    double vx, vy, vz; // Velocities

    Body(double m, double _x, double _y, double _z, double _vx, double _vy, double _vz)
            : mass(m), x(_x), y(_y), z(_z), vx(_vx), vy(_vy), vz(_vz) {}
};

void calculateForces(std::vector<Body>& bodies, int rank, int numProcesses) {
    const long massMultiplier = 1e12;
    int numBodies = bodies.size();
    int bodiesPerProcess = numBodies / numProcesses;
    int startIdx = rank * bodiesPerProcess;
    int endIdx = startIdx + bodiesPerProcess;

    for (int i = startIdx; i < endIdx; ++i) {
        for (int j = 0; j < numBodies; ++j) {
            if (i != j) {
                double dx = bodies[j].x - bodies[i].x;
                double dy = bodies[j].y - bodies[i].y;
                //double dz = bodies[j].z - bodies[i].z;
                //double dist = sqrt(dx * dx + dy * dy + dz * dz);
                double dist = sqrt(dx * dx + dy * dy );
                double epsilon = 1e-10;
                double force = (G * bodies[i].mass * massMultiplier * bodies[j].mass * massMultiplier) / (dist * dist * dist + epsilon);

                //double force = (G * bodies[i].mass * massMultiplier* bodies[j].mass * massMultiplier) / (dist * dist * dist);

                double fx = force * dx;
                double fy = force * dy;
                //double fz = force * dz;

                bodies[i].vx += fx / (bodies[i].mass * massMultiplier);
                bodies[i].vy += fy / (bodies[i].mass * massMultiplier);
                //bodies[i].vz += fz / (bodies[i].mass * massMultiplier);
            }
        }
    }
}

void updatePositions(std::vector<Body>& bodies, double dt) {
    for (auto& body : bodies) {
        body.x += body.vx * dt;
        body.y += body.vy * dt;
        //body.z += body.vz * dt;
    }
}

int main(int argc, char** argv) {
    srand(static_cast<unsigned>(time(nullptr)));

    int rank, numProcesses;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);

    std::vector<Body> bodies;

    bodies.push_back(Body(100, (rand() % 100) - 50, (rand() % 100) - 50, (rand() % 100) - 50, 0,0,0));


    for (int i = 0; i<500; i++){
        bodies.push_back(Body(1, (rand() % 100) - 50, (rand() % 100) - 50, (rand() % 100) - 50, 0,0,0));

    }

    //bodies.push_back(Body(10000000, 0, 0, 0, 0, 0, 0));
    //bodies.push_back(Body(1000000, 9500, -10000, 0, -1000, 1000, 0));
    //bodies.push_back(Body(100000, -10100, 9000, 0, 1000, -1000, 0));



    double dt = 0.01; // Time step

    std::ofstream csvFile("../data/nbody_simulation.csv");
    if (rank == 0) {
        csvFile << "step,x,y,z,mass\n";
    }

    for (int step = 0; step < 10000; ++step) {
        calculateForces(bodies, rank, numProcesses);

        // Gather and update positions
        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, bodies.data(), sizeof(Body) * bodies.size(), MPI_BYTE, MPI_COMM_WORLD);
        updatePositions(bodies, dt);

        if (rank == 0) {
            for (size_t i = 0; i < bodies.size(); ++i) {
                csvFile << step << "," << bodies[i].x << "," << bodies[i].y << "," << bodies[i].z << "," << bodies[i].mass << "\n";
            }
        }
    }

    csvFile.close();

    MPI_Finalize();

    return 0;
}
