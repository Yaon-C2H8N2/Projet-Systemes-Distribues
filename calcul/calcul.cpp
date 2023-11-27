//mpicxx -std=c++11 -o nbody_simulation calcul.cpp
//time mpiexec -np 10 nbody_simulation 5000 500
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <mpi.h>
#include <random>
#include <iomanip>

const double G = 6.67430e-11; // Gravitational constant

struct Body {
    double mass;
    double x, y, z; // Coordinates
    double vx, vy, vz; // Velocities

    Body() : mass(0), x(0), y(0), z(0), vx(0), vy(0), vz(0) {}
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

void initBodies(std::vector<Body>& bodies) {
    // masse au centre
    if (!bodies.empty()) {
        bodies[0].mass = 100000.0;
        bodies[0].x = 0;
        bodies[0].y = 0;
        bodies[0].z = 0;
        bodies[0].vx = 0;
        bodies[0].vy = 0;
        bodies[0].vz = 0;
    }

    int externalDiameter = 80000;
    int internalDiameter = 40000;

    for (size_t i = 1; i < bodies.size(); ++i) {

        // random position en anneau en fonction de externalDiameter et internalDiameter
        double r = (internalDiameter+(rand() % (externalDiameter-internalDiameter)));
        double theta = (rand() % 360) * (M_PI/180);
        double x = r * cos(theta);
        double y = r * sin(theta);

        bodies[i].mass = 1.0; // plus petite masse
        bodies[i].x = x;
        bodies[i].y = y;
        bodies[i].z = 0;
        bodies[i].vx = 10;
        bodies[i].vy = -10;
        bodies[i].vz = 0;
    }
}

void writeOutput(std::ofstream& file, const std::vector<Body>& bodies, int timestep) {
    for (const auto& body : bodies) {
        file << timestep << "," << body.x << "," << body.y << "," << body.z << "," << body.mass << "\n";
    }
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double dt = 0.1;

    int total_bodies = 1000;
    int num_steps = 100;
    int buffer_size = 10;

    // parse cli
    if (argc >= 3) {
        total_bodies = std::atoi(argv[1]);
        num_steps = std::atoi(argv[2]);
        if (rank == 0)
            std::cout << "Bodies: " << total_bodies << " | Steps: " << num_steps << std::endl;
    } else {
        std::cout << "Usage: " << argv[0] << " <total_bodies> <num_steps>" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    buffer_size = num_steps / 10;

    std::vector<Body> bodies(total_bodies);
    std::vector<Body> localBodies(total_bodies / size);  // Subset of bodies for each process
    std::vector<std::vector<Body>> buffer(0);

    if (rank == 0) {
        initBodies(bodies);
    }

    MPI_Bcast(bodies.data(), total_bodies * sizeof(Body), MPI_BYTE, 0, MPI_COMM_WORLD);

    std::ofstream outputFile;
    if (rank == 0) {
        outputFile.open("../data/nbody_simulation.csv");
        outputFile << "step,x,y,z,mass\n";
    }

    int startIdx = rank * (total_bodies / size);
    int endIdx = startIdx + (total_bodies / size);

    for (int t = 0; t < num_steps; ++t) {
        // Each process calculates forces for all bodies
        calculateForces(bodies, rank, size);

        // Update velocities and positions for the subset of bodies
        for (int i = startIdx; i < endIdx; ++i) {
            bodies[i].x += bodies[i].vx * dt;
            bodies[i].y += bodies[i].vy * dt;
            bodies[i].z += bodies[i].vz * dt;
        }

        // Synchronize the updated data
        MPI_Allgather(&bodies[startIdx], (total_bodies / size) * sizeof(Body), MPI_BYTE,
                      bodies.data(), (total_bodies / size) * sizeof(Body), MPI_BYTE,
                      MPI_COMM_WORLD);

        if (rank == 0) {
            buffer.push_back(bodies);
            if (num_steps % 10 == 0)
            {
                std::cout << "10 steps done "<< t << std::endl;
            }
            if(buffer.size() >= buffer_size) {
                printf("Dumping steps %d to %d into csv file...\n", (t+1) - buffer_size, (t + 1));
                for(int i = 0; i < buffer_size; ++i) {
                    writeOutput(outputFile, buffer[i], t - buffer_size + i);
                }
                buffer = std::vector<std::vector<Body>>(0);
            }
        }
    }

    if (rank == 0) {
        outputFile.close();
    }

    MPI_Finalize();
    return 0;
}



