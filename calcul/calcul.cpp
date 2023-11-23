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
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis_position(-50.0, 50.0); // Range for position
    std::uniform_real_distribution<> dis_velocity(-0.1, 0.1); // Range for velocity

    // masse au centre
    if (!bodies.empty()) {
        bodies[0].mass = 10000.0;
        bodies[0].x = 0;
        bodies[0].y = 0;
        bodies[0].z = 0;
        bodies[0].vx = 0;
        bodies[0].vy = 0;
        bodies[0].vz = 0;
    }


    for (size_t i = 1; i < bodies.size(); ++i) {
        bodies[i].mass = 1.0; // plus petite masse
        bodies[i].x = dis_position(gen);
        bodies[i].y = dis_position(gen);
        bodies[i].z = dis_position(gen);
        bodies[i].vx = dis_velocity(gen);
        bodies[i].vy = dis_velocity(gen);
        bodies[i].vz = dis_velocity(gen);
    }
}


int main(int argc, char** argv) {



    MPI_Init(&argc, &argv);



    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // params
    const double timeStep = 0.01;

    int total_bodies = 1000;
    int num_steps = 100;

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
    // distribuer
    int local_num_bodies = total_bodies / size;
    std::vector<Body> local_bodies(local_num_bodies);

    // local boides
    initBodies(local_bodies);

    std::ofstream csvFile;
    if (rank == 0) {
        csvFile.open("../data/nbody_simulation.csv");
        csvFile << "step,x,y,z,mass\n";
    }

    for (int step = 0; step < num_steps; ++step) {
        calculateForces(local_bodies, rank, size);
        updatePositions(local_bodies, timeStep);

        std::vector<Body> all_bodies;
        if (rank == 0) {
            all_bodies.resize(total_bodies);
        }

        MPI_Gather(local_bodies.data(), local_num_bodies * sizeof(Body), MPI_BYTE,
                   all_bodies.data(), local_num_bodies * sizeof(Body), MPI_BYTE, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            for (const auto& body : all_bodies) {
                csvFile << step << ","
                        << std::fixed << std::setprecision(6)
                        << body.x << ","
                        << body.y << ","
                        << body.z << ","
                        << body.mass << "\n";
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    if (rank == 0) {
        csvFile.close();
    }

    MPI_Finalize();
    return 0;
}
