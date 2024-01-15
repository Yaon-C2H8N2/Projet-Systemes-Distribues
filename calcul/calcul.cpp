//mpicxx -std=c++11 -o nbody_simulation calcul.cpp
//time mpiexec -np 10 nbody_simulation 5000 500
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <mpi.h>
#include <random>
#include <iomanip>
#include <thread>

const double G = 6.67430e-11; // Gravitational constant

struct Body {
    double mass;
    double x, y, z; // Coordinates
    double vx, vy, vz; // Velocities

    Body() : mass(0), x(0), y(0), z(0), vx(0), vy(0), vz(0) {}
    Body(double m, double _x, double _y, double _z, double _vx, double _vy, double _vz)
            : mass(m), x(_x), y(_y), z(_z), vx(_vx), vy(_vy), vz(_vz) {}
};

//Calculates the forces between the calling node's subset of bodies and all other bodies
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

//Updates the positions of the bodies
void updatePositions(std::vector<Body>& bodies, double dt, int rank, int numProcesses) {
    int numBodies = bodies.size();
    int bodiesPerProcess = numBodies / numProcesses;
    int startIdx = rank * bodiesPerProcess;
    int endIdx = startIdx + bodiesPerProcess;

    for (int i = startIdx; i < endIdx; ++i) {
        bodies[i].x += bodies[i].vx * dt;
        bodies[i].y += bodies[i].vy * dt;
        //bodies[i].z += bodies[i].vz * dt;
    }
}

//Initializes the bodies with random positions
void initBodies(std::vector<Body>& bodies) {
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

        //random position en anneau en fonction de externalDiameter et internalDiameter
        double r = (internalDiameter+(rand() % (externalDiameter-internalDiameter)));
        double theta = (rand() % 360) * (M_PI/180);
        double x = r * cos(theta);
        double y = r * sin(theta);

        bodies[i].mass = 1.0;
        bodies[i].x = x;
        bodies[i].y = y;
        bodies[i].z = 0;
        bodies[i].vx = 10;
        bodies[i].vy = -10;
        bodies[i].vz = 0;
    }
}

//Writes a single step into the output file
void writeOutput(std::ofstream& file, const std::vector<Body>& bodies, int timestep) {
    for (const auto& body : bodies) {
        file << timestep << "," << body.x << "," << body.y << "," << body.z << "," << body.mass << "\n";
    }
}

//Dumps the steps from the buffer into the output file
void dumpSteps(const std::vector<std::vector<Body>> buffer, int buffer_size, int t, std::ofstream& outputFile){
    printf("Dumping steps %d to %d into csv file...\n", (t+1) - buffer_size, (t + 1));
    for(int i = 0; i < buffer_size; ++i) {
        writeOutput(outputFile, buffer[i], t - buffer_size + i + 1);
    }
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size, hostname_length;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Get_processor_name(hostname, &hostname_length);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

    //Timescale of a single step
    double dt = 0.1;

    //Default settings
    int total_bodies = 1000;
    int num_steps = 100;
    int buffer_size = 10;

    //Parsing CLI arguments
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

    //Creating the initial bodies
    std::vector<Body> bodies(total_bodies);
    if (rank == 0) {
        initBodies(bodies);
    }

    //Broadcast the initial bodies to all nodes
    MPI_Bcast(bodies.data(), total_bodies * sizeof(Body), MPI_BYTE, 0, MPI_COMM_WORLD);

    //Output file initialization
    std::ofstream outputFile;
    if (rank == 0) {
        outputFile.open("../data/nbody_simulation.csv");
        outputFile << "step,x,y,z,mass\n";
    }

    //Initializing the output buffer and thread to limit the number of write operations
    std::vector<std::vector<Body>> buffer(0);
    std::thread dumper;

    //Calculating the start index of the subset of bodies for each node
    int startIdx = rank * (total_bodies / size);

    for (int t = 0; t < num_steps; ++t) {
        //Each node calculates forces for all bodies
        calculateForces(bodies, rank, size);

        //Update velocities and positions for the subset of bodies
        updatePositions(bodies, dt, rank, size);

        //Synchronize the updated data between all nodes
        MPI_Allgather(&bodies[startIdx], (total_bodies / size) * sizeof(Body), MPI_BYTE,
                      bodies.data(), (total_bodies / size) * sizeof(Body), MPI_BYTE,
                      MPI_COMM_WORLD);

        //Root node specific operations
        if (rank == 0) {
            //Push the updated data to the buffer
            buffer.push_back(bodies);
            //Print progress
            if (t % 10 == 0)
            {
                std::cout << "10 steps done "<< t << std::endl;
            }
            //Stop the dumper thread if it is running
            if(dumper.joinable()){
                printf("Joining dumper thread...\n");
                dumper.join();
            }
            //Start the dumper thread if the buffer is full
            if(buffer.size() >= buffer_size) {
                printf("Starting dumper thread...\n");
                dumper = std::thread(dumpSteps, buffer, buffer_size, t, std::ref(outputFile));
                buffer = std::vector<std::vector<Body>>(0);
            }
        }
    }

    //Finalize the dumper thread
    if (rank == 0) {
        outputFile.close();
        printf("Joining dumper thread...\n");
        dumper.join();
    }

    MPI_Finalize();
    return 0;
}



