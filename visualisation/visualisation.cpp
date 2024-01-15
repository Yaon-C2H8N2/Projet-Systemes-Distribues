#ifdef __APPLE__

#include <GLUT/glut.h> /* Pour Mac OS X */

#else
#include <GL/glut.h>
#endif

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <thread>


int NB_STEP = 300;
int NB_POINTS = 10000;

char presse;
int anglex, angley, x, y, xold, yold, zoom;
double ***points;

void affichage();

void reshape(int x, int y);

void mouse(int button, int state, int x, int y);

void readInputCSV();

int main(int argc, char **argv) {
    //Parsing CLI arguments
    if (argc >= 3) {
        NB_POINTS = std::atoi(argv[1]);
        NB_STEP = std::atoi(argv[2]);
    } else {
        std::cout << "Usage: " << argv[0] << " <total_bodies> <num_steps>" << std::endl;
        exit(1);
    }

    //Allocating memory for the points
    points = new double **[NB_STEP];
    for (int i = 0; i < NB_STEP; i++) {
        points[i] = new double *[NB_POINTS];
        for (int j = 0; j < NB_POINTS; j++) {
            points[i][j] = new double[2];
        }
    }

    //Starting file input thread to parallelize the reading of the csv
    std::thread reader(readInputCSV);
    //GLUT initialization and window creation
    zoom = 10000;
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowPosition(200, 200);
    glutInitWindowSize(500, 500);
    glutCreateWindow("N-Body Simulation");

    //OpenGL initialization
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glColor3f(1.0, 1.0, 1.0);
    glPointSize(5.0);
    glEnable(GL_DEPTH_TEST);

    //Callback function definitions
    glutDisplayFunc(affichage);
    glutReshapeFunc(reshape);
    glutMouseFunc(mouse);

    //Starting the main loop
    glutMainLoop();
    reader.join();

    //Freeing the memory
    for (int i = 0; i < NB_STEP; i++) {
        for (int j = 0; j < NB_POINTS; j++) {
            free(points[i][j]);
        }
        free(points[i]);
    }
    free(points);
    return 0;
}

int displayedStep = 0;

//Display function
void affichage() {
    //Clearing the window and setting shading model
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glShadeModel(GL_SMOOTH);

    //Camera position and orientation
    glLoadIdentity();
    glRotatef(0, 1.0, 0.0, 0.0);
    glRotatef(0, 0.0, 1.0, 0.0);
    glColor3f(1.0, 1.0, 1.0);

    //Showing each body as a point
    glBegin(GL_POINTS);
    glColor3f(1.0, 1.0, 1.0);
    //Dividing by 10000 to scale down the points spreading area
    for(int i = 0; i<NB_POINTS; i++){
        glVertex2d(points[displayedStep][i][0]/zoom, points[displayedStep][i][1]/zoom);
    }
    glEnd();

    glFlush();
    glutSwapBuffers();
    glutPostRedisplay();

    displayedStep = (displayedStep + 1) % NB_STEP;
}

//Function called when the window is resized
void reshape(int x, int y) {
    if (x < y)
        glViewport(0, (y - x) / 2, x, x);
    else
        glViewport((x - y) / 2, 0, y, y);
}

void mouse(int button, int state, int x, int y) {
    switch (button) {
        case 3:
            zoom = zoom - 1000;
            glutPostRedisplay();
            break;
        case 4:
            zoom = zoom + 1000;
            glutPostRedisplay();
            break;
        default:
            break;
    }
}

void readInputCSV() {
    //Filepath
    auto *inputFile = new std::ifstream("../data/nbody_simulation.csv");
    if (inputFile->is_open()) {
        std::cout << "Lecture des données..." << std::endl;
        std::string line;
        getline(*inputFile, line); // ligne meta
        //Reading step by step
        for (int step = 0; step < NB_STEP; step++) {
            //Reading body by body for each step
            for (int particule = 0; particule < NB_POINTS; particule++) {
                getline(*inputFile, line);
                //Debug var
                int stepCSV = 0;
                try {
                    stepCSV = stoi(line.substr(0, line.find(',')));
                } catch (std::invalid_argument &e) {
                    std::cout << "Error while parsing step : " << e.what() << std::endl;
                    std::cout << "step : " << step << " particule : " << particule << " detail : " << line << std::endl;
                    exit(1);
                }
                //Parsing the line
                line = line.substr(line.find(',') + 1);
                double xCSV = stod(line.substr(0, line.find(',')));
                line = line.substr(line.find(',') + 1);
                double yCSV = stod(line.substr(0, line.find(',')));

                //Storing the parsed data
                points[step][particule][0] = xCSV;
                points[step][particule][1] = yCSV;
            }
        }
        inputFile->close();
        std::cout << "Fin de la lecture des données" << std::endl;
    }
}