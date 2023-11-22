#ifdef __APPLE__

#include <GLUT/glut.h> /* Pour Mac OS X */

#else
#include <GL/glut.h>
#endif

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>

#define NB_STEP 10000
#define NB_POINTS 2


char presse;
int anglex, angley, x, y, xold, yold;
double points[NB_STEP][NB_POINTS][2];

void affichage();

void reshape(int x, int y);

void readInputCSV();

int main(int argc, char **argv) {
    std::cout << "Lecture des données..." << std::endl;
    readInputCSV();
    //todo : print le load time
    std::cout << "Fin de la lecture des données" << std::endl;
    //init glut et creation de la fenetre
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowPosition(200, 200);
    glutInitWindowSize(500, 500);
    glutCreateWindow("N-Body Simulation");

    //init opengl
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glColor3f(1.0, 1.0, 1.0);
    glPointSize(5.0);
    glEnable(GL_DEPTH_TEST);

    //def des fonctions de callbacks
    glutDisplayFunc(affichage);
    glutReshapeFunc(reshape);

    glutMainLoop();
    return 0;
}

int displayedStep = 0;
void affichage() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glShadeModel(GL_SMOOTH);

    glLoadIdentity();
    glRotatef(0, 1.0, 0.0, 0.0);
    glRotatef(0, 0.0, 1.0, 0.0);
    glColor3f(1.0, 1.0, 1.0);

    //lecture des prochaines positions des points depuis le fichier csv
    glBegin(GL_POINTS);
    //affichage des points avec glVertex2d(x,y)
    glColor3f(1.0, 1.0, 1.0);
    //todo: remplacer le 5000 par la coord x et y la plus grande pour que toutes les valeurs soient entre 0 et 1
    for(int i = 0; i<NB_POINTS; i++){
        glVertex2d(points[displayedStep][i][0]/10000, points[displayedStep][i][1]/10000);
    }
    glEnd();

    glFlush();
    glutSwapBuffers();
    glutPostRedisplay();

    displayedStep = (displayedStep+1)%NB_STEP;
}

void reshape(int x, int y) {
    if (x < y)
        glViewport(0, (y - x) / 2, x, x);
    else
        glViewport((x - y) / 2, 0, y, y);
}

void readInputCSV() {
    auto *inputFile = new std::ifstream("../data/nbody_simulation.csv");
    if(inputFile->is_open()){
        std::string line;
        getline(*inputFile, line); // ligne meta
        for(int step = 0; step<NB_STEP; step++){
            for(int particule = 0; particule<NB_POINTS; particule++){
                getline(*inputFile, line);
                //std::cout << "step : " << step << " particule : " << particule << " detail : " << line << std::endl;
                //todo parse la ligne et remplir le tableau
                int stepCSV = stoi(line.substr(0, line.find(',')));
                line = line.substr(line.find(',')+1);
                double xCSV = stod(line.substr(0, line.find(',')));
                line = line.substr(line.find(',')+1);
                double yCSV = stod(line.substr(0, line.find(',')));
                //std::cout << "parsedStep : " << stepCSV << " parsedX : " << xCSV << " parsedY : " << yCSV << std::endl;
                //todo : si step = stepcsv alors remplir le tableau sinon exit tout avec erreur de lecture du csv
                points[step][particule][0] = xCSV;
                points[step][particule][1] = yCSV;
            }
        }
        inputFile->close();
    }
}