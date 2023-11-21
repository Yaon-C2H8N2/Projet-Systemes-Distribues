#ifdef __APPLE__

#include <GLUT/glut.h> /* Pour Mac OS X */

#else
#include <GL/glut.h>
#endif

#include <cmath>
#include <cstdio>

char presse;
int anglex, angley, x, y, xold, yold;

void affichage();

void reshape(int x, int y);

int main(int argc, char **argv) {
    //init glut et creation de la fenetre
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowPosition(200, 200);
    glutInitWindowSize(500, 500);
    glutCreateWindow("N-Body Simulation");

    //init opengl
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glColor3f(1.0, 1.0, 1.0);
    glPointSize(3.0);
    glEnable(GL_DEPTH_TEST);

    //def des fonctions de callbacks
    glutDisplayFunc(affichage);
    glutReshapeFunc(reshape);

    glutMainLoop();
    return 0;
}

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
    glEnd();

    glFlush();
    glutSwapBuffers();
    glutPostRedisplay();
}

void reshape(int x, int y) {
    if (x < y)
        glViewport(0, (y - x) / 2, x, x);
    else
        glViewport((x - y) / 2, 0, y, y);
}