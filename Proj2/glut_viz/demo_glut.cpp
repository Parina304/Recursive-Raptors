#include <GL/freeglut.h>
using namespace std;
bool flag_fullscreen = false;

void render_cb(void){

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glBegin(GL_TRIANGLES);
    glVertex3f(-2, -2, .5);
    glVertex3f(2, 0, -5);
    glVertex3f(0, 2, -5);
    glEnd();
    glutSwapBuffers();
}

void changeSize(int w, int h){
    if(h == 0){
        h = 1;
    }
    float ratio = 1.0 * w / h;

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glViewport(0, 0, w, h);
    gluPerspective(45, ratio, 1, 100);
    glMatrixMode(GL_MODELVIEW);
}

void keyboard_cb(unsigned char key, int x, int y){
    switch (key)
    {
    case 27:
        exit(0);
        break;

    case 'q':
        exit(0);
        break;
    
    case 'f':
        glutFullScreenToggle();

    default:
        break;
    }
}

int main(void){
    int argc = 0;
    char* argv = {nullptr};
    int w = 1600, h = 900;
    glutInit(&argc, &argv);
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(w, h);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
    glutCreateWindow("Amazing GLUT");

    glutDisplayFunc(render_cb);
    glutReshapeFunc(changeSize);
    glutKeyboardFunc(keyboard_cb);

    
    glutMainLoop();
    return 0;
}