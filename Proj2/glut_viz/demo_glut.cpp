#include <iostream>
#include <GL/freeglut.h>
#include <cppcolormap.h>
#include <xtensor/xio.hpp>

using namespace std;
bool flag_fullscreen = false;

#define RED 1
#define GREEN 2
#define BLUE 3
#define ORANGE 4

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

void ProcessMenuEvents(int option){
    
    cout << option << endl;
}

void CreateMenu(){
    int menu;
    menu = glutCreateMenu(ProcessMenuEvents);

    glutAddMenuEntry("red", RED);
    glutAddMenuEntry("green", GREEN);
    glutAddMenuEntry("blue", BLUE);
    glutAddMenuEntry("orange", ORANGE);

    glutAttachMenu(GLUT_RIGHT_BUTTON);

}

int main(void){
    int argc = 0;
    char* argv = {nullptr};
    int w = 1600, h = 900;
    // glutInit(&argc, &argv);
    // glutInitWindowPosition(0, 0);
    // glutInitWindowSize(w, h);
    // glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
    // glutCreateWindow("Amazing GLUT");

    // glutDisplayFunc(render_cb);
    // glutReshapeFunc(changeSize);
    // glutKeyboardFunc(keyboard_cb);

    
    // glutMainLoop();
    xt::xtensor<double, 1> tensor = {0.4};

    auto map = cppcolormap::colormap("Reds", 10);
    auto c = cppcolormap::as_colors(tensor, cppcolormap::colormap("jet"), 0, 1);
    // for (auto& f: map){
    //     std::cout << f << " ";
    //     // for (auto& d: f){
    //     //     std::cout << typeid(d).name() << std::endl;
    //     // }
    // }
    std::cout << c;
    // std::cout << cppcolormap::colorcycle("tue", 256) << std::endl;
    return 0;
}