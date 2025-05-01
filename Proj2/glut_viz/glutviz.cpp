#ifdef _WIN32
#include <freeglut.h>
#else
#include <GL/freeglut.h>
#endif

#include "mesh_lib.h"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include "cppcolormap.h"
// #include "ColormapWrapper.h"
#include <xtensor/xio.hpp>

bool loadSecondFile = false;  // Flag to check if user wants to load 2 files

float angleX = 0.0f, angleY = 0.0f, zoom = 1.0f, panX = 0.0f, panY = 0.0f;
float cutPlaneZ = 0.0f;
float minZ = 0.0f, maxZ = 0.0f;
int window_width = 1600, window_height = 900;
// Globals for mouse interaction
bool flag_left_mouse_pressed = false;
bool flag_right_mouse_pressed = false;
int last_mouse_x = 0, last_mouse_y = 0;
bool flag_color_calced = false;

Mesh humanoid;

void DrawMesh(const Mesh& m, float translateX = 0.0f){
    glColor3f(0.0f, 1.0f, 1.0f);
    glPushMatrix();
    glTranslatef(translateX, 0.0f, 0.0f);
    float z_mean;
    xt::xtensor<double, 1> color_idx;
    
    for (const auto& f: m.faces) {
        glColor3f(f.r, f.g, f.b);
        glBegin(GL_POLYGON);
        for (const auto& idx: f.vertexIndices){
            glNormal3f(m.vertices[idx].x, m.vertices[idx].y, m.vertices[idx].z);
            glVertex3f(m.vertices[idx].x, m.vertices[idx].y, m.vertices[idx].z);
        }
        glEnd();
    }
    glPopMatrix();
}

void CalcMeshColors(Mesh& m){
    float y_mean = 0;
    xt::xtensor<double, 1> color_idx;
    // Color c;

    for (auto& f: m.faces){
        y_mean = 0;
        for (const auto& idx: f.vertexIndices){
            y_mean += m.vertices[idx].y;
        }
        y_mean /= f.vertexIndices.size();

        color_idx = {y_mean};
        auto c = cppcolormap::as_colors(color_idx, cppcolormap::viridis(), m.y_min, m.y_max);
        f.r = static_cast<float>(c(0, 0));
        f.g = static_cast<float>(c(0, 1));
        f.b = static_cast<float>(c(0, 2));
    }
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    glTranslatef(panX, panY, -5.0f);
    glScalef(zoom, zoom, zoom);
    glRotatef(angleX, 1, 0, 0);
    glRotatef(angleY, 0, 1, 0);
    // drawCuttingPlane();
    float color1[3] = { 0.0f, 1.0f, 1.0f };
    DrawMesh(humanoid);
    glutSwapBuffers();
}
bool wireframe = false;

void keyboard(unsigned char key, int x, int y) {
    switch (key) {
    case 'm': wireframe = !wireframe;
        if (wireframe)
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        else
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        break;
    case 'w': angleX -= 5.0f; break;
    case 's': angleX += 5.0f; break;
    case 'a': angleY -= 5.0f; break; 
    case 'd': angleY += 5.0f; break;
    case '+': zoom += 0.1f; break;
    case '-': zoom -= 0.1f; break;
    case 'j': panX -= 0.1f; break;  // Left
    case 'l': panX += 0.1f; break;  // Right
    case 'i': panY += 0.1f; break;  // Up
    case 'k': panY -= 0.1f; break;  // Down
    case 'z': cutPlaneZ += 0.1f; break;
    case 'x': cutPlaneZ -= 0.1f; break;
    
    case 'f': glutFullScreenToggle(); break;
    case 27: 
    case 'q': exit(0);
    }
    glutPostRedisplay();
}

void reshape(int width, int height) {
    // Prevent division by zero
    if (height == 0) height = 1;

    // Update the viewport to cover the new window dimensions
    glViewport(0, 0, width, height);

    // Update the projection matrix
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0, (double)width / (double)height, 1.0, 100.0);
    glMatrixMode(GL_MODELVIEW);
}

void initOpenGL() {
    glEnable(GL_DEPTH_TEST);

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);

    // Set light properties
    GLfloat ambient[] = { 0.5f, 0.5f, 0.5f, 1.0f };  // Ambient light
    GLfloat diffuse[] = { 0.05f, 0.05, 0.05f, 1.0f };  // Diffuse light
    GLfloat specular[] = { .0f, .0f, .0f, 1.0f }; // Specular light
    GLfloat position[] = { 1.0f, 1.0f, 1.0f, 0.0f }; // Light position

    glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
    glLightfv(GL_LIGHT0, GL_POSITION, position);

    // Enable color material mode
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

    // Set material properties
    GLfloat mat_specular[] = { .0f, .0f, .0f, 1.0f };
    GLfloat mat_shininess[] = {.0f }; // Shininess factor
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0, 1.0, 1.0, 100.0);
    glMatrixMode(GL_MODELVIEW);
}

void mouse(int button, int state, int x, int y) {
    if (button == GLUT_LEFT_BUTTON) {
        flag_left_mouse_pressed = (state == GLUT_DOWN);
    } else if (button == GLUT_RIGHT_BUTTON) {
        flag_right_mouse_pressed = (state == GLUT_DOWN);
    } else if (button == 3){
    // std::cout << "scrolled\n";
        if(state == GLUT_DOWN){
            zoom += .1f;
        }
    } else if (button == 4 && state == GLUT_DOWN){
        zoom -= .1f;
    }
    last_mouse_x = x;
    last_mouse_y = y;

    glutPostRedisplay();
}

void motion(int x, int y) {
    int dx = x - last_mouse_x;
    int dy = y - last_mouse_y;

    if (flag_left_mouse_pressed) {
        // Rotate the model
        angleX += dy * 0.5f;
        angleY += dx * 0.5f;
    } else if (flag_right_mouse_pressed) {
        // Pan the model
        panX += dx * 0.01f;
        panY -= dy * 0.01f;
    }

    last_mouse_x = x;
    last_mouse_y = y;

    glutPostRedisplay();
}

int main(int argc, char** argv) {
    int fake_argc = 1;
    char* fake_argv[] = { argv[0], nullptr };
    int fileCount = 1;
    
    std::string file1 = "../../assets/obj/humanoid_robot.obj";
    humanoid.loadOBJ(file1);
    humanoid.printMeshStats();
    if (!flag_color_calced){
        CalcMeshColors(humanoid);
        flag_color_calced = true;
    }

    glutInit(&fake_argc, fake_argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(window_width, window_height);
    
    glutCreateWindow("OBJ Viewer");
    initOpenGL();
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    // glutMouseWheelFunc(mouseWheel);  // Requires freeglut
    glutMainLoop();
    return 0;
}