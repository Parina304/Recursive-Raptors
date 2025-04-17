#ifdef _WIN32
#include <freeglut.h>
#else
#include <GL/freeglut.h>
#endif

#include "mesh_lib.h"
#include "mesh_merge.h"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>

bool loadSecondFile = false;  // Flag to check if user wants to load 2 files

float angleX = 0.0f, angleY = 0.0f, zoom = 1.0f, panX = 0.0f, panY = 0.0f;
float cutPlaneZ = 0.0f;
float minZ = 0.0f, maxZ = 0.0f;
int window_width = 1600, window_height = 900;
// Globals for mouse interaction
bool flag_left_mouse_pressed = false;
bool flag_right_mouse_pressed = false;
int last_mouse_x = 0, last_mouse_y = 0;

Mesh humanoid;

void computeZBounds(const std::vector<float>& vertices) {
    if (vertices.empty()) return;
    minZ = maxZ = vertices[2];
    for (size_t i = 2; i < vertices.size(); i += 3) {
        if (vertices[i] < minZ) minZ = vertices[i];
        if (vertices[i] > maxZ) maxZ = vertices[i];
    }
    cutPlaneZ = (minZ + maxZ) / 2.0f;
}

void drawCuttingPlane() {
    float planeSize = 2.0f;
    glDisable(GL_LIGHTING);
    glColor4f(0.0f, 0.0f, 1.0f, 0.3f);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glBegin(GL_QUADS);
    glVertex3f(-planeSize, -planeSize, cutPlaneZ);
    glVertex3f(planeSize, -planeSize, cutPlaneZ);
    glVertex3f(planeSize, planeSize, cutPlaneZ);
    glVertex3f(-planeSize, planeSize, cutPlaneZ);
    glEnd();
    glDisable(GL_BLEND);
    glEnable(GL_LIGHTING);
}

void drawMesh(const std::vector<float>& vertices, const std::vector<int>& faces, float color[3], float translateX = 0.0f) {
    glColor3f(0.0f, 1.0f, 1.0f);
    glPushMatrix();
    glTranslatef(translateX, 0.0f, 0.0f);
    glBegin(GL_TRIANGLES);
    for (size_t i = 0; i < faces.size(); i += 3) {
        int idx1 = faces[i] * 3;
        int idx2 = faces[i + 1] * 3;
        int idx3 = faces[i + 2] * 3;
        if (vertices[idx1 + 2] < cutPlaneZ && vertices[idx2 + 2] < cutPlaneZ && vertices[idx3 + 2] < cutPlaneZ) continue;
        glVertex3f(vertices[idx1], vertices[idx1 + 1], vertices[idx1 + 2]);
        glVertex3f(vertices[idx2], vertices[idx2 + 1], vertices[idx2 + 2]);
        glVertex3f(vertices[idx3], vertices[idx3 + 1], vertices[idx3 + 2]);
    }
    glEnd();
    glPopMatrix();
}

void DrawMesh(const Mesh& m, float translateX = 0.0f){
    glColor3f(0.0f, 1.0f, 1.0f);
    glPushMatrix();
    glTranslatef(translateX, 0.0f, 0.0f);
    
    for (const auto& f: m.faces) {
        glBegin(GL_POLYGON);
        for (const auto& idx: f.vertexIndices){
            glVertex3f(m.vertices[idx].x, m.vertices[idx].y, m.vertices[idx].z);
        }
        glEnd();
    }
    
    glPopMatrix();
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    glTranslatef(panX, panY, -5.0f);
    glScalef(zoom, zoom, zoom);
    glRotatef(angleX, 1, 0, 0);
    glRotatef(angleY, 0, 1, 0);
    drawCuttingPlane();
    float color1[3] = { 0.0f, 1.0f, 1.0f };
    // drawMesh(vertices1, faces1, color1, loadSecondFile ? -1.5f : 0.0f);
    DrawMesh(humanoid);
    if (loadSecondFile) {
        float color2[3] = { 1.0f, 0.0f, 1.0f };
        // drawMesh(vertices2, faces2, color2, 1.5f);
    }
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

    GLfloat ambient[] = { 1.0f, 1.0f, 1.0f, 1.0f };  // Bright white ambient
    glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);  // Ambient only (no diffuse/specular)

    // Enable color material mode
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

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
    
    // std::string file1 = getFileFromUser("Enter path to first OBJ file: ");
    std::string file1 = "../../assets/obj/humanoid_robot.obj";
    // loadOBJ(file1.c_str(), vertices1, faces1);
    humanoid.loadOBJ(file1);
    humanoid.printMeshStats();
    // computeZBounds(vertices1);
    // if (fileCount == 2) {
    //     loadSecondFile = true;
    //     std::string file2 = getFileFromUser("Enter path to second OBJ file: ");
    //     loadOBJ(file2.c_str(), vertices2, faces2);
    // }
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