#ifdef _WIN32
#include <glut.h>
#else
#include <GL/glut.h>
#endif

#include "mesh_lib.h"
#include "mesh_merge.h"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>

// Globals for OBJ files
std::vector<float> vertices1;
std::vector<int> faces1;
std::vector<float> vertices2;
std::vector<int> faces2;

bool loadSecondFile = false;  // Flag to check if user wants to load 2 files

float angleX = 0.0f, angleY = 0.0f, zoom = 1.0f, panX = 0.0f, panY = 0.0f;
float cutPlaneZ = 0.0f;
float minZ = 0.0f, maxZ = 0.0f;

std::string getFileFromUser(const std::string& promptMessage) {
    std::string filePath;
    std::cout << promptMessage;
    std::getline(std::cin, filePath);
    return filePath;
}

int getNumberOfFiles() {
    std::string input;
    int count = 0;
    while (true) {
        std::cout << "How many OBJ files do you want to load (1 or 2)? ";
        std::getline(std::cin, input);  // Read as string to handle any input

        if (input == "1" || input == "2") {
            count = std::stoi(input);  // Safe to convert now
            break;  // Valid input, break loop
        }
        else {
            std::cout << "Invalid input. Please enter only 1 or 2.\n";
        }
    }
    return count;
}


void loadOBJ(const char* filename, std::vector<float>& vertices, std::vector<int>& faces) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Cannot open file: " << filename << std::endl;
        return;
    }
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        if (line.substr(0, 2) == "v ") {
            float x, y, z;
            iss.ignore(2);
            iss >> x >> y >> z;
            vertices.push_back(x);
            vertices.push_back(y);
            vertices.push_back(z);
        }
        else if (line.substr(0, 2) == "f ") {
            int a, b, c;
            iss.ignore(2);
            iss >> a >> b >> c;
            faces.push_back(a - 1);
            faces.push_back(b - 1);
            faces.push_back(c - 1);
        }
    }
    file.close();
}

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

void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    glTranslatef(panX, panY, -5.0f);
    glScalef(zoom, zoom, zoom);
    glRotatef(angleX, 1, 0, 0);
    glRotatef(angleY, 0, 1, 0);
    drawCuttingPlane();
    float color1[3] = { 0.0f, 1.0f, 1.0f };
    drawMesh(vertices1, faces1, color1, loadSecondFile ? -1.5f : 0.0f);
    if (loadSecondFile) {
        float color2[3] = { 1.0f, 0.0f, 1.0f };
        drawMesh(vertices2, faces2, color2, 1.5f);
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
    case 27: exit(0);
    }
    glutPostRedisplay();
}

void initOpenGL() {
    glEnable(GL_DEPTH_TEST);

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);

    GLfloat ambient[] = { 1.0f, 1.0f, 1.0f, 1.0f };  // Bright white ambient
    glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);  // Ambient only (no diffuse/specular)

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0, 1.0, 1.0, 100.0);
    glMatrixMode(GL_MODELVIEW);
}


/*
void initOpenGL() {
    glEnable(GL_DEPTH_TEST);  // Keep depth for 3D

    // Optional: If you want to keep wireframe view
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);  // Wireframe for better shape visualization

    // Optional: If you want to hide inner triangles
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0, 1.0, 1.0, 100.0);
    glMatrixMode(GL_MODELVIEW);

    // Optional: Simple ambient light to keep everything bright

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    GLfloat ambient[] = { 1.0f, 1.0f, 1.0f, 1.0f };
    glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);

}
*/


// int main(int argc, char** argv) {
//     int fake_argc = 1;
//     char* fake_argv[] = { argv[0], nullptr };
//     int fileCount = getNumberOfFiles();
//     std::string file1 = getFileFromUser("Enter path to first OBJ file: ");
//     loadOBJ(file1.c_str(), vertices1, faces1);
//     computeZBounds(vertices1);
//     if (fileCount == 2) {
//         loadSecondFile = true;
//         std::string file2 = getFileFromUser("Enter path to second OBJ file: ");
//         loadOBJ(file2.c_str(), vertices2, faces2);
//     }
//     glutInit(&fake_argc, fake_argv);
//     glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
//     glutInitWindowSize(800, 600);
//     glutCreateWindow("OBJ Viewer");
//     initOpenGL();
//     glutDisplayFunc(display);
//     glutKeyboardFunc(keyboard);
//     glutMainLoop();
//     return 0;
// }
// ...existing code...

// Globals for mouse interaction
bool isLeftMousePressed = false;
bool isRightMousePressed = false;
int lastMouseX = 0, lastMouseY = 0;

void mouse(int button, int state, int x, int y) {
    if (button == GLUT_LEFT_BUTTON) {
        isLeftMousePressed = (state == GLUT_DOWN);
    } else if (button == GLUT_RIGHT_BUTTON) {
        isRightMousePressed = (state == GLUT_DOWN);
    }
    lastMouseX = x;
    lastMouseY = y;
}

void motion(int x, int y) {
    int dx = x - lastMouseX;
    int dy = y - lastMouseY;

    if (isLeftMousePressed) {
        // Rotate the model
        angleX += dy * 0.5f;
        angleY += dx * 0.5f;
    } else if (isRightMousePressed) {
        // Pan the model
        panX += dx * 0.01f;
        panY -= dy * 0.01f;
    }

    lastMouseX = x;
    lastMouseY = y;

    glutPostRedisplay();
}

void mouseWheel(int button, int dir, int x, int y) {
    if (dir > 0) {
        zoom += 0.1f;  // Zoom in
    } else {
        zoom -= 0.1f;  // Zoom out
    }
    glutPostRedisplay();
}

int main(int argc, char** argv) {
    int fake_argc = 1;
    char* fake_argv[] = { argv[0], nullptr };
    int fileCount = getNumberOfFiles();
    std::string file1 = getFileFromUser("Enter path to first OBJ file: ");
    loadOBJ(file1.c_str(), vertices1, faces1);
    computeZBounds(vertices1);
    if (fileCount == 2) {
        loadSecondFile = true;
        std::string file2 = getFileFromUser("Enter path to second OBJ file: ");
        loadOBJ(file2.c_str(), vertices2, faces2);
    }
    glutInit(&fake_argc, fake_argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(800, 600);
    glutCreateWindow("OBJ Viewer");
    initOpenGL();
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    // glutMouseWheelFunc(mouseWheel);  // Requires freeglut
    glutMainLoop();
    return 0;
}