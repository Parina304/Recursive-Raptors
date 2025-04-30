#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include "mesh_lib.h"
#include "mesh_merge.h"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include "cppcolormap.h"
#include "ColormapWrapper.h"
#include <xtensor/xio.hpp>

// Global variables
float angleX = 0.0f, angleY = 0.0f, zoom = 1.0f, panX = 0.0f, panY = 0.0f;
float cutPlaneZ = 0.0f;
int window_width = 1600, window_height = 900;
bool flag_left_mouse_pressed = false;
bool flag_right_mouse_pressed = false;
int last_mouse_x = 0, last_mouse_y = 0;
bool flag_color_calced = false;

Mesh humanoid;

// Function prototypes
void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods);
void cursor_position_callback(GLFWwindow* window, double xpos, double ypos);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void DrawMesh(const Mesh& m, float translateX = 0.0f);
void CalcMeshColors(Mesh& m);
void initOpenGL();

void DrawMesh(const Mesh& m, float translateX = 0.0f) {
    glColor3f(0.0f, 1.0f, 1.0f);
    glPushMatrix();
    glTranslatef(translateX, 0.0f, 0.0f);

    for (const auto& f : m.faces) {
        glColor3f(f.r, f.g, f.b);
        glBegin(GL_POLYGON);
        for (const auto& idx : f.vertexIndices) {
            glNormal3f(m.vertices[idx].x, m.vertices[idx].y, m.vertices[idx].z);
            glVertex3f(m.vertices[idx].x, m.vertices[idx].y, m.vertices[idx].z);
        }
        glEnd();
    }
    glPopMatrix();
}

void CalcMeshColors(Mesh& m) {
    float y_mean = 0;
    xt::xtensor<double, 1> color_idx;

    for (auto& f : m.faces) {
        y_mean = 0;
        for (const auto& idx : f.vertexIndices) {
            y_mean += m.vertices[idx].y;
        }
        y_mean /= f.vertexIndices.size();

        color_idx = { y_mean };
        auto c = cppcolormap::as_colors(color_idx, cppcolormap::viridis(), m.y_min, m.y_max);
        f.r = static_cast<float>(c(0, 0));
        f.g = static_cast<float>(c(0, 1));
        f.b = static_cast<float>(c(0, 2));
    }
}

void initOpenGL() {
    glEnable(GL_DEPTH_TEST);

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);

    // Set light properties
    GLfloat ambient[] = { 0.5f, 0.5f, 0.5f, 1.0f };
    GLfloat diffuse[] = { 0.05f, 0.05f, 0.05f, 1.0f };
    GLfloat specular[] = { 0.0f, 0.0f, 0.0f, 1.0f };
    GLfloat position[] = { 1.0f, 1.0f, 1.0f, 0.0f };

    glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
    glLightfv(GL_LIGHT0, GL_POSITION, position);

    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

    GLfloat mat_specular[] = { 0.0f, 0.0f, 0.0f, 1.0f };
    GLfloat mat_shininess[] = { 0.0f };
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0, (double)window_width / (double)window_height, 1.0, 100.0);
    glMatrixMode(GL_MODELVIEW);
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    if (height == 0) height = 1;
    glViewport(0, 0, width, height);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0, (double)width / (double)height, 1.0, 100.0);
    glMatrixMode(GL_MODELVIEW);
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    if (action == GLFW_PRESS || action == GLFW_REPEAT) {
        switch (key) {
        case GLFW_KEY_W: angleX -= 5.0f; break;
        case GLFW_KEY_S: angleX += 5.0f; break;
        case GLFW_KEY_A: angleY -= 5.0f; break;
        case GLFW_KEY_D: angleY += 5.0f; break;
        case GLFW_KEY_EQUAL: zoom += 0.1f; break;
        case GLFW_KEY_MINUS: zoom -= 0.1f; break;
        case GLFW_KEY_J: panX -= 0.1f; break;
        case GLFW_KEY_L: panX += 0.1f; break;
        case GLFW_KEY_I: panY += 0.1f; break;
        case GLFW_KEY_K: panY -= 0.1f; break;
        case GLFW_KEY_ESCAPE: glfwSetWindowShouldClose(window, true); break;
        }
    }
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        flag_left_mouse_pressed = (action == GLFW_PRESS);
    } else if (button == GLFW_MOUSE_BUTTON_RIGHT) {
        flag_right_mouse_pressed = (action == GLFW_PRESS);
    }
}

void cursor_position_callback(GLFWwindow* window, double xpos, double ypos) {
    int dx = xpos - last_mouse_x;
    int dy = ypos - last_mouse_y;

    if (flag_left_mouse_pressed) {
        angleX += dy * 0.5f;
        angleY += dx * 0.5f;
    } else if (flag_right_mouse_pressed) {
        panX += dx * 0.01f;
        panY -= dy * 0.01f;
    }

    last_mouse_x = xpos;
    last_mouse_y = ypos;
}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
    zoom += yoffset * 0.1f;
}

int main() {
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        return -1;
    }

    GLFWwindow* window = glfwCreateWindow(window_width, window_height, "OBJ Viewer", nullptr, nullptr);
    if (!window) {
        std::cerr << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);
    glewInit();

    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetKeyCallback(window, key_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetCursorPosCallback(window, cursor_position_callback);
    glfwSetScrollCallback(window, scroll_callback);

    initOpenGL();

    humanoid.loadOBJ("../../assets/obj/humanoid_robot.obj");
    humanoid.printMeshStats();
    if (!flag_color_calced) {
        CalcMeshColors(humanoid);
        flag_color_calced = true;
    }

    while (!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glLoadIdentity();
        glTranslatef(panX, panY, -5.0f);
        glScalef(zoom, zoom, zoom);
        glRotatef(angleX, 1, 0, 0);
        glRotatef(angleY, 0, 1, 0);

        DrawMesh(humanoid);

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}