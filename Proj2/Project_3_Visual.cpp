#include "mesh_lib.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <algorithm>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "implot.h"

// View modes
enum ViewMode { MODE_FACE = 0, MODE_WIREFRAME, MODE_MATERIAL };
static ViewMode currentMode = MODE_FACE;

// Camera & interaction
static float camYaw = 0, camPitch = 0, camDistance = 5, panX = 0, panY = 0;
static bool leftMouse = false, rightMouse = false;
static double lastX = 0, lastY = 0;

// Cutting plane
static float cutPlaneZ = 0.0f;

// Thickness profiles & ranges
static std::vector<float> carbon, glue, steel;
static float thickMin = 1e9f, thickMax = -1e9f;
static std::vector<float> heights, steelR, glueR, carbonR;
static const float robotLength = 2.5f;

// Load two-column CSV (position,value)
std::vector<float> loadProfile(const std::string& fname) {
    std::vector<float> vals;
    std::ifstream in(fname);
    if (!in) { std::cerr << "Failed to open " << fname << "\n"; return vals; }
    std::string line; std::getline(in, line);
    while (std::getline(in, line)) {
        std::istringstream ss(line);
        float x, v; char comma;
        ss >> x >> comma >> v;
        vals.push_back(v);
        thickMin = glm::min(thickMin, v);
        thickMax = glm::max(thickMax, v);
    }
    return vals;
}

// Linear interpolation lookup
float lookup(const std::vector<float>& prof, float s) {
    if (prof.empty()) return 0.f;
    float idx = s * (prof.size() - 1);
    int i0 = (int)floor(idx), i1 = glm::min(i0 + 1, (int)prof.size() - 1);
    float f = idx - i0;
    return prof[i0] * (1 - f) + prof[i1] * f;
}

// GLFW callbacks
void key_callback(GLFWwindow*, int key, int, int action, int) {
    if (action == GLFW_PRESS && key == GLFW_KEY_W) {
        currentMode = ViewMode((currentMode + 1) % 3);
        GLenum m = (currentMode == MODE_WIREFRAME ? GL_LINE : GL_FILL);
        glPolygonMode(GL_FRONT_AND_BACK, m);
    }
}
void mouse_button_callback(GLFWwindow*, int b, int a, int) {
    if (b == GLFW_MOUSE_BUTTON_LEFT)  leftMouse = (a == GLFW_PRESS);
    if (b == GLFW_MOUSE_BUTTON_RIGHT) rightMouse = (a == GLFW_PRESS);
}
void cursor_position_callback(GLFWwindow*, double x, double y) {
    double dx = x - lastX, dy = y - lastY; lastX = x; lastY = y;
    if (leftMouse) { camYaw += dx * 0.3f; camPitch += dy * 0.3f; camPitch = glm::clamp(camPitch, -89.9f, 89.9f); }
    if (rightMouse) { float ps = 0.002f * camDistance; panX -= dx * ps; panY += dy * ps; }
}
void scroll_callback(GLFWwindow*, double, double yoff) {
    camDistance *= (1.0f - (float)yoff * 0.1f);
    camDistance = glm::max(camDistance, 0.1f);
}

// Shaders with clipping
const char* vSrc = R"(
#version 330 core
layout(location=0) in vec3 aPos;
layout(location=1) in vec3 aColor;
uniform mat4 uMVP;
uniform vec4 uClipPlane;
out vec3 fragColor;
void main(){
    gl_Position = uMVP * vec4(aPos,1.0);
    fragColor = aColor;
    gl_ClipDistance[0] = dot(vec4(aPos,1.0), uClipPlane);
}
)";
const char* fSrc = R"(
#version 330 core
in vec3 fragColor;
out vec4 outColor;
void main(){ outColor = vec4(fragColor,1.0); }
)";

GLuint compileShader(GLenum t, const char* src) {
    GLuint s = glCreateShader(t);
    glShaderSource(s, 1, &src, nullptr);
    glCompileShader(s);
    GLint ok; glGetShaderiv(s, GL_COMPILE_STATUS, &ok);
    if (!ok) { char b[512]; glGetShaderInfoLog(s, 512, nullptr, b); std::cerr << b; }
    return s;
}
GLuint createProgram() {
    GLuint vs = compileShader(GL_VERTEX_SHADER, vSrc), fs = compileShader(GL_FRAGMENT_SHADER, fSrc), p = glCreateProgram();
    glAttachShader(p, vs); glAttachShader(p, fs); glLinkProgram(p);
    GLint ok; glGetProgramiv(p, GL_LINK_STATUS, &ok);
    if (!ok) { char b[512]; glGetProgramInfoLog(p, 512, nullptr, b); std::cerr << b; }
    glDeleteShader(vs); glDeleteShader(fs); return p;
}

int main() {
    if (!glfwInit()) return 1;
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    GLFWwindow* w = glfwCreateWindow(1280, 720, "Mesh Viewer", nullptr, nullptr);
    glfwMakeContextCurrent(w);
    gladLoadGLLoader((GLADloadproc)glfwGetProcAddress);
    glfwSetKeyCallback(w, key_callback);
    glfwSetMouseButtonCallback(w, mouse_button_callback);
    glfwSetCursorPosCallback(w, cursor_position_callback);
    glfwSetScrollCallback(w, scroll_callback);
    glfwGetCursorPos(w, &lastX, &lastY);

    IMGUI_CHECKVERSION(); ImGui::CreateContext(); ImPlot::CreateContext();
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(w, true);
    ImGui_ImplOpenGL3_Init("#version 330 core");

    // Enable alpha blending for translucent plots
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // Load CSV profiles from assets/csv
    carbon = loadProfile("assets/csv/carbon_thickness.csv");
    glue = loadProfile("assets/csv/glue_thickness.csv");
    steel = loadProfile("assets/csv/steel_thickness.csv");
    int N = (int)carbon.size();
    heights.resize(N); steelR.resize(N); glueR.resize(N); carbonR.resize(N);
    for (int i = 0;i < N;++i) {
        float s = i / (float)(N - 1);
        heights[i] = s * robotLength;
        steelR[i] = steel[i] / 100.0f;
        glueR[i] = steelR[i] + glue[i] / 100.0f;
        carbonR[i] = glueR[i] + carbon[i] / 100.0f;
    }

    Mesh mesh; if (!mesh.loadOBJ("assets/obj/humanoid_robot_2d.obj")) return 1;
    float minZ = 1e9f, maxZ = -1e9f;
    for (auto& v : mesh.vertices) { minZ = glm::min(minZ, v.z); maxZ = glm::max(maxZ, v.z); }
    cutPlaneZ = 0.0f;

    auto& V = mesh.vertices;
    std::vector<unsigned int> I; I.reserve(mesh.faces.size() * 3);
    for (auto& F : mesh.faces) { I.push_back(F.v[0]);I.push_back(F.v[1]);I.push_back(F.v[2]); }
    GLuint VAO, VBO, CBO, EBO;
    glGenVertexArrays(1, &VAO); glGenBuffers(1, &VBO); glGenBuffers(1, &CBO); glGenBuffers(1, &EBO);
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, V.size() * sizeof(glm::vec3), V.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0); glEnableVertexAttribArray(0);
    std::vector<glm::vec3> cols(V.size(), glm::vec3(0.8f));
    glBindBuffer(GL_ARRAY_BUFFER, CBO);
    glBufferData(GL_ARRAY_BUFFER, cols.size() * sizeof(glm::vec3), cols.data(), GL_DYNAMIC_DRAW);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0); glEnableVertexAttribArray(1);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, I.size() * sizeof(unsigned int), I.data(), GL_STATIC_DRAW);
    glBindVertexArray(0);

    GLuint prog = createProgram();
    GLint locMVP = glGetUniformLocation(prog, "uMVP"), locClip = glGetUniformLocation(prog, "uClipPlane");
    glEnable(GL_DEPTH_TEST); glEnable(GL_CLIP_DISTANCE0);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    while (!glfwWindowShouldClose(w)) {
        glfwPollEvents(); int W, H; glfwGetFramebufferSize(w, &W, &H); H = glm::max(H, 1);

        ImGui_ImplOpenGL3_NewFrame(); ImGui_ImplGlfw_NewFrame(); ImGui::NewFrame();
        // Set wider menu
        ImGui::SetNextWindowSize(ImVec2(260, 0), ImGuiCond_Once);
        ImGui::Begin("View Mode", nullptr, ImGuiWindowFlags_NoCollapse);
        ImGui::Text("(Press W to toggle)");
        if (ImGui::RadioButton("Face", currentMode == MODE_FACE)) { currentMode = MODE_FACE; glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); }
        if (ImGui::RadioButton("Wireframe", currentMode == MODE_WIREFRAME)) { currentMode = MODE_WIREFRAME; glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); }
        if (ImGui::RadioButton("Material", currentMode == MODE_MATERIAL)) { currentMode = MODE_MATERIAL; glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); }
        ImGui::SliderFloat("Cut Z", &cutPlaneZ, minZ, maxZ);
        ImGui::End();

        if (currentMode != MODE_MATERIAL) {
            // reset to default gray
            std::vector<glm::vec3> def(V.size(), glm::vec3(0.8f));
            glBindBuffer(GL_ARRAY_BUFFER, CBO);
            glBufferSubData(GL_ARRAY_BUFFER, 0, def.size() * sizeof(glm::vec3), def.data());
        }

        if (currentMode == MODE_MATERIAL) {
            ImGui::Begin("Thickness Profile");
            if (ImPlot::BeginPlot("##thickness", ImVec2(-1, 200))) {
                ImPlot::SetupAxisLimits(ImAxis_X1, 0.0, carbonR.back(), ImPlotCond_Always);
                ImPlot::SetupAxisLimits(ImAxis_Y1, robotLength, 0.0, ImPlotCond_Always);
                ImPlot::SetupAxis(ImAxis_X1, "Radius (m)");
                ImPlot::SetupAxis(ImAxis_Y1, "Height (m)", ImPlotAxisFlags_NoDecorations);
                ImPlot::PushStyleColor(ImPlotCol_Fill, IM_COL32(50, 100, 200, 100));
                ImPlot::PlotShaded("Steel", steelR.data(), heights.data(), N);
                ImPlot::PopStyleColor();
                ImPlot::PushStyleColor(ImPlotCol_Fill, IM_COL32(200, 130, 50, 100));
                ImPlot::PlotShaded("Glue", glueR.data(), heights.data(), N);
                ImPlot::PopStyleColor();
                ImPlot::PushStyleColor(ImPlotCol_Fill, IM_COL32(50, 200, 50, 100));
                ImPlot::PlotShaded("Carbon", carbonR.data(), heights.data(), N);
                ImPlot::PopStyleColor();
                ImPlot::EndPlot();
            }
            ImGui::End();
        }

        ImGui::Render();

        glm::vec3 tgt(panX, panY, 0);
        float yR = glm::radians(camYaw), pR = glm::radians(camPitch);
        glm::vec3 pos = tgt + glm::vec3(camDistance * cos(pR) * sin(yR), camDistance * sin(pR), camDistance * cos(pR) * cos(yR));
        glm::mat4 view = glm::lookAt(pos, tgt, glm::vec3(0, 1, 0));
        glm::mat4 proj = glm::perspective(glm::radians(60.0f), W / (float)H, 0.1f, 1000.0f);
        glm::mat4 MVP = proj * view;

        // Material coloring update
        if (currentMode == MODE_MATERIAL) {
            std::vector<glm::vec3> matCols(V.size());
            float minX = 1e9f, maxX = -1e9f;
            for (auto& v : V) { minX = glm::min(minX, v.x); maxX = glm::max(maxX, v.x); }
            for (auto& F : mesh.faces) {
                float s = ((V[F.v[0]].x + V[F.v[1]].x + V[F.v[2]].x) / 3.0f - minX) / (maxX - minX);
                float tv = 0;
                switch (F.mat) { case 0:tv = lookup(steel, s);break; case 1:tv = lookup(glue, s);break; case 2:tv = lookup(carbon, s);break; }
                                       float tn = glm::clamp((tv - thickMin) / (thickMax - thickMin), 0.0f, 1.0f);
                                       glm::vec3 col = (tn < 0.5f ? glm::mix(glm::vec3(0, 0, 1), glm::vec3(0, 1, 0), tn * 2)
                                           : glm::mix(glm::vec3(0, 1, 0), glm::vec3(1, 0, 0), (tn - 0.5f) * 2));
                                       matCols[F.v[0]] = col; matCols[F.v[1]] = col; matCols[F.v[2]] = col;
            }
            glBindBuffer(GL_ARRAY_BUFFER, CBO);
            glBufferSubData(GL_ARRAY_BUFFER, 0, matCols.size() * sizeof(glm::vec3), matCols.data());
        }

        glViewport(0, 0, W, H);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glUseProgram(prog);
        glUniformMatrix4fv(locMVP, 1, GL_FALSE, glm::value_ptr(MVP));
        glUniform4f(locClip, 0, 0, 1, -cutPlaneZ);
        glBindVertexArray(VAO);
        glDrawElements(GL_TRIANGLES, (GLsizei)I.size(), GL_UNSIGNED_INT, 0);
        glBindVertexArray(0);

        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(w);
    }
    return 0;
}
