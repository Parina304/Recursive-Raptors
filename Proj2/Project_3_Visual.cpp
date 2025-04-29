// mesh_viewer.cpp – single‑file demo that toggles between a 3‑D and 2‑D OBJ
// model and still supports Face, Wireframe, and Material colouring modes.
//
// Build:
//   g++ mesh_viewer.cpp -I<path‑to‑glad> -I<path‑to‑glm> -lglfw -ldl -lGL -std=c++17 -o viewer
// (plus the ImGui / ImPlot / glad sources in your build system)
//
// Controls:
//   • Radio buttons   : pick 3‑D vs 2‑D mesh, view mode, and cut plane.
//   • W key           : cycles Face/Wireframe/Material if the GUI is hidden.
//   • Mouse (LMB)     : orbit
//   • Mouse (RMB)     : pan
//   • Scroll          : dolly zoom
//--------------------------------------------------------------------------

#ifdef _WIN32
#define BASE_PATH "./"
#else
#define BASE_PATH "../"
#endif

#include "mesh_lib.h"          // simple OBJ + Face data structure (user‑supplied)

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
// #include <unistd.h>

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "implot.h"

//--------------------------------------------------------------------------
// ENUMERATIONS
//--------------------------------------------------------------------------

// 3 visual styles (shading) that apply no matter which mesh is loaded
enum ViewMode { MODE_FACE = 0, MODE_WIREFRAME, MODE_MATERIAL };
static ViewMode currentViewMode = MODE_FACE;

// Which mesh is currently resident on the GPU
enum MeshType { MESH_3D = 0, MESH_2D };
static MeshType currentMeshType = MESH_3D;

//--------------------------------------------------------------------------
// CAMERA / INTERACTION STATE
//--------------------------------------------------------------------------
static float camYaw = 0.0f, camPitch = 0.0f, camDistance = 5.0f;
static float panX = 0.0f, panY = 0.0f;
static bool  leftMouse = false, rightMouse = false;
static double lastX = 0.0, lastY = 0.0;

//--------------------------------------------------------------------------
// CUTTING PLANE
//--------------------------------------------------------------------------
static float cutPlaneZ = 0.0f;   // in model space (metres)

//--------------------------------------------------------------------------
// THICKNESS PROFILES (loaded from CSV)
//--------------------------------------------------------------------------
static std::vector<float> carbon, glue, steel;
static float thickMin = 1e9f, thickMax = -1e9f;
static std::vector<float> heights, thermal, steelR, glueR, carbonR, thermalR;
static const float robotLength = 2.5f;                   // m
static std::vector<float> thermZ, thermT;                // thermal profile

//--------------------------------------------------------------------------
// OPENGL RESOURCES THAT NEED TO BE RECREATED WHENEVER WE SWITCH MESH
//--------------------------------------------------------------------------
static GLuint VAO = 0, VBO = 0, CBO = 0, EBO = 0;
static Mesh   mesh;                       // CPU‑side representation (faces, verts)
static std::vector<glm::vec3>   V;
static std::vector<Face>        F;
static std::vector<unsigned int>I;

//--------------------------------------------------------------------------
// FORWARD DECLARATIONS
//--------------------------------------------------------------------------
bool  loadThicknessCSV();                     // called once at start‑up
bool  loadMeshToGPU(const std::string& path); // used every time mesh changes
GLuint compileShader(GLenum, const char*);
GLuint createProgram();

//--------------------------------------------------------------------------
// CSV HELPERS (single‑column and two‑column)
//--------------------------------------------------------------------------

// makes excutable compatible in ./ or ./build  
std::string PrependBasePath (const std::string& path){
    return BASE_PATH + path;
}

// Load two-column CSV (position,value)
std::vector<float> loadProfile1(const std::string& fname) {
    std::vector<float> vals;
    std::ifstream in(fname);
    if (!in) {
        std::cerr << "[CSV] Cannot open " << fname << "\n";
        return vals;
    }
    std::string line; std::getline(in, line);              // skip header
    while (std::getline(in, line)) {
        std::istringstream ss(line);
        float s, v; char comma;
        ss >> s >> comma >> v;             // s is the 0…1 normalised height
        vals.push_back(v);
        thickMin = std::min(thickMin, v);
        thickMax = std::max(thickMax, v);
    }
    return vals;
}

static void loadThermal(const std::string& fname)
{
    std::ifstream in(fname);
    if (!in) { std::cerr << "[CSV] Cannot open " << fname << "\n"; return; }
    std::string line; std::getline(in, line);
    while (std::getline(in, line)) {
        std::istringstream ss(line);
        float z, t; char comma;
        ss >> z >> comma >> t;
        thermZ.push_back(z);
        thermT.push_back(t);
    }
}

// Linear look‑ups for material and thermal profiles
static float interp1D(const std::vector<float>& arr, float s)
{
    if (arr.empty()) return 0.0f;
    float idx = s * (arr.size() - 1);
    int   i0 = (int)std::floor(idx);
    int   i1 = std::min(i0 + 1, (int)arr.size() - 1);
    float f = idx - i0;
    return arr[i0] * (1.0f - f) + arr[i1] * f;
}

static float interpThermal(float h)
{
    if (thermZ.empty()) return 0.0f;
    if (h <= thermZ.front()) return thermT.front();
    if (h >= thermZ.back())  return thermT.back();
    int lo = 0, hi = (int)thermZ.size() - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (thermZ[mid] <= h) lo = mid; else hi = mid;
    }
    float f = (h - thermZ[lo]) / (thermZ[hi] - thermZ[lo]);
    return thermT[lo] * (1.0f - f) + thermT[hi] * f;
}

//--------------------------------------------------------------------------
// MESH LOADING + GPU BUFFER UPLOAD
//--------------------------------------------------------------------------

bool loadMeshToGPU(const std::string& path)
{
    // Load OBJ into CPU structures
    if (!mesh.loadOBJ(path)) {
        std::cerr << "[OBJ] Failed to load " << path << "\n";
        return false;
    }
    V = mesh.vertices;
    F = mesh.faces;

    // Build (or rebuild) index buffer
    I.clear(); I.reserve(F.size() * 3);
    for (const auto& f : F) {
        I.push_back(f.v[0]); I.push_back(f.v[1]); I.push_back(f.v[2]);
    }

    // Lazy create GL objects once
    if (VAO == 0) {
        glGenVertexArrays(1, &VAO);
        glGenBuffers(1, &VBO);
        glGenBuffers(1, &CBO);
        glGenBuffers(1, &EBO);
    }

    // ── Upload vertex positions ────────────────────────────────────────────────
    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, V.size() * sizeof(glm::vec3), V.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
    glEnableVertexAttribArray(0);

    // ── Default (grey) colours – updated every frame when needed ───────────────
    std::vector<glm::vec3> defCols(V.size(), glm::vec3(0.8f));
    glBindBuffer(GL_ARRAY_BUFFER, CBO);
    glBufferData(GL_ARRAY_BUFFER, defCols.size() * sizeof(glm::vec3), defCols.data(), GL_DYNAMIC_DRAW);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
    glEnableVertexAttribArray(1);

    // ── Indices ────────────────────────────────────────────────────────────────
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, I.size() * sizeof(unsigned int), I.data(), GL_STATIC_DRAW);

    glBindVertexArray(0);
    return true;
}

//--------------------------------------------------------------------------
// SHADERS (single clip plane)
//--------------------------------------------------------------------------
static const char* vSrc = R"(
#version 330 core
layout(location=0) in vec3 aPos;
layout(location=1) in vec3 aColor;

uniform mat4 uMVP;
uniform vec4 uClipPlane;

out vec3 vColor;

void main() {
    gl_Position = uMVP * vec4(aPos, 1.0);
    vColor      = aColor;
    gl_ClipDistance[0] = dot(vec4(aPos,1.0), uClipPlane);
})";

static const char* fSrc = R"(
#version 330 core
in  vec3 vColor;
out vec4 FragColor;
void main() { FragColor = vec4(vColor, 1.0); }
)";

GLuint compileShader(GLenum type, const char* src)
{
    GLuint s = glCreateShader(type);
    glShaderSource(s, 1, &src, nullptr);
    glCompileShader(s);
    GLint ok; glGetShaderiv(s, GL_COMPILE_STATUS, &ok);
    if (!ok) {
        char log[1024];
        glGetShaderInfoLog(s, 1024, nullptr, log);
        std::cerr << "[GLSL] " << log << "\n";
    }
    return s;
}

GLuint createProgram()
{
    GLuint vs = compileShader(GL_VERTEX_SHADER, vSrc);
    GLuint fs = compileShader(GL_FRAGMENT_SHADER, fSrc);
    GLuint p = glCreateProgram();
    glAttachShader(p, vs);
    glAttachShader(p, fs);
    glLinkProgram(p);
    GLint ok; glGetProgramiv(p, GL_LINK_STATUS, &ok);
    if (!ok) {
        char log[1024];
        glGetProgramInfoLog(p, 1024, nullptr, log);
        std::cerr << "[GLSL] " << log << "\n";
    }
    glDeleteShader(vs);
    glDeleteShader(fs);
    return p;
}

//--------------------------------------------------------------------------
// GLFW CALLBACKS
//--------------------------------------------------------------------------
static void key_cb(GLFWwindow* window, int key, int, int action, int)
{
    if (action == GLFW_PRESS) {
        if(key == GLFW_KEY_W){
            currentViewMode = ViewMode((currentViewMode + 1) % 3);
            GLenum mode = (currentViewMode == MODE_WIREFRAME ? GL_LINE : GL_FILL);
            glPolygonMode(GL_FRONT_AND_BACK, mode);
        } else if (key == GLFW_KEY_ESCAPE || key == GLFW_KEY_Q){
            glfwSetWindowShouldClose(window, GLFW_TRUE); // Signal to close the window
        }
    }
}

static void mouse_btn_cb(GLFWwindow*, int button, int action, int)
{
    if (button == GLFW_MOUSE_BUTTON_LEFT)  leftMouse = (action == GLFW_PRESS);
    if (button == GLFW_MOUSE_BUTTON_RIGHT) rightMouse = (action == GLFW_PRESS);
}

static void cursor_cb(GLFWwindow*, double x, double y)
{
    double dx = x - lastX, dy = y - lastY;
    lastX = x; lastY = y;

    if (leftMouse) {
        camYaw += float(dx) * 0.3f;
        camPitch += float(dy) * 0.3f;
        camPitch = std::clamp(camPitch, -89.9f, 89.9f);
    }
    if (rightMouse) {
        float s = 0.002f * camDistance;
        panX -= float(dx) * s;
        panY += float(dy) * s;
    }
}

static void scroll_cb(GLFWwindow*, double, double yoff)
{
    camDistance *= (1.0f - float(yoff) * 0.1f);
    camDistance = std::max(camDistance, 0.1f);
}

//--------------------------------------------------------------------------
// MAIN
//--------------------------------------------------------------------------
int main()
{
    if (!glfwInit()) return 1;
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    GLFWwindow* win = glfwCreateWindow(1280, 720, "Mesh Viewer 2D/3D", nullptr, nullptr);
    glfwMakeContextCurrent(win);

    gladLoadGLLoader((GLADloadproc)glfwGetProcAddress);

    glfwSetKeyCallback(win, key_cb);
    glfwSetMouseButtonCallback(win, mouse_btn_cb);
    glfwSetCursorPosCallback(win, cursor_cb);
    glfwSetScrollCallback(win, scroll_cb);
    glfwGetCursorPos(win, &lastX, &lastY);

    // ── ImGui + ImPlot ─────────────────────────────────────────────────────────
    IMGUI_CHECKVERSION(); ImGui::CreateContext(); ImPlot::CreateContext();
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(win, true);
    ImGui_ImplOpenGL3_Init("#version 330 core");

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    //------------------------------------------------------------------ CSV ----
    carbon = loadProfile1(PrependBasePath("assets/csv/carbon_thickness.csv"));
    glue = loadProfile1(PrependBasePath("assets/csv/glue_thickness.csv"));
    steel = loadProfile1(PrependBasePath("assets/csv/steel_thickness.csv"));
    loadThermal(PrependBasePath("assets/csv/Thickness_1_hr.csv"));

    // Load CSV profiles from assets/csv
    // carbon = loadProfile(PrependBasePath("assets/csv/carbon_thickness.csv"));
    // glue = loadProfile(PrependBasePath("assets/csv/glue_thickness.csv"));
    // steel = loadProfile(PrependBasePath("assets/csv/steel_thickness.csv"));
    int N = (int)carbon.size();
    heights.resize(N);
    thermal.resize(N);
    steelR.resize(N); glueR.resize(N); carbonR.resize(N); thermalR.resize(N);

    for (int i = 0; i < N; ++i) {
        float s = i / float(N - 1);
        heights[i] = s * robotLength;
        thermal[i] = interpThermal(heights[i]);
        steelR[i] = steel[i] / 100.0f;
        glueR[i] = steelR[i] + glue[i] / 100.0f;
        carbonR[i] = glueR[i] + carbon[i] / 100.0f;
        thermalR[i] = carbonR[i] + thermal[i] / 100.0f;
    }

    //----------------------------------------------------------------- OBJ ----
    const std::string OBJ_3D = PrependBasePath("assets/obj/humanoid_robot.obj");
    const std::string OBJ_2D = PrependBasePath("assets/obj/humanoid_robot_2d.obj");

    if (!loadMeshToGPU(OBJ_3D)) return 1;          // start with 3‑D model

    //----------------------------------------------------------------- GLSL ----

    GLuint prog = createProgram();
    GLint  locMVP = glGetUniformLocation(prog, "uMVP");
    GLint  locClip = glGetUniformLocation(prog, "uClipPlane");

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CLIP_DISTANCE0);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    //---------------------------------------------------------------- MAIN LOOP -
    while (!glfwWindowShouldClose(win)) {
        glfwPollEvents();
        int fbW, fbH; glfwGetFramebufferSize(win, &fbW, &fbH);
        fbH = std::max(fbH, 1);

        // ── ImGui frame ───────────────────────────────────────────────────────
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        ImGui::SetNextWindowSize(ImVec2(280, 0), ImGuiCond_Once);
        ImGui::Begin("Controls");

        // --- Mesh selection ---------------------------------------------------
        ImGui::Text("Mesh");
        if (ImGui::RadioButton("3D", currentMeshType == MESH_3D)) {
            if (currentMeshType != MESH_3D) {
                loadMeshToGPU(OBJ_3D);
                currentMeshType = MESH_3D;
            }
        }
        ImGui::SameLine();
        if (ImGui::RadioButton("2D", currentMeshType == MESH_2D)) {
            if (currentMeshType != MESH_2D) {
                loadMeshToGPU(OBJ_2D);
                currentMeshType = MESH_2D;
            }
        }

        ImGui::Separator();
        ImGui::Text("View Mode (or press W)");
        if (ImGui::RadioButton("Face", currentViewMode == MODE_FACE)) { currentViewMode = MODE_FACE;      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); }
        if (ImGui::RadioButton("Wireframe", currentViewMode == MODE_WIREFRAME)) { currentViewMode = MODE_WIREFRAME; glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); }
        if (ImGui::RadioButton("Material", currentViewMode == MODE_MATERIAL)) { currentViewMode = MODE_MATERIAL;  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); }

        ImGui::Separator();
        ImGui::SliderFloat("Cut Z", &cutPlaneZ, -1.0f, 1.0f, "%.2f m");
        ImGui::End();

        // ---------------------------------------------------------------------
        // Thickness stacked‑area plot if Material mode is active
        if (currentViewMode == MODE_MATERIAL) {
            ImGui::Begin("Thickness Profile");
            if (ImPlot::BeginPlot("##profile", ImVec2(-1, 300))) {
                ImPlot::SetupAxisLimits(ImAxis_Y1, 0.0f, *std::max_element(thermalR.begin(), thermalR.end()), ImPlotCond_Always);
                ImPlot::SetupAxisLimits(ImAxis_X1, robotLength, 0.0f, ImPlotCond_Always);
                ImPlot::SetupAxis(ImAxis_Y1, "Radius (m)");
                // ImPlot::SetupAxis(ImAxis_X1, "Height (m)", ImPlotAxisFlags_NoDecorations);
                ImPlot::SetupAxis(ImAxis_X1, "Height (m)");

                ImPlot::PushStyleColor(ImPlotCol_Fill, IM_COL32(200,50, 200,100));   ImPlot::PlotShaded("TPS", heights.data(), thermalR.data(), N);   ImPlot::PopStyleColor();
                ImPlot::PushStyleColor(ImPlotCol_Fill, IM_COL32(50, 200,50, 100));   ImPlot::PlotShaded("Carbon", heights.data(), carbonR.data(), N);   ImPlot::PopStyleColor();
                ImPlot::PushStyleColor(ImPlotCol_Fill, IM_COL32(200,130,50, 100));   ImPlot::PlotShaded("Glue", heights.data(), glueR.data(), N);   ImPlot::PopStyleColor();
                ImPlot::PushStyleColor(ImPlotCol_Fill, IM_COL32(50, 100,200,100));   ImPlot::PlotShaded("Steel", heights.data(), steelR.data(), N);   ImPlot::PopStyleColor();
                ImPlot::EndPlot();
            }
            ImGui::End();
            
            ImGui::Begin("Thickness Profile");
            if (ImPlot::BeginPlot("##profile", ImVec2(-1, 300))) {
                // Set up axes
                ImPlot::SetupAxisLimits(ImAxis_Y1, 0.0f, *std::max_element(thermal.begin(), thermal.end()), ImPlotCond_Always);
                ImPlot::SetupAxisLimits(ImAxis_X1, robotLength, 0.0f, ImPlotCond_Always);
                ImPlot::SetupAxis(ImAxis_Y1, "Radius (m)");
                ImPlot::SetupAxis(ImAxis_X1, "Height (m)");

                // Plot steel as a line
                ImPlot::PushStyleColor(ImPlotCol_Line, IM_COL32(50, 100, 200, 255)); // Blue
                ImPlot::PlotLine("Steel", heights.data(), steel.data(), N);
                ImPlot::PopStyleColor();

                // Plot glue as a line
                ImPlot::PushStyleColor(ImPlotCol_Line, IM_COL32(200, 130, 50, 255)); // Orange
                ImPlot::PlotLine("Glue", heights.data(), glue.data(), N);
                ImPlot::PopStyleColor();

                // Plot carbon as a line
                ImPlot::PushStyleColor(ImPlotCol_Line, IM_COL32(50, 200, 50, 255)); // Green
                ImPlot::PlotLine("Carbon", heights.data(), carbon.data(), N);
                ImPlot::PopStyleColor();

                // Plot thermal as a line
                ImPlot::PushStyleColor(ImPlotCol_Line, IM_COL32(200, 50, 200, 255)); // Purple
                ImPlot::PlotLine("TPS", heights.data(), thermal.data(), N);
                ImPlot::PopStyleColor();

                ImPlot::EndPlot();
            }
            ImGui::End();

            // ImGui::Begin("Thickness Profile");
            // if (ImPlot::BeginPlot("##profile", ImVec2(-1, 400))) {
            //     // Set up axes
            //     ImPlot::SetupAxisLimits(ImAxis_X1, 0.0f, 1, ImPlotCond_Always);
            //     ImPlot::SetupAxisLimits(ImAxis_Y1, robotLength, 0.0f, ImPlotCond_Always);
            //     ImPlot::SetupAxis(ImAxis_X1, "Radius (m)");
            //     ImPlot::SetupAxis(ImAxis_Y1, "Height (m)");

            //     // Plot steel as a line
            //     ImPlot::PushStyleColor(ImPlotCol_Line, IM_COL32(50, 100, 200, 255)); // Blue
            //     ImPlot::PlotLine("Steel", steelR.data(), heights.data(), N);
            //     ImPlot::PopStyleColor();

            //     // Plot glue as a line
            //     ImPlot::PushStyleColor(ImPlotCol_Line, IM_COL32(200, 130, 50, 255)); // Orange
            //     ImPlot::PlotLine("Glue", glueR.data(), heights.data(), N);
            //     ImPlot::PopStyleColor();

            //     // Plot carbon as a line
            //     ImPlot::PushStyleColor(ImPlotCol_Line, IM_COL32(50, 200, 50, 255)); // Green
            //     ImPlot::PlotLine("Carbon", carbonR.data(), heights.data(), N);
            //     ImPlot::PopStyleColor();

            //     // Plot thermal as a line
            //     ImPlot::PushStyleColor(ImPlotCol_Line, IM_COL32(200, 50, 200, 255)); // Purple
            //     ImPlot::PlotLine("Thermal", thermalR.data(), heights.data(), N);
            //     ImPlot::PopStyleColor();

            //     ImPlot::EndPlot();
            // }
            // ImGui::End();

        }

        ImGui::Render();

        // ── Camera matrices ───────────────────────────────────────────────────
        glm::vec3 tgt(panX, panY, 0.0f);
        float     yR = glm::radians(camYaw);
        float     pR = glm::radians(camPitch);
        glm::vec3 eye = tgt + glm::vec3(camDistance * cos(pR) * sin(yR),
            camDistance * sin(pR),
            camDistance * cos(pR) * cos(yR));
        glm::mat4 view = glm::lookAt(eye, tgt, glm::vec3(0, 1, 0));
        glm::mat4 proj = glm::perspective(glm::radians(60.0f), fbW / float(fbH), 0.1f, 1000.0f);
        glm::mat4 MVP = proj * view;

        // ── Material colours ---------------------------------------------------
        if (currentViewMode == MODE_MATERIAL) {
            std::vector<glm::vec3> matCols(V.size());
            float minX = 1e9f, maxX = -1e9f;
            for (const auto& v : V) { minX = std::min(minX, v.x); maxX = std::max(maxX, v.x); }

            for (const auto& f : F) {
                float s = ((V[f.v[0]].x + V[f.v[1]].x + V[f.v[2]].x) / 3.0f - minX) / (maxX - minX);
                float tv = 0.0f;
                switch (f.mat) {
                case 0: tv = interp1D(steel, s); break;
                case 1: tv = interp1D(glue, s); break;
                case 2: tv = interp1D(carbon, s); break;
                default: break;
                }
                float tn = std::clamp((tv - thickMin) / (thickMax - thickMin), 0.0f, 1.0f);
                glm::vec3 col = (tn < 0.5f)
                    ? glm::mix(glm::vec3(0, 0, 1), glm::vec3(0, 1, 0), tn * 2.0f)
                    : glm::mix(glm::vec3(0, 1, 0), glm::vec3(1, 0, 0), (tn - 0.5f) * 2.0f);
                matCols[f.v[0]] = matCols[f.v[1]] = matCols[f.v[2]] = col;
            }
            glBindBuffer(GL_ARRAY_BUFFER, CBO);
            glBufferSubData(GL_ARRAY_BUFFER, 0, matCols.size() * sizeof(glm::vec3), matCols.data());
        }
        else {
            // plain grey for Face / Wireframe
            std::vector<glm::vec3> grey(V.size(), glm::vec3(0.8f));
            glBindBuffer(GL_ARRAY_BUFFER, CBO);
            glBufferSubData(GL_ARRAY_BUFFER, 0, grey.size() * sizeof(glm::vec3), grey.data());
        }

        // ── Drawing ───────────────────────────────────────────────────────────
        glViewport(0, 0, fbW, fbH);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glUseProgram(prog);
        glUniformMatrix4fv(locMVP, 1, GL_FALSE, glm::value_ptr(MVP));
        glUniform4f(locClip, 0.0f, 0.0f, 1.0f, -cutPlaneZ);

        glBindVertexArray(VAO);
        glDrawElements(GL_TRIANGLES, (GLsizei)I.size(), GL_UNSIGNED_INT, 0);
        glBindVertexArray(0);

        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(win);
    }

    return 0;
}
