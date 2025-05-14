// mesh_viewer.cpp – single‑file demo that toggles between a 3‑D and 2‑D OBJ.
// model and still supports Face, Wireframe, and Material colouring modes.
//
// Build:.
// g++ mesh_viewer.cpp -I<path‑to‑glad> -I<path‑to‑glm> -lglfw -ldl -lGL -std=c++17 -o viewer.
// (plus the ImGui / ImPlot / glad sources in your build system).
//
// Controls:.
// • Radio buttons   : pick 3‑D vs 2‑D mesh, view mode, and cut plane.
// • W key           : cycles Face/Wireframe/Material if the GUI is hidden.
// • Mouse (LMB)     : orbit.
// • Mouse (RMB)     : pan.
// • Scroll          : dolly zoom.
// --------------------------------------------------------------------------.

#ifdef _WIN32
#define BASE_PATH "./"
#define FONT_SIZE 1.0f
#else
#define BASE_PATH "../"
#define FONT_SIZE 1.5f
#endif

#include "mesh_lib.h"          // simple OBJ + Face data structure (user‑supplied)

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
// #include <unistd.h>.

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "implot.h"

#include "cppcolormap.h"
#include <xtensor/xio.hpp>


// --------------------------------------------------------------------------.
// Enumerations.
// --------------------------------------------------------------------------.

const std::string sim_time = "1";
// 3 visual styles (shading) that apply no matter which mesh is loaded.
enum ViewMode { MODE_FACE = 0, MODE_WIREFRAME, MODE_MATERIAL };
static ViewMode currentViewMode = MODE_FACE;

// Which mesh is currently resident on the GPU.
enum MeshType { MESH_3D = 0, MESH_2D, MESH_TPS };
static MeshType currentMeshType = MESH_3D;

// --------------------------------------------------------------------------.
// Camera / interaction state.
// --------------------------------------------------------------------------.
static float camYaw = 0.0f, camPitch = 0.0f, camDistance = 5.0f;
static float panX = 0.0f, panY = 0.0f;
static bool  leftMouse = false, rightMouse = false;
static double lastX = 0.0, lastY = 0.0;

// --------------------------------------------------------------------------.
// Cutting plane.
// --------------------------------------------------------------------------.
static float cutPlaneZ = 0.0f;   // in model space (metres)

// --------------------------------------------------------------------------.
// THICKNESS PROFILES (loaded from CSV).
// --------------------------------------------------------------------------.
static std::vector<float> carbon, glue, steel;
static float thickMin = 1e9f, thickMax = -1e9f;
static std::vector<float> heights, tps,
steel_metric, glue_metric, carbon_metric, tps_metric,
glue_metric_stacked, carbon_metric_stacked, tps_metric_stacked;
static const float robotLength = 2.5f;                   // m
static std::vector<float> thermZ, thermT;                // thermal profile

// --------------------------------------------------------------------------.
// Opengl resources that need to be recreated whenever we switch mesh.
// --------------------------------------------------------------------------.
static GLuint VAO = 0, VBO = 0, CBO = 0, EBO = 0;
static Mesh   mesh;                       // CPU‑side representation (faces, verts)
static std::vector<glm::vec3>   V;
static std::vector<Face>        F;
static std::vector<unsigned int>I;

// --------------------------------------------------------------------------.
// Forward declarations.
// --------------------------------------------------------------------------.
bool  loadThicknessCSV();                     // called once at start‑up
bool  loadMeshToGPU(const std::string& path); // used every time mesh changes
GLuint compileShader(GLenum, const char*);
GLuint createProgram();

// Really short Functions.
#define Min(arr) *std::min_element(arr.begin(), arr.end());
#define Max(arr) *std::max_element(arr.begin(), arr.end());
#define Normalize(elem, arr) std::clamp((elem - *std::min_element(arr.begin(), arr.end())) / (*std::max_element(arr.begin(), arr.end()) - *std::min_element(arr.begin(), arr.end())), 0.f, 1.f);

// --------------------------------------------------------------------------.
// CSV HELPERS (single‑column and two‑column).
// --------------------------------------------------------------------------.

// makes excutable compatible in ./ or ./build.
std::string PrependBasePath(const std::string& path) {
    return BASE_PATH + path;
}

// Load two-column CSV (position,value).
// std::vector<float> loadProfile1(const std::string& fname) {
//     std::vector<float> vals;
//     std::ifstream in(fname);
//     if (!in) {
//         std::cerr << "[CSV] Cannot open " << fname << "\n";
//         return vals;
//     }
//     std::string line; std::getline(in, line);              // skip header
//     while (std::getline(in, line)) {
//         std::istringstream ss(line);
//         float s, v; char comma;
//         ss >> s >> comma >> v;             // s is the 0…1 normalised height
//         vals.push_back(v);
//         thickMin = std::min(thickMin, v);
//         thickMax = std::max(thickMax, v);
//     }
//     return vals;
// }
std::vector<float> loadProfile1(const std::string& fname) {
    std::vector<float> vals;
    std::ifstream in(fname);
    if (!in) {
        std::cerr << "[CSV] ERROR: Cannot open file " << fname << "\n";
        return vals;
    }

    std::string line;
    if (!std::getline(in, line)) {
        std::cerr << "[CSV] ERROR: File " << fname << " is empty or invalid.\n";
        return vals;
    }

    while (std::getline(in, line)) {
        try {
            std::istringstream ss(line);
            float s, v;
            char comma;

            // Ensure the line has the correct format (e.g., "s,v").
            if (!(ss >> s >> comma >> v) || comma != ',') {
                throw std::runtime_error("Invalid CSV format (expected 's,v').");
            }

            vals.push_back(v);
            thickMin = std::min(thickMin, v);
            thickMax = std::max(thickMax, v);
        }
        catch (const std::exception& e) {
            std::cerr << "[CSV] WARNING: Skipping invalid line: " << line << " (" << e.what() << ")\n";
            continue;
        }
    }

    if (vals.empty()) {
        std::cerr << "[CSV] ERROR: No valid data found in file " << fname << "\n";
    }

    return vals;
}

// static void loadThermal(const std::string& fname)
// {
//     std::ifstream in(fname);
//     if (!in) { std::cerr << "[CSV] Cannot open " << fname << "\n"; return; }
//     std::string line; std::getline(in, line);
//     while (std::getline(in, line)) {
//         std::istringstream ss(line);
//         float z, t; char comma;
//         ss >> z >> comma >> t;
//         thermZ.push_back(z);
//         thermT.push_back(t);
//     }
// }

static void loadThermal(const std::string& fname) {
    std::ifstream in(fname);
    if (!in) {
        std::cerr << "[CSV] ERROR: Cannot open file " << fname << "\n";
        return;
    }

    std::string line;
    if (!std::getline(in, line)) {
        std::cerr << "[CSV] ERROR: File " << fname << " is empty or invalid.\n";
        return;
    }

    while (std::getline(in, line)) {
        try {
            std::istringstream ss(line);
            float z, t;
            char comma;

            // Ensure the line has the correct format (e.g., "z,t").
            if (!(ss >> z >> comma >> t) || comma != ',') {
                throw std::runtime_error("Invalid CSV format (expected 'z,t').");
            }

            thermZ.push_back(z);
            thermT.push_back(t);
        } catch (const std::exception& e) {
            std::cerr << "[CSV] WARNING: Skipping invalid line: " << line << " (" << e.what() << ")\n";
            continue;
        }
    }

    if (thermZ.empty() || thermT.empty()) {
        std::cerr << "[CSV] ERROR: No valid data found in file " << fname << "\n";
    }
}

// Linear look‑ups for material and thermal profiles.
// s shold be float between 0.0 and 1.0, and will be clamped to that interval if it exceeds.
static float interp1D(const std::vector<float>& arr, float s)
{
    s = std::clamp(s, 0.f, 1.f);
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

// Color helper functions.
glm::vec3 ImU32ToVec3(ImU32 color) {
    float r = ((color >> IM_COL32_R_SHIFT) & 0xFF) / 255.0f; // Extract red
    float g = ((color >> IM_COL32_G_SHIFT) & 0xFF) / 255.0f; // Extract green
    float b = ((color >> IM_COL32_B_SHIFT) & 0xFF) / 255.0f; // Extract blue
    return glm::vec3(r, g, b); // Return as glm::vec3
}

ImU32 Vec3ToImU32(const glm::vec3& color, float alpha = 1.0f) {
    int r = static_cast<int>(std::clamp(color.r, 0.0f, 1.0f) * 255.0f);
    int g = static_cast<int>(std::clamp(color.g, 0.0f, 1.0f) * 255.0f);
    int b = static_cast<int>(std::clamp(color.b, 0.0f, 1.0f) * 255.0f);
    int a = static_cast<int>(std::clamp(alpha, 0.0f, 1.0f) * 255.0f);
    return IM_COL32(r, g, b, a);
}

static glm::vec3 InterpolateColormapToVec3(const xt::xtensor<double, 2>& colormap, float s) {
    // Clamp s to the range [0.0, 1.0].
    s = std::clamp(s, 0.0f, 1.0f);

    // Map s to an index in the colormap.
    size_t index = static_cast<size_t>(s * (colormap.shape(0) - 1));

    // Extract the RGB values (normalized to [0.0, 1.0]).
    float r = static_cast<float>(colormap(index, 0));
    float g = static_cast<float>(colormap(index, 1));
    float b = static_cast<float>(colormap(index, 2));

    // Convert to ImU32 (RGBA format, alpha = 255).
    return glm::vec3(r, g, b);
}

// --------------------------------------------------------------------------.
// Mesh loading + gpu buffer upload.
// --------------------------------------------------------------------------.

bool loadMeshToGPU(const std::string& path)
{
    // 1) Load into your CPU‐side mesh object.
    if (!mesh.loadOBJ(path)) {
        std::cerr << "[OBJ] Failed to load " << path << "\n";
        return false;
    }

    // 2) Copy raw vertices & faces.
    V = mesh.vertices;
    F = mesh.faces;

    // ──────────────────────────────────────────────────────────────────────────.
    // 3) FILTER OUT BAD FACES (guard against out-of-range indices).
    std::vector<Face> filteredF;
    filteredF.reserve(F.size());
    for (auto& f : F) {
        if (f.v[0] < V.size() && f.v[1] < V.size() && f.v[2] < V.size()) {
            filteredF.push_back(f);
        }
        else {
            std::cerr << "[loadMeshToGPU] dropping face with bad indices: ("
                << f.v[0] << "," << f.v[1] << "," << f.v[2] << ")\n";
        }
    }
    F.swap(filteredF);
    // ──────────────────────────────────────────────────────────────────────────.

    // 4) Rebuild your index buffer from the sanitized face list.
    I.clear();
    I.reserve(F.size() * 3);
    for (auto& f : F) {
        I.push_back(f.v[0]);
        I.push_back(f.v[1]);
        I.push_back(f.v[2]);
    }

    // 5) (Re)upload VBO, CBO, EBO exactly as before….
    if (VAO == 0) {
        glGenVertexArrays(1, &VAO);
        glGenBuffers(1, &VBO);
        glGenBuffers(1, &CBO);
        glGenBuffers(1, &EBO);
    }

    glBindVertexArray(VAO);

    // vertex positions.
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, V.size() * sizeof(glm::vec3), V.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
    glEnableVertexAttribArray(0);

    // default grey colours.
    std::vector<glm::vec3> defCols(V.size(), glm::vec3(0.8f));
    glBindBuffer(GL_ARRAY_BUFFER, CBO);
    glBufferData(GL_ARRAY_BUFFER, defCols.size() * sizeof(glm::vec3), defCols.data(), GL_DYNAMIC_DRAW);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
    glEnableVertexAttribArray(1);

    // triangle indices.
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, I.size() * sizeof(unsigned int), I.data(), GL_STATIC_DRAW);

    glBindVertexArray(0);
    return true;
}

void frameMeshOnLoad()
{
    // compute axis-aligned bounding box of the newly loaded mesh.
    glm::vec3 mn(1e9f), mx(-1e9f);
    for (auto& v : mesh.vertices) {
        mn = glm::min(mn, v);
        mx = glm::max(mx, v);
    }
    // center in X,Y (we treat Z as depth for 3D).
    glm::vec3 center = (mn + mx) * 0.5f;
    panX = center.x;
    panY = center.y;
    // radius = half the diagonal.
    float radius = glm::length(mx - mn) * 0.5f;
    if (radius < 1e-3f) radius = 1.0f;
    camDistance = radius * 3.0f;    // 3× so you see it comfortably
    camYaw = 0.0f;
    camPitch = 0.0f;
}

// --------------------------------------------------------------------------.
// SHADERS (single clip plane).
// --------------------------------------------------------------------------.
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
const std::string OBJ_TPS = PrependBasePath("assets/obj/TPS_" + sim_time + "_hr.obj");

// --------------------------------------------------------------------------.
// Glfw callbacks.
// --------------------------------------------------------------------------.
static void key_cb(GLFWwindow* window, int key, int /*scancode*/, int action, int /*mods*/)
{
    if (action != GLFW_PRESS)
        return;

    if (key == GLFW_KEY_W)
    {
        // 1) cycle Face→Wireframe→Material.
        currentViewMode = ViewMode((currentViewMode + 1) % 3);

        // 2) set OpenGL fill/line.
        GLenum mode = (currentViewMode == MODE_WIREFRAME ? GL_LINE : GL_FILL);
        glPolygonMode(GL_FRONT_AND_BACK, mode);

        // 3) if we just switched into Material, force-load the TPS-2D mesh.
        if (currentViewMode == MODE_MATERIAL && currentMeshType != MESH_TPS)
        {
            if (loadMeshToGPU(OBJ_TPS))
            {
                currentMeshType = MESH_TPS;
                frameMeshOnLoad();
            }
            else
            {
                std::cerr << "[Mesh] failed to load TPS-2D on W→Material\n";
            }
        }
    }
    else if (key == GLFW_KEY_F)
    {
        frameMeshOnLoad();
    }
    else if (key == GLFW_KEY_ESCAPE || key == GLFW_KEY_Q)
    {
        glfwSetWindowShouldClose(window, GLFW_TRUE);
    }
}


static void mouse_btn_cb(GLFWwindow*, int button, int action, int)
{
    // if mouse is over any ImGui window or plot, do not start/stop camera drag.
    if (ImGui::GetIO().WantCaptureMouse)
        return;

    if (button == GLFW_MOUSE_BUTTON_LEFT)
        leftMouse = (action == GLFW_PRESS);
    if (button == GLFW_MOUSE_BUTTON_RIGHT)
        rightMouse = (action == GLFW_PRESS);
}
static void cursor_cb(GLFWwindow*, double x, double y)
{
    // ignore camera motion when interacting with UI.
    if (ImGui::GetIO().WantCaptureMouse)
        return;

    double dx = x - lastX;
    double dy = y - lastY;
    lastX = x;
    lastY = y;

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
    // ignore zoom when over UI.
    if (ImGui::GetIO().WantCaptureMouse)
        return;

    camDistance *= (1.0f - float(yoff) * 0.1f);
    camDistance = std::max(camDistance, 0.1f);
}


// --------------------------------------------------------------------------.
// Main.
// --------------------------------------------------------------------------.
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

    // ── ImGui + ImPlot ─────────────────────────────────────────────────────────.
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.FontGlobalScale = FONT_SIZE;
    ImPlot::CreateContext();
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(win, true);
    ImGui_ImplOpenGL3_Init("#version 330 core");

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // ------------------------------------------------------------------ csv ----.
    carbon = loadProfile1(PrependBasePath("assets/csv/carbon_thickness.csv"));
    glue = loadProfile1(PrependBasePath("assets/csv/glue_thickness.csv"));
    steel = loadProfile1(PrependBasePath("assets/csv/steel_thickness.csv"));
    loadThermal(PrependBasePath("assets/csv/Thickness_" + sim_time + "_hr.csv"));

    // Check material vectors
    if (carbon.empty() || glue.empty() || steel.empty()) {
        std::cerr << "[CSV] ERROR: One or more thickness CSV files are empty.\n";
        return 1; // Exit or handle the error
    }
    if (carbon.size() != glue.size() || glue.size() != steel.size()) {
        std::cerr << "[CSV] ERROR: Thickness CSV files have inconsistent sizes.\n";
        return 1; // Exit or handle the error
    }

    
    // Load CSV profiles from assets/csv.
    int N = (int)carbon.size();
    heights.resize(N);
    tps.resize(N);
    steel_metric.resize(N);
    glue_metric.resize(N);
    carbon_metric.resize(N);
    tps_metric.resize(N);
    glue_metric_stacked.resize(N);
    carbon_metric_stacked.resize(N);
    tps_metric_stacked.resize(N);

    // calculate thicknesses.
    for (int i = 0; i < N; ++i) {
        float s = i / float(N - 1);
        heights[i] = s * robotLength;
        tps[i] = interpThermal(heights[i]);
        steel_metric[i] = steel[i] / 100.0f;
        glue_metric[i] = glue[i] / 100.0f;
        carbon_metric[i] = carbon[i] / 100.0f;
        tps_metric[i] = tps[i] / 100.0f;
        glue_metric_stacked[i] = steel_metric[i] + glue_metric[i];
        carbon_metric_stacked[i] = glue_metric_stacked[i] + carbon_metric[i];
        tps_metric_stacked[i] = carbon_metric_stacked[i] + tps_metric[i];
    }

    // cmap constants.
    static const std::string cmap_name = "GnBu_r";
    static const  xt::xtensor<double, 2> cmap = cppcolormap::colormap(cmap_name, 256);
    // ----------------------------------------------------------------- obj ----.
    const std::string OBJ_3D = PrependBasePath("assets/obj/humanoid_robot.obj");
    const std::string OBJ_2D = PrependBasePath("assets/obj/humanoid_robot_2d.obj");
    const std::string OBJ_TPS = PrependBasePath("assets/obj/TPS_" + sim_time + "_hr.obj");

    if (!loadMeshToGPU(OBJ_3D)) return 1;          // start with 3‑D model

    // ----------------------------------------------------------------- glsl ----.

    GLuint prog = createProgram();
    GLint  locMVP = glGetUniformLocation(prog, "uMVP");
    GLint  locClip = glGetUniformLocation(prog, "uClipPlane");

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CLIP_DISTANCE0);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    // ---------------------------------------------------------------- main loop -.
    while (!glfwWindowShouldClose(win)) {
        glfwPollEvents();
        int fbW, fbH; glfwGetFramebufferSize(win, &fbW, &fbH);
        fbH = std::max(fbH, 1);

        // ── ImGui frame ───────────────────────────────────────────────────────.
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        ImGui::SetNextWindowSize(ImVec2(280, 0), ImGuiCond_Once);
        ImGui::Begin("Controls");

        // --- Mesh selection ---------------------------------------------------.

        ImGui::BeginDisabled(currentViewMode == MODE_MATERIAL);
        ImGui::Text("Mesh");
        if (ImGui::RadioButton("3D", currentMeshType == MESH_3D)) {
            loadMeshToGPU(OBJ_3D);
            currentMeshType = MESH_3D;
        }
        ImGui::SameLine();
        if (ImGui::RadioButton("2D", currentMeshType == MESH_2D)) {
            loadMeshToGPU(OBJ_2D);
            currentMeshType = MESH_2D;
        }
        ImGui::SameLine();
        if (ImGui::RadioButton("Visualize TPS", currentMeshType == MESH_TPS)) {
            loadMeshToGPU(OBJ_TPS);
            currentMeshType = MESH_TPS;
        }
        ImGui::EndDisabled();

        ImGui::Separator();
        if (ImGui::Button("Zoom to Fit (F)")) { frameMeshOnLoad(); }
        ImGui::Text("View Mode (Switch using 'W')");
        if (ImGui::RadioButton("Face", currentViewMode == MODE_FACE)) {
            currentViewMode = MODE_FACE;
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        }
        if (ImGui::RadioButton("Wireframe", currentViewMode == MODE_WIREFRAME)) {
            currentViewMode = MODE_WIREFRAME;
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        }
        if (ImGui::RadioButton("Material", currentViewMode == MODE_MATERIAL)) {
            currentViewMode = MODE_MATERIAL;
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            // automatically switch mesh to TPS-2D:.
            if (currentMeshType != MESH_TPS) {
                loadMeshToGPU(OBJ_TPS);
                currentMeshType = MESH_TPS;
                frameMeshOnLoad();
            }
        }
        ImGui::Separator();
        ImGui::SliderFloat("Cut Z", &cutPlaneZ, -1.0f, 1.0f, "%.2f m");
        ImGui::End();

        if (currentViewMode == MODE_MATERIAL) {
            ImGuiIO& io = ImGui::GetIO();
            float windowWidth = io.DisplaySize.x;
            float windowHeight = io.DisplaySize.y;

            // Dynamically position the colorbar on the right side.
            ImGui::SetNextWindowPos(ImVec2(windowWidth - 100, 50), ImGuiCond_Always); // 100px from the right
            ImGui::SetNextWindowSize(ImVec2(80, windowHeight - 100), ImGuiCond_Always); // 80px wide, dynamic height

            ImGui::Begin("Colorbar", nullptr, ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoCollapse);
            // ImGui::SetNextWindowPos(ImVec2(50, 50), ImGuiCond_Once); // Set position.
            // ImGui::Begin("Colorbar");.


            // Define the size and position of the colorbar.
            ImVec2 barSize = ImVec2(40, 200); // Width: 20px, Height: 200px
            ImVec2 barPos = ImGui::GetCursorScreenPos(); // Position at the current cursor
            barPos.y += 2 * ImGui::GetTextLineHeight();

            // Draw the gradient colorbar.
            ImDrawList* drawList = ImGui::GetWindowDrawList();
            int gradient_segments = 256; // Increase for smoother gradient

            for (int i = 0; i < gradient_segments; ++i) {
                float t = 1 - (float)i / (gradient_segments - 1); // Normalized value [0, 1]
                ImU32 color = Vec3ToImU32(InterpolateColormapToVec3(cmap, t));

                float yStart = barPos.y + i * (barSize.y / gradient_segments);
                float yEnd = barPos.y + (i + 1) * (barSize.y / gradient_segments);
                drawList->AddRectFilled(ImVec2(barPos.x, yStart), ImVec2(barPos.x + barSize.x, yEnd), color);

            }

            // Draw a border around the colorbar.
            // drawList->AddRect(barPos, ImVec2(barPos.x + barSize.x, barPos.y + barSize.y), IM_COL32(255, 255, 255, 255), 0.0f, 0, 2.0f);.

            // Add labels for min and max values.
            ImGui::SetCursorScreenPos(ImVec2(barPos.x, barPos.y - 2 * ImGui::GetTextLineHeight()));
            ImGui::Text("Max:\n%.2f", *std::max_element(tps_metric_stacked.begin(), tps_metric_stacked.end()));
            ImGui::SetCursorScreenPos(ImVec2(barPos.x, barPos.y + barSize.y));
            ImGui::Text("Min:\n%.2f", *std::min_element(tps_metric_stacked.begin(), tps_metric_stacked.end()));

            ImGui::End();

            // ---------------------------------------------------------------------.
            // Thickness stacked‑area plot if Material mode is active.
            ImU32 magenta = IM_COL32(200, 50, 200, 255);
            ImU32 green = IM_COL32(50, 200, 50, 255);
            ImU32 orange = IM_COL32(200, 130, 50, 255);
            ImU32 blue = IM_COL32(50, 100, 200, 255);

            int plot_height = 250;


            ImGui::SetNextWindowPos(ImVec2(0, windowHeight - (2 * plot_height + 50)), ImGuiCond_Appearing); // 100px from the right
            ImGui::SetNextWindowSize(ImVec2(windowWidth, (2 * plot_height + 50)), ImGuiCond_Appearing);
            // ImGui::SetNextWindowSize(ImVec2(windowWidth, 2 * plot_height), ImGuiCond_Once);.
            ImGui::Begin("Thickness Profile");
            if (ImPlot::BeginPlot("Thickness Stacked", ImVec2(-1, plot_height))) {


                ImPlot::SetupAxisLimits(ImAxis_Y1, 0.0f, *std::max_element(tps_metric_stacked.begin(), tps_metric_stacked.end()), ImPlotCond_Always);
                ImPlot::SetupAxisLimits(ImAxis_X1, robotLength, 0.0f, ImPlotCond_Always);
                ImPlot::SetupAxis(ImAxis_Y1, "Thickness (m)");
                // ImPlot::SetupAxis(ImAxis_X1, "Height (m)", ImPlotAxisFlags_NoDecorations);.
                ImPlot::SetupAxis(ImAxis_X1, "Height (m)");

                ImPlot::PushStyleColor(ImPlotCol_Fill, magenta);
                ImPlot::PlotShaded("TPS", heights.data(), tps_metric_stacked.data(), N);
                ImPlot::PopStyleColor();

                ImPlot::PushStyleColor(ImPlotCol_Fill, green);
                ImPlot::PlotShaded("Carbon", heights.data(), carbon_metric_stacked.data(), N);
                ImPlot::PopStyleColor();

                ImPlot::PushStyleColor(ImPlotCol_Fill, orange);
                ImPlot::PlotShaded("Glue", heights.data(), glue_metric_stacked.data(), N);
                ImPlot::PopStyleColor();

                ImPlot::PushStyleColor(ImPlotCol_Fill, blue);
                ImPlot::PlotShaded("Steel", heights.data(), steel_metric.data(), N);
                ImPlot::PopStyleColor();

                ImPlot::EndPlot();
            }
            ImGui::End();

            ImGui::Begin("Thickness Profile");
            if (ImPlot::BeginPlot("Thickness Separate", ImVec2(-1, plot_height))) {
                // Set up axes.
                ImPlot::SetupAxisLimits(ImAxis_Y1, 0.0f, *std::max_element(tps_metric.begin(), tps_metric.end()), ImPlotCond_Always);
                ImPlot::SetupAxisLimits(ImAxis_X1, robotLength, 0.0f, ImPlotCond_Always);
                ImPlot::SetupAxis(ImAxis_Y1, "Thickness (m)");
                ImPlot::SetupAxis(ImAxis_X1, "Height (m)");

                float line_thickness = 2.0f;
                ImPlot::PushStyleVar(ImPlotStyleVar_LineWeight, line_thickness);
                // Plot thermal as a line.
                ImPlot::PushStyleColor(ImPlotCol_Line, magenta); // Purple
                ImPlot::PlotLine("TPS", heights.data(), tps_metric.data(), N);
                ImPlot::PopStyleColor();

                // Plot carbon as a line.
                ImPlot::PushStyleColor(ImPlotCol_Line, green); // Green
                ImPlot::PlotLine("Carbon", heights.data(), carbon_metric.data(), N);
                ImPlot::PopStyleColor();

                // Plot glue as a line.
                ImPlot::PushStyleColor(ImPlotCol_Line, orange); // Orange
                ImPlot::PlotLine("Glue", heights.data(), glue_metric.data(), N);
                ImPlot::PopStyleColor();

                // Plot steel as a line.
                ImPlot::PushStyleColor(ImPlotCol_Line, blue); // Blue
                ImPlot::PlotLine("Steel", heights.data(), steel_metric.data(), N);
                ImPlot::PopStyleColor();
                ImPlot::PopStyleVar();

                ImPlot::EndPlot();
            }
            ImGui::End();
        }

        ImGui::Render();

        // ── Camera matrices ───────────────────────────────────────────────────.
        glm::vec3 tgt(panX, panY, 0.0f);
        float     yR = glm::radians(camYaw);
        float     pR = glm::radians(camPitch);
        glm::vec3 eye = tgt + glm::vec3(camDistance * cos(pR) * sin(yR),
            camDistance * sin(pR),
            camDistance * cos(pR) * cos(yR));
        glm::mat4 view = glm::lookAt(eye, tgt, glm::vec3(0, 1, 0));
        glm::mat4 proj = glm::perspective(glm::radians(60.0f), fbW / float(fbH), 0.1f, 1000.0f);
        glm::mat4 MVP = proj * view;

        // ——————————————————————————————————————————————.
        // MATERIAL COLOURS (replace your entire original block).
        // ——————————————————————————————————————————————.
        if (currentViewMode == MODE_MATERIAL) {
            std::vector<glm::vec3> vertex_colors(V.size());
            float minX = 1e9f, maxX = -1e9f;
            float minY = 1e9f, maxY = -1e9f;
            float minZ = 1e9f, maxZ = -1e9f;
            for (const auto& v : V) {
                minX = std::min(minX, v.x);
                maxX = std::max(maxX, v.x);
                minY = std::min(minY, v.y);
                maxY = std::max(maxY, v.y);
                minZ = std::min(minZ, v.z);
                maxZ = std::max(maxX, v.z);
            }

            for (int i = 0; i < V.size(); i++) {
                const auto& v = V[i];
                float normalized_y_at_v = (v.y - minY) / (maxY - minY);
                float tps_at_v = interp1D(tps_metric_stacked, normalized_y_at_v);

                float normalized_tps_at_v = Normalize(tps_at_v, tps_metric_stacked);
                vertex_colors[i] = InterpolateColormapToVec3(cmap, normalized_tps_at_v);
            }

            glBindBuffer(GL_ARRAY_BUFFER, CBO);
            glBufferSubData(GL_ARRAY_BUFFER, 0, vertex_colors.size() * sizeof(glm::vec3), vertex_colors.data());
        }
        else {
            // your grey‐only fallback (unchanged).
            std::vector<glm::vec3> grey(V.size(), glm::vec3(0.8f));
            glBindBuffer(GL_ARRAY_BUFFER, CBO);
            glBufferSubData(
                GL_ARRAY_BUFFER,
                0,
                (GLsizeiptr)(grey.size() * sizeof(glm::vec3)),
                grey.data()
            );
        }


        // ── Drawing Mesh ───────────────────────────────────────────────────────────.
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
