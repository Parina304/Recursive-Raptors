#include "mesh_lib.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <glm/glm.hpp>

// Helper to trim whitespace from a string (for parsing)
static std::string trim(const std::string& s) {
    size_t start = s.find_first_not_of(" \t\r\n");
    size_t end = s.find_last_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    return s.substr(start, end - start + 1);
}

bool Mesh::loadOBJ(const std::string& filePath) {
    std::ifstream file(filePath);
    if (!file) {
        std::cerr << "ERROR: Could not open OBJ file " << filePath << std::endl;
        return false;
    }
    std::string line;
    std::string currentMatName = "";
    MaterialType currentMat = STEEL; // default to steel if none specified
    // Temporary arrays for parsing
    std::vector<glm::vec3> tempVerts;
    // (We ignore texture coords and normals in this loader for simplicity)
    minX = minY = minZ = 1e9f;
    maxX = maxY = maxZ = -1e9f;
    while (std::getline(file, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;
        if (line.rfind("v ", 0) == 0) {
            // Vertex line
            std::istringstream iss(line.substr(2));
            glm::vec3 v;
            iss >> v.x >> v.y >> v.z;
            tempVerts.push_back(v);
            // update bounds
            if (v.x < minX) minX = v.x;
            if (v.x > maxX) maxX = v.x;
            if (v.y < minY) minY = v.y;
            if (v.y > maxY) maxY = v.y;
            if (v.z < minZ) minZ = v.z;
            if (v.z > maxZ) maxZ = v.z;
        }
        else if (line.rfind("usemtl ", 0) == 0) {
            // Material usage for subsequent faces
            currentMatName = trim(line.substr(7));
            // Convert material name to enum
            currentMat = materialTypeFromName(currentMatName);
        }
        else if (line.rfind("f ", 0) == 0) {
            // Face line (triangle assumed). OBJ faces are 1-indexed.
            // We also handle faces with texture/normal: "f v/t/n v/t/n v/t/n"
            std::istringstream iss(line.substr(2));
            std::string v1, v2, v3;
            if (!(iss >> v1 >> v2 >> v3)) {
                continue; // skip if parse fails
            }
            // Parse each vertex index group up to '/', ignoring texture/normal indices.
            auto parseIndex = [](const std::string& token) {
                // token might be "int" or "int/int/int" etc.
                size_t slashPos = token.find('/');
                int idx = 0;
                if (slashPos == std::string::npos) {
                    // no slash
                    idx = std::stoi(token);
                }
                else {
                    idx = std::stoi(token.substr(0, slashPos));
                }
                return idx;
                };
            int i1 = parseIndex(v1);
            int i2 = parseIndex(v2);
            int i3 = parseIndex(v3);
            Face face;
            face.v[0] = i1 - 1; // convert to 0-index
            face.v[1] = i2 - 1;
            face.v[2] = i3 - 1;
            face.mat = currentMat;
            // Assign an initial temperature based on material as a simple heuristic:
            float baseTemp = 300.0f; // base ambient temperature in K
            switch (face.mat) {
            case INSULATION: face.temp = baseTemp + 10.0f; break;
            case CFRP:       face.temp = baseTemp + 40.0f; break;  // outer layer hotter
            case GLUE:       face.temp = baseTemp + 30.0f; break;
            case STEEL:      face.temp = baseTemp + 20.0f; break;
            }
            faces.push_back(face);
        }
        // (We ignore other OBJ statements like mtllib, normals, etc., for brevity)
    }
    file.close();
    // Move unique vertices into Mesh::vertices
    vertices = std::move(tempVerts);
    // Determine initial global minTemp and maxTemp
    minTemp = 1e9f;
    maxTemp = -1e9f;
    for (const Face& f : faces) {
        if (f.temp < minTemp) minTemp = f.temp;
        if (f.temp > maxTemp) maxTemp = f.temp;
    }
    // Prepare flat vertex/color arrays for rendering
    createFlatArrays();
    CalcStats();
    CalcFaceCentroid();
    CalcFaceNormal();

    return true;
}

void Mesh::createFlatArrays() {
    flatVertices.clear();
    flatColors.clear();
    flatVertices.reserve(faces.size() * 3);
    flatColors.reserve(faces.size() * 3);
    for (size_t i = 0; i < faces.size(); ++i) {
        const Face& face = faces[i];
        // Color will be set later per mode, but initialize to white (or material color as default)
        glm::vec3 faceColor = glm::vec3(1.0f, 1.0f, 1.0f);
        // Or we could initialize flatColor to the material color initially:
        faceColor = getMaterialColor(face.mat);
        for (int j = 0; j < 3; ++j) {
            int vIndex = face.v[j];
            flatVertices.push_back(vertices[vIndex]);
            flatColors.push_back(faceColor);
        }
    }
}

MaterialType Mesh::materialTypeFromName(const std::string& name) const {
    // Case-insensitive match for material names
    std::string lower;
    lower.reserve(name.size());
    for (char c : name) lower.push_back(std::tolower(c));
    if (lower.find("insulation") != std::string::npos || lower.find("protection") != std::string::npos) {
        return INSULATION;
    }
    if (lower.find("glue") != std::string::npos) {
        return GLUE;
    }
    if (lower.find("steel") != std::string::npos) {
        return STEEL;
    }
    if (lower.find("cfrp") != std::string::npos || lower.find("carbon") != std::string::npos) {
        return CFRP;
    }
    // Default to steel if unknown
    return STEEL;
}

glm::vec3 Mesh::getMaterialColor(MaterialType mat) const {
    switch (mat) {
    case INSULATION: return glm::vec3(1.0f, 0.5f, 0.0f);   // orange
    case GLUE:       return glm::vec3(0.2f, 0.8f, 0.2f);   // green
    case STEEL:      return glm::vec3(0.2f, 0.2f, 0.8f);   // blue
    case CFRP:       return glm::vec3(0.8f, 0.8f, 0.2f);   // yellow
    }
    return glm::vec3(1.0f); // default white
}

glm::vec3 Mesh::getColorForTemperature(float value) const {
    // Normalize value between minTemp and maxTemp
    float t = 0.0f;
    if (maxTemp > minTemp) {
        t = (value - minTemp) / (maxTemp - minTemp);
        if (t < 0.0f) t = 0.0f;
        if (t > 1.0f) t = 1.0f;
    }
    // Use HSV interpolation from blue (HSV hue=240) to red (HSV hue=0)
    // S=1, V=1 for vivid colors.
    float hue = (1.0f - t) * 240.0f; // 0 -> red, 120 -> green, 240 -> blue
    float H = hue / 60.0f;
    int i = (int)floor(H);
    float f = H - i;
    float p = 0.0f;
    float q = 1.0f - f;
    float r, g, b;
    switch (i) {
    case 0: r = 1.0f; g = f;    b = 0.0f; break;      // red to yellow
    case 1: r = q;    g = 1.0f; b = 0.0f; break;      // yellow to green
    case 2: r = 0.0f; g = 1.0f; b = f;    break;      // green to cyan
    case 3: r = 0.0f; g = q;    b = 1.0f; break;      // cyan to blue
    case 4: r = f;    g = 0.0f; b = 1.0f; break;      // blue to (should be magenta, but we stop at 240 so i=4 covers 240->300 deg not used)
    default: r = 1.0f; g = 0.0f; b = 0.0f; break;     // fall-back to red
    }
    return glm::vec3(r, g, b);
}

void Mesh::increaseProtectionThickness() {
    // Simulate adding a thermal protection layer: turn all CFRP faces into Insulation
    for (Face& face : faces) {
        if (face.mat == CFRP) {
            face.mat = INSULATION;
            // Optionally adjust temperature down due to added insulation
            face.temp -= 10.0f;
            if (face.temp < minTemp) minTemp = face.temp;
            // (We keep maxTemp same or adjust if needed, but reducing a bit won't hurt)
        }
    }
    // Recompute temperature bounds after modification
    minTemp = 1e9f;
    maxTemp = -1e9f;
    for (const Face& f : faces) {
        if (f.temp < minTemp) minTemp = f.temp;
        if (f.temp > maxTemp) maxTemp = f.temp;
    }
    // Update flatColors to reflect new material assignment (if Material mode is active next render)
    for (size_t i = 0; i < faces.size(); ++i) {
        glm::vec3 newCol = getMaterialColor(faces[i].mat);
        // flatColors has 3 entries per face (one per vertex)
        flatColors[i * 3 + 0] = newCol;
        flatColors[i * 3 + 1] = newCol;
        flatColors[i * 3 + 2] = newCol;
    }
    std::cout << "Thermal protection layer added: outer CFRP faces now Insulation.\n";
}

void Mesh::solveThermal() {
    // Stub for simulation: modify temperatures to simulate heating/cooling.
    if (!faces.empty()) {
        bool protectionApplied = false;
        // Check if any face is labeled as Insulation to infer if protection was added
        for (const Face& f : faces) {
            if (f.mat == INSULATION) { protectionApplied = true; break; }
        }
        if (protectionApplied) {
            // If insulation is added, assume overall temperatures drop (improved thermal protection)
            for (Face& f : faces) {
                f.temp *= 0.9f; // reduce temperature by 10%
            }
            std::cout << "Solved thermal (with protection): temperatures reduced.\n";
        }
        else {
            // Without added protection, assume temperatures rise over time
            for (Face& f : faces) {
                // Increase temperatures of each face slightly, more so for outer materials
                if (f.mat == CFRP)       f.temp *= 1.10f;
                else if (f.mat == GLUE)  f.temp *= 1.05f;
                else if (f.mat == STEEL) f.temp *= 1.02f;
                // Insulation (if present initially) might remain cooler
            }
            std::cout << "Solved thermal (no protection): temperatures increased.\n";
        }
        // Recompute min/max temperature after solving
        minTemp = 1e9f;
        maxTemp = -1e9f;
        for (const Face& f : faces) {
            if (f.temp < minTemp) minTemp = f.temp;
            if (f.temp > maxTemp) maxTemp = f.temp;
        }
    }
}

// Saves the mesh in a Binary PLY format for efficient storage
bool Mesh::savePLY(const std::string& filename) const {
    std::ofstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "Error: Cannot write PLY file " << filename << "\n";
        return false;
    }
    file << "ply\n";
    file << "format binary_little_endian 1.0\n";
    file << "element vertex " << vertices.size() << "\n";
    file << "property float x\nproperty float y\nproperty float z\n";
    file << "element face " << faces.size() << "\n";
    file << "property list uchar int vertex_indices\n";
    file << "end_header\n";

    file.write(reinterpret_cast<const char*>(vertices.data()), vertices.size() * sizeof(Vertex));

    for (const auto& face : faces) {
        uint8_t vertexCount = static_cast<uint8_t>(face.vertexIndices.size());
        file.write(reinterpret_cast<const char*>(&vertexCount), sizeof(uint8_t));
        file.write(reinterpret_cast<const char*>(face.vertexIndices.data()), vertexCount * sizeof(uint32_t));
    }

    file.close();
    return true;
}

// Reads a Binary PLY file and reconstructs the mesh in memory
bool Mesh::loadPLY(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "Error: Cannot open PLY file " << filename << "\n";
        return false;
    }

    std::string line;
    bool headerParsed = false;
    uint32_t vertexCount = 0, faceCount = 0;

    // Read Header
    while (std::getline(file, line)) {
        if (line.find("element vertex") != std::string::npos) {
            vertexCount = std::stoi(line.substr(15));
        } else if (line.find("element face") != std::string::npos) {
            faceCount = std::stoi(line.substr(13));
        } else if (line == "end_header") {
            headerParsed = true;
            break;
        }
    }

    if (!headerParsed) {
        std::cerr << "Error: Invalid PLY file format\n";
        return false;
    }

    // Read Binary Vertex Data
    vertices.resize(vertexCount);
    file.read(reinterpret_cast<char*>(vertices.data()), vertexCount * sizeof(Vertex));

    // Read Binary Face Data
    faces.resize(faceCount);
    for (uint32_t i = 0; i < faceCount; i++) {
        uint8_t vertexCount;
        file.read(reinterpret_cast<char*>(&vertexCount), sizeof(uint8_t));

        faces[i].vertexIndices.resize(vertexCount);
        file.read(reinterpret_cast<char*>(faces[i].vertexIndices.data()), vertexCount * sizeof(uint32_t));
    }

    file.close();
    
    // Calculates centroid
    CalcStats();
    CalcFaceNormal();

    return true;
}

// Prints basic mesh statistics
void Mesh::printMeshStats() const {
    std::cout << "Mesh Statistics:\n";
    std::cout << "  - Vertices: " << vertices.size() << "\n";
    std::cout << "  - Faces: " << faces.size() << "\n";
    std::cout << "  - Centroid: x: " << centroid.x << " y: " << centroid.y << " z: " << centroid.z << "\n";
    std::cout << "  - Bounds:\n";
    std::cout << "      x_min: " << x_min << ", x_max: " << x_max << "\n";
    std::cout << "      y_min: " << y_min << ", y_max: " << y_max << "\n";
    std::cout << "      z_min: " << z_min << ", z_max: " << z_max << "\n";
}

void Mesh::CalcStats(){
    float cx = 0, cy = 0, cz = 0;
    for(const auto& v: vertices){
        cx += v.x;
        cy += v.y;
        cz += v.z;
        x_min = std::min(x_min, v.x);
        y_min = std::min(y_min, v.y);
        z_min = std::min(z_min, v.z);
        x_max = std::max(x_max, v.x);
        y_max = std::max(y_max, v.y);
        z_max = std::max(z_max, v.z);
    }
    centroid.x = cx / vertices.size();
    centroid.y = cy / vertices.size();
    centroid.z = cz / vertices.size();
}

// calculates the centroid for each face
void Mesh::CalcFaceCentroid(){
    Vertex centroid;
    for(auto& f: faces){
        for(const auto& v_idx: f.vertexIndices){
            f.centroid = f.centroid + vertices[v_idx];
        }
        f.centroid = f.centroid / f.vertexIndices.size();
    }
}