#include "mesh_lib.h"

Mesh::Mesh(){
    x_min = std::numeric_limits<float>::max();
    y_min = std::numeric_limits<float>::max();
    z_min = std::numeric_limits<float>::max();
    x_max = std::numeric_limits<float>::min();
    y_max = std::numeric_limits<float>::min();
    z_max = std::numeric_limits<float>::min();
}

// Reads an OBJ file and stores its mesh data efficiently
bool Mesh::loadOBJ(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error: Cannot open OBJ file " << filename << "\n";
        return false;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string type;
        iss >> type;

        if (type == "v") {
            Vertex v;
            iss >> v.x >> v.y >> v.z;
            vertices.push_back(v);
        } else if (type == "f") {
            Face f;
            std::string vertexData;
            while (iss >> vertexData) {
                std::istringstream vStream(vertexData);
                std::string vIdx;
                std::getline(vStream, vIdx);  
                f.vertexIndices.push_back(std::stoi(vIdx) - 1);
            }
            faces.push_back(f);
        }
    }

    file.close();
    CalcStats();
    CalcFaceNormal();

    return true;
}

bool Mesh::saveOBJ(const std::string& filename) const{
    std::ofstream file(filename, std::ios::out);
    if (!file) {
        std::cerr << "Error: Cannot open OBJ file for writing: " << filename << std::endl;
        return false;
    }

    file << "# Optimized OBJ file\n";
    file << "# Vertices: " << vertices.size() << "\n";
    file << "# Faces: " << faces.size() << "\n";

    // Write vertices efficiently
    for (const auto& v : vertices) {
        file << "v " << v.x << " " << v.y << " " << v.z << "\n";
    }

    // Write faces efficiently
    for (const auto& f : faces) {
        file << "f";
        for (const auto& idx : f.vertexIndices) {
            file << " " << (idx + 1); // OBJ indices are 1-based
        }
        file << "\n";
    }

    file.close();
    std::cout << "Successfully written OBJ file: " << filename << "\n";
    return true;
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