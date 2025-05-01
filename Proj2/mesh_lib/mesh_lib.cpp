#include "mesh_lib.h"

// Reads an OBJ file and stores its mesh data efficiently
bool Mesh::loadOBJ(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error: Cannot open OBJ file " << filename << "\n";
        return false;
    }

    std::string line;
    while (std::getline(file, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;
        if (line.rfind("v ", 0) == 0) {
            try{
                // Vertex line
                std::istringstream iss(line.substr(2));
                glm::vec3 v;
                if (!(iss >> v.x >> v.y >> v.z)){
                    throw std::runtime_error("Invalid vertex data");
                }
                tempVerts.push_back(v);
                // update bounds
                if (v.x < minX) minX = v.x;
                if (v.x > maxX) maxX = v.x;
                if (v.y < minY) minY = v.y;
                if (v.y > maxY) maxY = v.y;
                if (v.z < minZ) minZ = v.z;
                if (v.z > maxZ) maxZ = v.z;
            } catch (const std::exception& e) {
                std::cerr << "WARNING: Skipping invalid vertex line: " << line << " (" << e.what() << ")" << std::endl;
                continue;
            } 
        }
        else if (line.rfind("usemtl ", 0) == 0) {
            // Material usage for subsequent faces
            currentMatName = trim(line.substr(7));
            // Convert material name to enum
            currentMat = materialTypeFromName(currentMatName);
        }
        else if (line.rfind("f ", 0) == 0) {
            try{
                // Face line (triangle assumed). OBJ faces are 1-indexed.
                // We also handle faces with texture/normal: "f v/t/n v/t/n v/t/n"
                std::istringstream iss(line.substr(2));
                std::vector<std::string> vertex_tokens;
                std:: string token;
                while (iss >> token){
                    vertex_tokens.push_back(token);
                }

                // Parse each vertex index group up to '/', ignoring texture/normal indices.
                auto parseIndex = [](const std::string& token) {
                    // token might be "int" or "int/int/int" etc.
                    size_t slashPos = token.find('/');
                    int idx = 0;
                    try{
                        if (slashPos == std::string::npos) {
                            // no slash
                            idx = std::stoi(token);
                        }
                        else {
                            idx = std::stoi(token.substr(0, slashPos));
                        }
                    return idx;
                    } catch (const std::invalid_argument& e){
                        throw std::invalid_argument("Non-numverical vertex index");
                    } catch (const std::out_of_range& e){
                        throw std::out_of_range("Vertex index out of range.");
                    }
                };

                // Parse vertex tokens to indices
                std::vector<int> vertex_indices;
                for (const auto& vertex_token: vertex_tokens){
                    vertex_indices.push_back(parseIndex(vertex_token));
                }
                if (vertex_indices.size() < 3){
                    throw std::runtime_error("Face with fewer than 3 vertices.");
                }
                for (size_t i = 1; i < vertex_indices.size() - 1; ++i){
                    Face face;
                    face.v[0] = vertex_indices[0] - 1;
                    face.v[1] = vertex_indices[i] - 1;
                    face.v[2] = vertex_indices[i + 1] - 1;
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
            } catch (const std::exception& e){
                std::cerr << "WARNING: Skipping invalid face line: " << line << " (" << e.what() << ")" << std::endl;
                continue;
            }
        }
    }

    file.close();
    calc_centroid();
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
    calc_centroid();

    return true;
}

// Prints basic mesh statistics
void Mesh::printMeshStats() const {
    std::cout << "Mesh Statistics:\n";
    std::cout << "  - Vertices: " << vertices.size() << "\n";
    std::cout << "  - Faces: " << faces.size() << "\n";
    std::cout << "  - Centroid: x: " << centroid.x << " y: " << centroid.y << " z: " << centroid.z << "\n";
}

void Mesh::calc_centroid(){
    float cx = 0, cy = 0, cz = 0;
    for(const auto& v: vertices){
        cx += v.x;
        cy += v.y;
        cz += v.z;
    }
    centroid.x = cx / vertices.size();
    centroid.y = cy / vertices.size();
    centroid.z = cz / vertices.size();
}