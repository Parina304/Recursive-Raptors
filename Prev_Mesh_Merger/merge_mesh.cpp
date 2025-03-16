#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <windows.h>

struct Vertex {
    double x, y, z;
};

struct Face {
    std::vector<int> vertexIndices;
};

// Read an OBJ file
bool readOBJ(const std::string& filename, std::vector<Vertex>& vertices, 
             std::vector<Face>& faces, bool translate) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error: Cannot open file " << filename << "!\n";
        return false;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string type;
        iss >> type;

        if (type == "v") {  // Vertex
            Vertex v;
            iss >> v.x >> v.y >> v.z;
            if (translate)
                v.y += 0.7;
            vertices.push_back(v);
        } else if (type == "f") {  // Face
            Face f;
            std::string vertexData;
            while (iss >> vertexData) {
                std::istringstream vStream(vertexData);
                std::string vIdx;
                std::getline(vStream, vIdx, '/');
                f.vertexIndices.push_back(std::stoi(vIdx));
            }
            faces.push_back(f);
        }
    }
    file.close();
    return true;
}

// Merge two OBJ files and write the result
void mergeOBJ(const std::string& file1, const std::string& file2, const std::string& outputFile) {
    std::vector<Vertex> vertices1, vertices2;
    std::vector<Face> faces1, faces2;

    if (!readOBJ(file1, vertices1, faces1, false) || 
        !readOBJ(file2, vertices2, faces2, true)) {
        return;
    }

    // Compute offsets
    int vertexOffset = vertices1.size();

    // Offset face indices in the second file
    for (auto& f : faces2) {
        for (auto& v : f.vertexIndices) v += vertexOffset;
    }

    // Write merged file
    std::ofstream outFile(outputFile);
    if (!outFile) {
        std::cerr << "Error: Cannot write to file " << outputFile << "!\n";
        return;
    }

    for (const auto& v : vertices1) outFile << "v " << v.x << " " << v.y << " " << v.z << "\n";
    for (const auto& v : vertices2) outFile << "v " << v.x << " " << v.y << " " << v.z << "\n";

    for (const auto& f : faces1) {
        outFile << "f";
        for (size_t i = 0; i < f.vertexIndices.size(); ++i) {
            outFile << " " << f.vertexIndices[i];
        }
        outFile << "\n";
    }

    for (const auto& f : faces2) {
        outFile << "f";
        for (size_t i = 0; i < f.vertexIndices.size(); ++i) {
            outFile << " " << f.vertexIndices[i];
        }
        outFile << "\n";
    }

    outFile.close();
    std::cout << "Merged OBJ saved as: " << outputFile << "\n";
}

void ShowWarningPopup(const char* title, const char* message) {
    MessageBoxA(NULL, message, title, MB_OK | MB_ICONINFORMATION);
}

int main() {
    std::string file1 = "spherical_surface_smooth (1).obj";  // First modified mesh
    std::string file2 = "humanoid_robot_1.obj";      // Original second mesh
    std::string outputFile = "merged_humanoid.obj";    // Output merged mesh

    mergeOBJ(file1, file2, outputFile);
    ShowWarningPopup("Done!", "File Written Successfully!");
    return 0;
}
