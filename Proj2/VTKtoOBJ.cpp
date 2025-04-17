#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <chrono>
#include "mesh_lib.h" //Vertex and Face from here

//replace all Point and Cell references 

using namespace std;

// Function to parse VTK file (UNSTRUCTURED_GRID with hexahedral elements)
void parseVTK(const string& filename, vector<Vertex>& points, vector<Face>& faces) {
    ifstream file(filename);
    if (!file) {
        cerr << "Error opening VTK file: " << filename << endl;
        return;
    }
    string line;
    bool readingPoints = false, readingCells = false;
    int numPoints = 0, numCells = 0, cellDataSize = 0;

    while (getline(file, line)) {
        istringstream ss(line);
        string token;
        ss >> token;

        //reads in the points 
        if (token == "POINTS") { 
            ss >> numPoints;
            string dataType;
            ss >> dataType; // Read "float" or "double"

            points.resize(numPoints);
            for (int i = 0; i < numPoints; ++i) {
                file >> points[i].x >> points[i].y >> points[i].z;
            }
        } 
        //reads in the cells
        else if (token == "CELLS") { 
            ss >> numCells >> cellDataSize;
            faces.resize(numCells);

            for (int i = 0; i < numCells; ++i) {
                int numVertices;
                file >> numVertices; // First number is the vertex count
                faces[i].vertexIndices.resize(numVertices);
                for (int j = 0; j < numVertices; ++j) {
                    file >> faces[i].vertexIndices[j];
                }
            }
        }
    }
    file.close();
}

// Function to write the OBJ file with hexahedral faces
void writeOBJ(const string& filename, const vector<Vertex>& points, const vector<Face>& faces) {
    ofstream objFile(filename);
    if (!objFile) {
        cerr << "Error opening OBJ file for writing: " << filename << endl;
        return;
    }

    // Write vertex positions
    for (const auto& p : points) {
        objFile << "v " << p.x << " " << p.y << " " << p.z << "\n";
    }

    // Write hexahedral faces as quads
    for (const auto& face : faces) {
        if (face.vertexIndices.size() == 8) { // Ensure it's a hexahedron
            objFile << "f " << (face.vertexIndices[0] + 1) << " " << (face.vertexIndices[1] + 1) << " " << (face.vertexIndices[2] + 1) << " " << (face.vertexIndices[3] + 1) << "\n"; // Bottom Face
            objFile << "f " << (face.vertexIndices[4] + 1) << " " << (face.vertexIndices[5] + 1) << " " << (face.vertexIndices[6] + 1) << " " << (face.vertexIndices[7] + 1) << "\n"; // Top Face
            objFile << "f " << (face.vertexIndices[0] + 1) << " " << (face.vertexIndices[1] + 1) << " " << (face.vertexIndices[5] + 1) << " " << (face.vertexIndices[4] + 1) << "\n"; // Side Face 1
            objFile << "f " << (face.vertexIndices[1] + 1) << " " << (face.vertexIndices[2] + 1) << " " << (face.vertexIndices[6] + 1) << " " << (face.vertexIndices[5] + 1) << "\n"; // Side Face 2
            objFile << "f " << (face.vertexIndices[2] + 1) << " " << (face.vertexIndices[3] + 1) << " " << (face.vertexIndices[7] + 1) << " " << (face.vertexIndices[6] + 1) << "\n"; // Side Face 3
            objFile << "f " << (face.vertexIndices[3] + 1) << " " << (face.vertexIndices[0] + 1) << " " << (face.vertexIndices[4] + 1) << " " << (face.vertexIndices[7] + 1) << "\n"; // Side Face 4
        }
    }

    objFile.close();
    cout << "OBJ file saved: " << filename << endl;
}

// Main function
int main() {
    // record start time
    auto start = std::chrono::high_resolution_clock::now();

    string vtkFile = "../assets/robot/robot_2.vtk";  // Change this to your VTK file
    string objFile = "output2.obj";

    vector<Vertex> points;
    vector<Face> faces;

    parseVTK(vtkFile, points, faces);

    if (points.empty() || faces.empty()) {
        cerr << "Error: Failed to parse the VTK file.\n";
        return 1;
    }

    writeOBJ(objFile, points, faces);

    cout << "Conversion to OBJ complete. Open " << objFile << " in a 3D viewer." << endl;

    // record end time
    auto end = std::chrono::high_resolution_clock::now();
    // Calculate the duration
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    // Output the duration
    std::cout << "Function execution time: " << duration.count() << " microseconds" << std::endl;

    return 0;
}
