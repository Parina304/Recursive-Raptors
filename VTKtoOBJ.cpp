#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include "mesh_lib.h"

//replace all Point and Cell references 

using namespace std;

// Struct for 3D point coordinates
struct Point {
    double x, y, z;
};

// Struct for storing hexahedral cell connectivity
using Cell = vector<int>;

// Function to parse VTK file (UNSTRUCTURED_GRID with hexahedral elements)
void parseVTK(const string& filename, vector<Point>& points, vector<Cell>& cells) {
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

        if (token == "POINTS") { 
            ss >> numPoints;
            string dataType;
            ss >> dataType; // Read "float" or "double"

            points.resize(numPoints);
            for (int i = 0; i < numPoints; ++i) {
                file >> points[i].x >> points[i].y >> points[i].z;
            }
        } 
        else if (token == "CELLS") { 
            ss >> numCells >> cellDataSize;
            cells.resize(numCells);

            for (int i = 0; i < numCells; ++i) {
                int numVertices;
                file >> numVertices; // First number is the vertex count
                cells[i].resize(numVertices);
                for (int j = 0; j < numVertices; ++j) {
                    file >> cells[i][j];
                }
            }
        }
    }
    file.close();
}

// Function to write the OBJ file with hexahedral faces
void writeOBJ(const string& filename, const vector<Point>& points, const vector<Cell>& cells) {
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
    for (const auto& cell : cells) {
        if (cell.size() == 8) { // Ensure it's a hexahedron
            objFile << "f " << (cell[0] + 1) << " " << (cell[1] + 1) << " " << (cell[2] + 1) << " " << (cell[3] + 1) << "\n"; // Bottom Face
            objFile << "f " << (cell[4] + 1) << " " << (cell[5] + 1) << " " << (cell[6] + 1) << " " << (cell[7] + 1) << "\n"; // Top Face
            objFile << "f " << (cell[0] + 1) << " " << (cell[1] + 1) << " " << (cell[5] + 1) << " " << (cell[4] + 1) << "\n"; // Side Face 1
            objFile << "f " << (cell[1] + 1) << " " << (cell[2] + 1) << " " << (cell[6] + 1) << " " << (cell[5] + 1) << "\n"; // Side Face 2
            objFile << "f " << (cell[2] + 1) << " " << (cell[3] + 1) << " " << (cell[7] + 1) << " " << (cell[6] + 1) << "\n"; // Side Face 3
            objFile << "f " << (cell[3] + 1) << " " << (cell[0] + 1) << " " << (cell[4] + 1) << " " << (cell[7] + 1) << "\n"; // Side Face 4
        }
    }

    objFile.close();
    cout << "OBJ file saved: " << filename << endl;
}

// Main function
int main() {
    string vtkFile = "Robot_2.vtk";  // Change this to your VTK file
    string objFile = "output2.obj";

    vector<Point> points;
    vector<Cell> cells;

    parseVTK(vtkFile, points, cells);

    if (points.empty() || cells.empty()) {
        cerr << "Error: Failed to parse the VTK file.\n";
        return 1;
    }

    writeOBJ(objFile, points, cells);

    cout << "Conversion to OBJ complete. Open " << objFile << " in a 3D viewer." << endl;
    return 0;
}
