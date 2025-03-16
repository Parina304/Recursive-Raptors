#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <regex>

using namespace std;

struct Point {
    double x, y, z;
};

using Cell = vector<int>;

// Function to parse OpenFOAM points file
vector<Point> parsePoints(const string& filename) {
    ifstream file(filename);
    vector<Point> points;
    string line;

    if (!file) {
        cerr << "Error opening file: " << filename << endl;
        return points;
    }

    // Skip comments and header lines
    while (getline(file, line)) {
        if (!line.empty() && isdigit(line[0])) {
            break; // Found the number of points
        }
    }

    int numPoints = stoi(line);
    cout << "Detected " << numPoints << " points." << endl;

    // Skip the next line (usually "(" in OpenFOAM)
    getline(file, line);

    // Read point coordinates
    for (int i = 0; i < numPoints; ++i) {
        getline(file, line);
        line = regex_replace(line, regex("[()]"), ""); // Remove parentheses

        istringstream ss(line);
        double x, y, z;
        if (ss >> x >> y >> z) {
            points.push_back({x, y, z});
        } else {
            cerr << "Warning: Invalid point format at line " << i + 1 << ": " << line << endl;
        }
    }

    file.close();
    return points;
}

// Function to parse OpenFOAM faces/volume file
vector<Cell> parseCells(const string& filename) {
    ifstream file(filename);
    vector<Cell> cells;
    string line;

    if (!file) {
        cerr << "Error opening file: " << filename << endl;
        return cells;
    }

    // Skip comments and header lines
    while (getline(file, line)) {
        if (!line.empty() && isdigit(line[0])) {
            break; // Found the number of cells
        }
    }

    int numCells = stoi(line);
    cout << "Detected " << numCells << " cells." << endl;

    // Skip the next line (usually "(" in OpenFOAM)
    getline(file, line);

    // Read cell data
    for (int i = 0; i < numCells; ++i) {
        getline(file, line);

        // Match lines like "4(1087 1208 1209 1088)" and extract only the numbers inside parentheses
        smatch match;
        regex cellRegex(R"(\d+\(([^)]+)\))");  // Matches "4(...)" and captures contents inside "()"

        if (regex_search(line, match, cellRegex)) {
            string numbers = match[1];  // Extract the part inside parentheses
            istringstream ss(numbers);
            Cell cell;
            int vertexIndex;

            while (ss >> vertexIndex) {
                cell.push_back(vertexIndex - 1); // Convert to zero-based index
            }

            if (!cell.empty()) {
                cells.push_back(cell);
            } else {
                cerr << "Warning: Empty cell line: " << line << endl;
            }
        } else {
            cerr << "Warning: Skipping invalid cell line: " << line << endl;
        }
    }

    file.close();
    return cells;
}


// Function to write parsed volume mesh data into a VTK file
void writeVTK(const string& filename, const vector<Point>& points, const vector<Cell>& cells) {
    ofstream vtkFile(filename);
    if (!vtkFile) {
        cerr << "Error opening VTK file for writing: " << filename << endl;
        return;
    }

    vtkFile << "# vtk DataFile Version 2.0\n";
    vtkFile << "OpenFOAM Volume Mesh\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET POLYDATA\n";

    // Write points
    vtkFile << "POINTS " << points.size() << " float\n";
    for (const auto& p : points) {
        vtkFile << p.x << " " << p.y << " " << p.z << "\n";
    }

    // Count total indices (including cell size prefixes)
    int totalIndices = 0;
    for (const auto& cell : cells) {
        totalIndices += cell.size() + 1; // +1 for cell size prefix
    }

    // Write faces
    vtkFile << "POLYGONS " << cells.size() << " " << totalIndices << "\n";
    for (const auto& cell : cells) {
        vtkFile << cell.size();
        // Write the cell indices without an extra space at the end
        for (size_t i = 0; i < cell.size(); ++i) {
            vtkFile << " " << cell[i];
        }
        vtkFile << "\n";
    }


    vtkFile.close();
    cout << "VTK file saved: " << filename << endl;
}

// Main function
int main() {
    string pointsFile = "points";
    string cellsFile = "faces";  // Use "faces" if OpenFOAM stores volumes there
    string vtkFile = "output_volume.vtk";

    vector<Point> points = parsePoints(pointsFile);
    vector<Cell> cells = parseCells(cellsFile);

    if (points.empty() || cells.empty()) {
        cerr << "Error: Failed to parse input files.\n";
        return 1;
    }

    writeVTK(vtkFile, points, cells);

    cout << "VTK volume mesh generation complete. Open " << vtkFile << " in ParaView." << endl;
    return 0;
}
