#ifndef _MESH_LIB_H
#define _MESH_LIB_H

#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <sstream>
#include <cstdint>
#include <string>
#include <cmath>

// Structure to store a vertex with 3D coordinates
struct Vertex {
    float x, y, z;
};

// Structure to store a face with vertex indices
struct Face {
    std::vector<uint32_t> vertexIndices;
    std::vector<std::tuple<double, double>> materialData; // Example: (materialProp, thickness)
};

// Mesh class to handle loading, conversion, and writing
class Mesh {
public:
    std::vector<Vertex> vertices;
    std::vector<Face> faces;
    Vertex centroid;

public:
    // mesh i/o functions
    bool loadOBJ(const std::string& filename);
    bool saveOBJ(const std::string& filename) const;
    bool savePLY(const std::string& filename) const;
    bool loadPLY(const std::string& filename);
    void printMeshStats() const;
    void calc_centroid();

    // mesh operations
    void translate(float x, float y, float z);
    void rotate(float rad, int axis);
    void scale(float scale_fac);
    void to_origin();

};

#endif