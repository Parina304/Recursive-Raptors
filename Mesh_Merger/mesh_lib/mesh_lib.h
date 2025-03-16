#ifndef _MESH_LIB_H
#define _MESH_LIB_H

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cstdint>
#include <string>

// Structure to store a vertex with 3D coordinates
struct Vertex {
    float x, y, z;
};

// Structure to store a face with vertex indices
struct Face {
    std::vector<uint32_t> vertexIndices;
};

// Mesh class to handle loading, conversion, and writing
class Mesh {
public:
    std::vector<Vertex> vertices;
    std::vector<Face> faces;
    std::vector<float> centroid;

public:
    bool loadOBJ(const std::string& filename);
    bool saveOBJ(const std::string& filename) const;
    bool savePLY(const std::string& filename) const;
    bool loadPLY(const std::string& filename);
    void printMeshStats() const;
};

#endif