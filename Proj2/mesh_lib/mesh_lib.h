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
#include <regex>
#include <limits>

// Structure to store a vertex with 3D coordinates
struct Vertex {
    float x, y, z;
    Vertex(float x = 0, float y = 0, float z = 0){
        this->x = x;
        this->y = y;
        this->z = z;
    }
    Vertex operator + (const Vertex& v){
        return Vertex(this->x + v.x, this->y + v.y, this->z + v.z);
    }
    Vertex operator - (const Vertex& v){
        return Vertex(this->x - v.x, this->y - v.y, this->z - v.z);
    }
    Vertex operator / (float f){
        return Vertex(this->x / f, this->y / f, this->z / f);
    }
    float Norm(){
        return sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
    }
};

// struct Color{
//     float r, g, b;
//     Color(float r = 0, float g = 0, float b = 0){
//         this->r = r;
//         this->g = g;
//         this->b = b;
//     }
// };


// Structure to store a face with vertex indices
struct Face {
    std::vector<uint32_t> vertexIndices;
    Vertex centroid;
    float r, g, b;
};

// Mesh class to handle loading, conversion, and writing
class Mesh {
public:
    std::vector<Vertex> vertices;
    std::vector<Face> faces;
    std::vector<Vertex> face_normals;
    Vertex centroid;
    float x_min, y_min, z_min;
    float x_max, y_max, z_max;

public:
    Mesh();

    // mesh i/o functions
    bool loadOBJ(const std::string& filename);
    bool saveOBJ(const std::string& filename) const;
    bool savePLY(const std::string& filename) const;
    bool loadPLY(const std::string& filename);
    void printMeshStats() const;
    

    // mesh operations
    void translate(float x, float y, float z);
    void rotate(float rad, int axis);
    void scale(float scale_fac);
    void to_origin();
    void CalcStats();
    void CalcFaceNormal();
    void CalcFaceCentroid();

    Vertex operator - (Vertex &v1);
    Vertex operator + (Vertex &v1);

// private:
//     Vertex VectorCrossProduct(Vertex vec1, Vertex vec2);
//     float CalcNorm(Vertex v);
};

class VolMesh : public Mesh {
    
public:
    bool parsePoints(const std::string& filename);
    bool parseCells(const std::string& filename);
    bool parseVTK(const std::string& filename);
    bool writeVTK(const std::string& filename) const;
};

Vertex VectorCrossProduct(Vertex v1, Vertex v2);

#endif