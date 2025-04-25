#ifndef MESH_LIB_H
#define MESH_LIB_H

#include <vector>
#include <string>
#include <glm/glm.hpp>

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

// Define material types for faces
enum MaterialType {
    INSULATION,
    GLUE,
    STEEL,
    CFRP
};

// Face structure: holds indices of the vertices and associated material & temperature
struct Face {
    std::vector<uint32_t> vertexIndices;
    Vertex centroid;
    float r, g, b;

// Mesh class to handle loading, conversion, and writing
    int v[3];            // indices into the vertex list (for the original mesh vertices)
    MaterialType mat;    // material type of this face
    float temp;          // representative temperature of this face (for coloring)
};

class Mesh {
public:
    std::vector<Vertex> vertices;
    std::vector<Face> faces;
    std::vector<Vertex> face_normals;
    Vertex centroid;
    float x_min, y_min, z_min;
    float x_max, y_max, z_max;
    float minTemp, maxTemp;
    float minX, maxX;
    float minY, maxY;
    float minZ, maxZ;

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

    Vertex VectorCrossProduct(Vertex v1, Vertex v2);

    std::vector<glm::vec3> vertices;      // unique vertices from OBJ
    std::vector<Face> faces;             // faces with indices into vertices
    std::vector<glm::vec3> flatVertices; // flattened vertex array (each face's vertices duplicated) for GL drawing
    std::vector<glm::vec3> flatColors;   // parallel array of colors for each vertex in flatVertices

    // // Bounding values (for camera and cut plane)
    // float minX, maxX;
    // float minY, maxY;
    // float minZ, maxZ;
    // // Temperature range (for legend in temperature mode)
    // float minTemp, maxTemp;

    // Mesh() : minX(0), maxX(0), minY(0), maxY(0), minZ(0), maxZ(0), minTemp(0), maxTemp(0) {}

    // Load mesh from an OBJ file. Returns true on success.
    bool loadOBJ(const std::string& filePath);

    // After loading, prepare flattened buffers for rendering (populate flatVertices and flatColors)
    void createFlatArrays();

    // Utility: convert material name (from OBJ/MTL) to our MaterialType enum
    MaterialType materialTypeFromName(const std::string& name) const;

    // Get a color for a given material type (for Material view mode)
    glm::vec3 getMaterialColor(MaterialType mat) const;

    // Compute an RGB color (0-1 each) for a given temperature value using jet colormap (Temperature mode)
    glm::vec3 getColorForTemperature(float value) const;

    // Increase thermal protection layer thickness (simulate adding insulation)
    void increaseProtectionThickness();

    // “Solve” the thermal simulation: update face temperatures (this is a stub for actual simulation)
    void solveThermal();
};

#endif // MESH_LIB_H
