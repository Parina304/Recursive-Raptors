#ifndef MESH_LIB_H
#define MESH_LIB_H

#include <vector>
#include <string>
#include <glm/glm.hpp>

// Define material types for faces
enum MaterialType {
    INSULATION,
    GLUE,
    STEEL,
    CFRP
};

// Face structure: holds indices of the vertices and associated material & temperature
struct Face {
    int v[3];            // indices into the vertex list (for the original mesh vertices)
    MaterialType mat;    // material type of this face
    float temp;          // representative temperature of this face (for coloring)
};

class Mesh {
public:
    std::vector<glm::vec3> vertices;      // unique vertices from OBJ
    std::vector<Face> faces;             // faces with indices into vertices
    std::vector<glm::vec3> flatVertices; // flattened vertex array (each face's vertices duplicated) for GL drawing
    std::vector<glm::vec3> flatColors;   // parallel array of colors for each vertex in flatVertices

    // Bounding values (for camera and cut plane)
    float minX, maxX;
    float minY, maxY;
    float minZ, maxZ;
    // Temperature range (for legend in temperature mode)
    float minTemp, maxTemp;

    Mesh() : minX(0), maxX(0), minY(0), maxY(0), minZ(0), maxZ(0), minTemp(0), maxTemp(0) {}

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
