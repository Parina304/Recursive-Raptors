// Include necessary CGAL headers for mesh processing
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/transform.h>
#include <CGAL/Aff_transformation_3.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdint>
#include "mesh_lib.h"
#include <algorithm> 
#include <cctype>

namespace PMP = CGAL::Polygon_mesh_processing;

// Define kernel and mesh types for CGAL
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> CGALMesh;
typedef Kernel::Aff_transformation_3 Transformation;

// Function to convert custom Mesh class to CGAL::Surface_mesh
CGALMesh convertToCGALMesh(const Mesh& customMesh) {
    CGALMesh mesh;
    std::vector<CGALMesh::Vertex_index> v_indices;

    // Add vertices from the custom mesh
    for (const auto& v : customMesh.vertices) {
        v_indices.push_back(mesh.add_vertex(Kernel::Point_3(v.x, v.y, v.z)));
    }

    // Add faces using vertex indices
    for (const auto& f : customMesh.faces) {
        std::vector<CGALMesh::Vertex_index> face_vertices;
        for (uint32_t index : f.vertexIndices) {
            face_vertices.push_back(v_indices[index]);
        }
        mesh.add_face(face_vertices);
    }
    
    return mesh;
}

// Function to convert CGAL::Surface_mesh back to custom Mesh class
Mesh convertToNormalMesh(const CGALMesh& cgalMesh) {
    Mesh normalMesh;
    std::vector<Vertex> vertices;
    std::vector<Face> faces;
    std::map<CGALMesh::Vertex_index, uint32_t> vertexMap;
    
    // Extract vertices and create a mapping
    uint32_t index = 0;
    for (const auto v : cgalMesh.vertices()) {
        Kernel::Point_3 p = cgalMesh.point(v);
        Vertex vertex = { static_cast<float>(p.x()), static_cast<float>(p.y()), static_cast<float>(p.z()) };
        vertexMap[v] = index++;  // Map CGAL vertex index to new index
        vertices.push_back(vertex);
    }
    
    // Extract faces using the mapped indices
    for (const auto f : cgalMesh.faces()) {
        Face face;
        for (const auto v : CGAL::vertices_around_face(cgalMesh.halfedge(f), cgalMesh)) {
            face.vertexIndices.push_back(vertexMap[v]);
        }
        faces.push_back(face);
    }
    
    normalMesh.vertices = std::move(vertices);
    normalMesh.faces = std::move(faces);
    return normalMesh;
}

// Case-insensitive string comparison function
bool caseInsensitiveCompare(const std::string& a, const std::string& b) {
    std::string a_lower = a, b_lower = b;
    std::transform(a_lower.begin(), a_lower.end(), a_lower.begin(), ::tolower);
    std::transform(b_lower.begin(), b_lower.end(), b_lower.begin(), ::tolower);
    return a_lower == b_lower;
}

int main() {
    // Declare custom mesh objects
    Mesh customSphere, customHumanoid;
    
    // Initialize file paths
    std::string customHumanoidPath = "";
    std::string customSpherePath = "";
    
    // Get user input for mesh file paths
    std::cout << "Please enter the path of mesh file 1 i.e., sphere (type default for default) >> ";
    std::cin >> customSpherePath;
    std::cout << "Please enter the path of mesh file 2 i.e., humanoid (type default for default) >> ";
    std::cin >> customHumanoidPath;

    // Set default paths if needed
    if(caseInsensitiveCompare(customSpherePath,"default")) {
        customSpherePath = "../assets/obj/spherical_surface_smooth.obj";
    }
    if(caseInsensitiveCompare(customHumanoidPath,"default")) {
        customHumanoidPath = "../assets/obj/humanoid_robot.obj";
    }

    // Load the meshes from files
    if (!customSphere.loadOBJ(customSpherePath) || !customHumanoid.loadOBJ(customHumanoidPath)) {
        std::cerr << "Error reading mesh files!" << std::endl;
        return EXIT_FAILURE;
    }

    // Convert custom meshes to CGAL::Surface_mesh
    CGALMesh sphere = convertToCGALMesh(customSphere);
    CGALMesh humanoid = convertToCGALMesh(customHumanoid);

    // Ensure the humanoid mesh is triangulated
    PMP::triangulate_faces(humanoid);

    // Apply user-defined translation to the humanoid mesh
    double displacement = 0.0;
    std::cout << "\nPlease enter the desired displacement of humanoid relative to the sphere >> ";
    std::cin >> displacement;
    Transformation translate(CGAL::TRANSLATION, Kernel::Vector_3(0, displacement, 0));
    PMP::transform(translate, humanoid);

    // Perform Boolean union operation between the meshes
    CGALMesh resultCGAL;
    if (!PMP::corefine_and_compute_union(sphere, humanoid, resultCGAL)) {
        std::cerr << "Boolean union failed!" << std::endl;
        return EXIT_FAILURE;
    }
    
    // Convert the merged CGAL mesh back to the custom Mesh class
    Mesh Result = convertToNormalMesh(resultCGAL);

    // Get output filename from the user
    std::string outFileName = "";
    std::cout << "\nPlease enter the output file name ([name].obj) >> ";
    std::cin >> outFileName;
    
    // Save the merged mesh to a file
    if (!CGAL::IO::write_polygon_mesh(outFileName, resultCGAL)) {
        std::cerr << "Error saving the merged mesh!" << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "Meshes successfully merged and saved as " << outFileName << "!" << std::endl;
    return EXIT_SUCCESS;
}
