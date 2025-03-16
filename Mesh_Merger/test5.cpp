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

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> CGALMesh;
typedef Kernel::Aff_transformation_3 Transformation;


// Function to convert our Mesh class to CGAL::Surface_mesh
CGALMesh convertToCGALMesh(const Mesh& customMesh) {
    CGALMesh mesh;
    std::vector<CGALMesh::Vertex_index> v_indices;

    // Add vertices
    for (const auto& v : customMesh.vertices) {
        v_indices.push_back(mesh.add_vertex(Kernel::Point_3(v.x, v.y, v.z)));
    }

    // Add faces
    for (const auto& f : customMesh.faces) {
        std::vector<CGALMesh::Vertex_index> face_vertices;
        for (uint32_t index : f.vertexIndices) {
            face_vertices.push_back(v_indices[index]);
        }
        mesh.add_face(face_vertices);
    }
    
    return mesh;
}


// Function to convert CGAL::Surface_mesh back to our Mesh class
Mesh convertToNormalMesh(const CGALMesh& cgalMesh) {
    Mesh normalMesh;
    std::vector<Vertex> vertices;
    std::vector<Face> faces;
    std::map<CGALMesh::Vertex_index, uint32_t> vertexMap;
    
    // Extract vertices
    uint32_t index = 0;
    for (const auto v : cgalMesh.vertices()) {
        Kernel::Point_3 p = cgalMesh.point(v);
        Vertex vertex = { static_cast<float>(p.x()), static_cast<float>(p.y()), static_cast<float>(p.z()) };
        vertexMap[v] = index++;  // Map CGAL vertex index to new index
        vertices.push_back(vertex);
    }
    
    // Extract faces
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


int main() {
    Mesh customSphere, customHumanoid;
    
    // Load the meshes using test4_Nishi's Mesh class
    if (!customSphere.loadOBJ("D:/TAMU_2nd_Sem/Computing/Group_Assignments/Project2/spherical_surface_smooth (1).obj") ||
        !customHumanoid.loadOBJ("D:/TAMU_2nd_Sem/Computing/Group_Assignments/Project2/humanoid_robot_1.obj")) {
        std::cerr << "Error reading mesh files!" << std::endl;
        return EXIT_FAILURE;
    }

    // Convert loaded data to CGAL::Surface_mesh
    CGALMesh sphere = convertToCGALMesh(customSphere);
    CGALMesh humanoid = convertToCGALMesh(customHumanoid);

    // Triangulate the humanoid if it has quadrilateral faces
    PMP::triangulate_faces(humanoid);

    // Translate the humanoid upwards
    double displacement = 0.0;
    std::cout << "Please enter the desired displacement of humanoid relative to the sphere >> ";
    std::cin >> displacement;
    Transformation translate(CGAL::TRANSLATION, Kernel::Vector_3(0, displacement, 0));
    PMP::transform(translate, humanoid);

    // Perform Boolean Union
    CGALMesh resultCGAL;
    if (!PMP::corefine_and_compute_union(sphere, humanoid, resultCGAL)) {
        std::cerr << "Boolean union failed!" << std::endl;
        return EXIT_FAILURE;
    }

    Mesh Result = convertToNormalMesh(resultCGAL);
    // Save the final merged mesh
    if (!CGAL::IO::write_polygon_mesh("merged.obj", resultCGAL)) {
        std::cerr << "Error saving the merged mesh!" << std::endl;
        return EXIT_FAILURE;
    }


    std::cout << "Meshes successfully merged and saved as merged.obj!" << std::endl;
    return EXIT_SUCCESS;
}