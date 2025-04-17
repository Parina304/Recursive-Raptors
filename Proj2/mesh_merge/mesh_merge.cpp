#include "mesh_merge.h"

namespace PMP = CGAL::Polygon_mesh_processing;

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

    normalMesh.calc_centroid();

    return normalMesh;
}

Mesh meshUnion(const Mesh& mesh1, const Mesh& mesh2){

    // Convert custom meshes to CGAL::Surface_mesh
    CGALMesh cmesh1 = convertToCGALMesh(mesh1);
    CGALMesh cmesh2 = convertToCGALMesh(mesh2);

    // Ensure the humanoid mesh is triangulated
    PMP::triangulate_faces(cmesh1);
    PMP::triangulate_faces(cmesh2);

    // Perform Boolean union operation between the meshes
    CGALMesh resultCGAL;
    if (!PMP::corefine_and_compute_union(cmesh1, cmesh2, resultCGAL)) {
        std::cerr << "Boolean union failed!" << std::endl;
        // return EXIT_FAILURE;
    }
    
    // Convert the merged CGAL mesh back to the custom Mesh class
    Mesh Result = convertToNormalMesh(resultCGAL);

    return Result;
}

// Case-insensitive string comparison function
bool caseInsensitiveCompare(const std::string& a, const std::string& b) {
    std::string a_lower = a, b_lower = b;
    std::transform(a_lower.begin(), a_lower.end(), a_lower.begin(), ::tolower);
    std::transform(b_lower.begin(), b_lower.end(), b_lower.begin(), ::tolower);
    return a_lower == b_lower;
}

std::string getInputWithDefault(const std::string& prompt, const std::string& defaultValue) {
    std::string input;
    std::cout << prompt << " [ no input defaults to " << defaultValue << "]: ";
    std::getline(std::cin, input);
    if (input.empty()) {
        return defaultValue;
    }
    return input;
}

