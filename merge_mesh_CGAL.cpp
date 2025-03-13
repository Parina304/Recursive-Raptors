// CGAL includes
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/transform.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/IO/polygon_mesh_io.h>

// standard includes
#include <iostream>
namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;
typedef Kernel::Aff_transformation_3 Transformation;

int main() {
    Mesh sphere, humanoid;

    // Load meshes
    if (!CGAL::IO::read_polygon_mesh("/spherical_surface_smooth (1).obj", sphere) ||
        !CGAL::IO::read_polygon_mesh("/humanoid_robot_1.obj", humanoid)) {
        std::cerr << "Error reading mesh files!" << std::endl;
        return EXIT_FAILURE;
    }

    // Triangulate the humanoid if it has quadrilateral faces
    PMP::triangulate_faces(humanoid);

    // Translate the humanoid upwards based on user input
    double displacement = 0.0; 
    std::cout << "Please enter the desired displacement of humanoid relative to the sphere >> " ;
    std::cin >> displacement;
    Transformation translate(CGAL::TRANSLATION, Kernel::Vector_3(0, displacement, 0));
    PMP::transform(translate, humanoid);

    // Perform Boolean Union
    Mesh result;
    if (!PMP::corefine_and_compute_union(sphere, humanoid, result)) {
        std::cerr << "Boolean union failed!" << std::endl;
        return EXIT_FAILURE;
    }

    // Save the final merged mesh
    if (!CGAL::IO::write_polygon_mesh("merged.obj", result)) {
        std::cerr << "Error saving the merged mesh!" << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "Meshes successfully merged and saved as merged.obj!" << std::endl;
    return EXIT_SUCCESS;
}
