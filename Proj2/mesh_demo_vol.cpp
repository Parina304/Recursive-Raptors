#include "mesh_lib.h"

// Main function
int main() {
    std::string fname_foam_in_points = "../assets/sphere/constant/polyMesh/points";
    std::string fname_foam_in_faces = "../assets/sphere/constant/polyMesh/faces";  // Use "faces" if OpenFOAM stores volumes there
    std::string fname_vtk_in = "../assets/robot/robot_2.vtk";  // Use "faces" if OpenFOAM stores volumes there
    std::string fname_out_sphere_vtk = "out_sphere_vtk.vtk";
    std::string fname_out_sphere_obj = "out_sphere_obj.obj";
    std::string fname_out_robot_vtk = "out_robot_vtk.vtk";
    std::string fname_out_robot_obj = "out_robot_obj.obj";
    VolMesh sphere, robot;

    sphere.parsePoints(fname_foam_in_points);
    sphere.parseCells(fname_foam_in_faces);

    if (sphere.vertices.empty() || sphere.faces.empty()) {
        std::cerr << "Error: Failed to parse input files.\n";
        return -1;
    }
    robot.parseVTK(fname_vtk_in);
    sphere.saveOBJ(fname_out_sphere_obj);
    sphere.writeVTK(fname_out_sphere_vtk);
    robot.saveOBJ(fname_out_robot_obj);
    robot.writeVTK(fname_out_robot_vtk);


    std::cout << "VTK volume mesh generation complete." << std::endl;
    return 0;
}