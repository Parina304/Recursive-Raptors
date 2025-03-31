#include "mesh_lib.h"
#include "mesh_merge.h"

// Main function to test OBJ to PLY conversion and reading
int main() {
    Mesh sphere, humanoid, head_out, feet_in;
    std::string fname_in_1 = "../assets/obj/humanoid_robot.obj";
    std::string fname_in_2 = "../assets/obj/spherical_surface_smooth.obj";
    std::string fname_out_head = "head_out.obj";
    std::string fname_out_feet = "feet_in.obj";


    // Load the OBJ file
    if (!humanoid.loadOBJ(fname_in_1)) {
        return -1;
    }

    if (!sphere.loadOBJ(fname_in_2)) {
        return -1;
    }
    // Print statistics and save as Binary PLY
    humanoid.printMeshStats();
    humanoid.to_origin();
    humanoid.translate(0, 1.4, 0);
    humanoid.printMeshStats();
    sphere.scale(2.5);
    // sphere.saveOBJ("spehre_out.obj");
    head_out = meshUnion(humanoid, sphere);
    head_out.printMeshStats();
    if (!head_out.saveOBJ(fname_out_head)) {
        return -1;
    }
    std::cout << ".obj saved at " << fname_out_head << "\n";
    humanoid.translate(0, 4, 0);
    feet_in = meshUnion(humanoid, sphere);
    feet_in.printMeshStats();
    if (!feet_in.saveOBJ(fname_out_feet)) {
        return -1;
    }
    std::cout << ".obj saved at " << fname_out_feet << "\n";
    return 0;
}