#include "mesh_lib.h"
#include "mesh_merge.h"

// Main function to test OBJ to PLY conversion and reading
int main(int argc, char* argv []) {
    Mesh obj1, obj2, mesh_union;
    std::string fname_in_1 = "../assets/obj/spherical_surface_smooth.obj";
    std::string fname_in_2 = "../assets/obj/humanoid_robot.obj";
    std::string fname_out = "out_union.obj";

    // Load the OBJ file
    std::cout << "This program accepts two .obj files, operates on the 2nd object, then outputs their union.\n";
    if (argc > 1){
        fname_in_1 = argv[1];
    } else {
        fname_in_1 = getInputWithDefault("Enter path to 1st OBJ file", fname_in_1);
    }
    if (!obj1.loadOBJ(fname_in_1)) {
        return -1;
    }
    if (argc > 2){
        fname_in_2 = argv[2];
    } else {
        fname_in_2 = getInputWithDefault("Enter path to 2nd OBJ file", fname_in_2);
    }
    if (!obj2.loadOBJ(fname_in_2)) {
        return -1;
    }

    // Repeatedly prompt the user for input to perform operations on obj2
    std::string userInput;
    while (true) {
        std::cout << "Enter operation for obj2 (h for help): ";
        std::getline(std::cin, userInput);

        std::istringstream iss(userInput);
        char command;
        iss >> command;

        if (command == 't') {
            float x, y, z;
            iss >> x >> y >> z;
            obj2.translate(x, y, z);
        } else if (command == 'r') {
            float angle;
            int axis;
            iss >> angle >> axis;
            obj2.rotate(angle, axis);
        } else if (command == 's') {
            float scaleFactor;
            iss >> scaleFactor;
            if(!scaleFactor){
                std::cout << "Warning: Can't scale by 0! Scaling aborted.\n";
            } else{
                obj2.scale(scaleFactor);
            }
        }else if (command == 'i'){
            obj2.printMeshStats();
        } else if (command == 'u') {
            break; // Exit the loop to perform the union operation
        } else if (command == 'h'){
            std::cout << "Available commands:\n";
            std::cout << "  t x y z  - Translate by x, y, z\n";
            std::cout << "  r a k    - Rotate by angle a (in radians) around axis k (0 for z, 1 for x, 2 for y)\n";
            std::cout << "  s c      - Scale by factor c\n";
            std::cout << "  i        - Print obj2 stats\n";
            std::cout << "  u        - Union the two objects\n";
            std::cout << "  h        - Display this help message\n";
        } else {
            std::cout << "Invalid operation. Please enter t x y z, r a k, s c, or u." << std::endl;
        }
    }
    if (argc > 3){
        fname_out = argv[3];
    } else {
        fname_out = getInputWithDefault("Enter path to output OBJ file", fname_out);
    }
    mesh_union = meshUnion(obj1, obj2);
    mesh_union.printMeshStats();
    if (!mesh_union.saveOBJ(fname_out)) {
        return -1;
    }
    // std::cout << ".obj saved at " << fname_out_head << "\n";
    // humanoid.translate(0, 4, 0);
    // feet_in = meshUnion(humanoid, sphere);
    // feet_in.printMeshStats();
    // if (!feet_in.saveOBJ(fname_out_feet)) {
    //     return -1;
    // }
    // std::cout << ".obj saved at " << fname_out_feet << "\n";
    return 0;
}