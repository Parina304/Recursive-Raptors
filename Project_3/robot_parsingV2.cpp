// Copyright 2025 Recursive Raptors
/**
 * @file robot_thermal_mesh_processor.cpp
 * @brief Processes robot mesh files for thermal simulation analysis
 * 
 * This program reads VTK mesh files of a humanoid robot, assigns material properties
 * based on position-dependent thickness profiles, and outputs an OBJ file with
 * material information for subsequent thermal simulations.
 */

 #include <cmath>
 #include <chrono>
 #include <fstream>
 #include <iostream>
 #include <limits>
 #include <sstream>
 #include <string>
 #include <unordered_map>
 #include <vector>
 
 #include "mesh_lib.h"  // Vertex and Face
 
 namespace {
 
 constexpr double kRobotLength = 2.5;  // meters
 
 // Material property structure
 struct MaterialProperties {
   double thermal_conductivity;
   double specific_heat;
   double glass_transition_temp;
   double density;
 };
 
 // Material map
 const std::unordered_map<std::string, MaterialProperties> kMaterialDB = {
     {"Steel", {50.2, 490, 600, 7850}},
     {"Glue", {0.2, 1200, 373, 1100}},
     {"Carbon", {8.0, 710, 400, 1600}},
     {"ThermalProtection", {0.05, 1300, 550, 350}}};
 
 // Thickness profiles as a function of normalized position l/L
 double GetSteelThickness(double l_norm) {
   return 0.01 + 0.005 * std::sin(M_PI * l_norm);
 }
 
 double GetGlueThickness(double l_norm) {
   return 0.002 + 0.001 * std::cos(2 * M_PI * l_norm);
 }
 
 double GetCarbonThickness(double l_norm) {
   return 0.015 + 0.005 * std::fabs(std::sin(2 * M_PI * l_norm));
 }
 
 // Assigns material stack to a face based on radial distance and axial location
 std::vector<std::string> AssignMaterialStack(double r, double l_norm) {
   double steel = GetSteelThickness(l_norm);
   double glue = GetGlueThickness(l_norm);
   double carbon = GetCarbonThickness(l_norm);
 
   std::vector<std::string> stack;
 
   if (r <= steel) stack.push_back("Steel");
   if (r > steel && r <= steel + glue) stack.push_back("Glue");
   if (r > steel + glue && r <= steel + glue + carbon) stack.push_back("Carbon");
 
   if (l_norm < 0.1 || l_norm > 0.9) stack.push_back("ThermalProtection");
 
   return stack;
 }
 
 // Parses VTK mesh file and assigns material stacks to each face
 void ParseVTK(const std::string& filename, std::vector<Vertex>& points,
               std::vector<Face>& faces, std::vector<double>& normalizedZ,
               std::vector<std::vector<std::string>>& faceMaterialStacks) {
   std::ifstream file(filename);
   if (!file) {
     std::cerr << "Error opening VTK file: " << filename << std::endl;
     return;
   }
 
   std::string line;
   int num_points = 0, num_cells = 0, cell_data_size = 0;
 
   double z_min = std::numeric_limits<double>::max();
   double z_max = std::numeric_limits<double>::lowest();
 
   while (std::getline(file, line)) {
     std::istringstream ss(line);
     std::string token;
     ss >> token;
 
     if (token == "POINTS") {
       ss >> num_points;
       std::string data_type;
       ss >> data_type;
 
       points.resize(num_points);
       for (int i = 0; i < num_points; ++i) {
         std::string point_line;
         if (!std::getline(file, point_line)) {
           std::cerr << "Unexpected end of file while reading points." << std::endl;
           return;
         }
 
         std::istringstream linestream(point_line);
         if (!(linestream >> points[i].x >> points[i].y >> points[i].z)) {
           std::cerr << "Error parsing point " << i << ": '" << point_line << "'" << std::endl;
           return;
         }
 
         z_min = std::min(z_min, points[i].z);
         z_max = std::max(z_max, points[i].z);
       }
 
       if (points.size() != static_cast<size_t>(num_points)) {
         std::cerr << "Mismatch in point count." << std::endl;
         return;
       }
     } else if (token == "CELLS") {
       ss >> num_cells >> cell_data_size;
       faces.resize(num_cells);
       normalizedZ.resize(num_cells);
       faceMaterialStacks.resize(num_cells);
 
       for (int i = 0; i < num_cells; ++i) {
         int num_vertices;
         if (!(file >> num_vertices)) {
           std::cerr << "Error reading number of vertices for cell " << i << std::endl;
           return;
         }
 
         faces[i].vertexIndices.resize(num_vertices);
         double z_avg = 0.0, x_avg = 0.0, y_avg = 0.0;
 
         for (int j = 0; j < num_vertices; ++j) {
           if (!(file >> faces[i].vertexIndices[j])) {
             std::cerr << "Error reading vertex index for face " << i << std::endl;
             return;
           }
 
           if (faces[i].vertexIndices[j] >= static_cast<int>(points.size())) {
             std::cerr << "Vertex index out of range." << std::endl;
             return;
           }
 
           const Vertex& v = points[faces[i].vertexIndices[j]];
           x_avg += v.x;
           y_avg += v.y;
           z_avg += v.z;
         }
 
         x_avg /= num_vertices;
         y_avg /= num_vertices;
         z_avg /= num_vertices;
 
         if (kRobotLength <= 0.0) {
           std::cerr << "Invalid robot length." << std::endl;
           return;
         }
 
         double l_norm = (z_avg - z_min) / kRobotLength;
         double r = std::sqrt(x_avg * x_avg + y_avg * y_avg);
 
         normalizedZ[i] = l_norm;
         faceMaterialStacks[i] = AssignMaterialStack(r, l_norm);
       }
 
       if (faces.size() != static_cast<size_t>(num_cells)) {
         std::cerr << "Mismatch in cell count." << std::endl;
         return;
       }
     }
   }
 
   file.close();
 }
 
 // Writes mesh data to an OBJ file
 void WriteOBJ(const std::string& filename, const std::vector<Vertex>& points,
               const std::vector<Face>& faces) {
   std::ofstream obj_file(filename);
   if (!obj_file) {
     std::cerr << "Error opening OBJ file for writing: " << filename << std::endl;
     return;
   }
 
   for (const auto& p : points) {
     obj_file << "v " << p.x << " " << p.y << " " << p.z << "\n";
   }
 
   for (const auto& face : faces) {
     if (face.vertexIndices.size() == 8) {
       obj_file << "f " << (face.vertexIndices[0] + 1) << " " << (face.vertexIndices[1] + 1)
                << " " << (face.vertexIndices[2] + 1) << " " << (face.vertexIndices[3] + 1) << "\n";
       obj_file << "f " << (face.vertexIndices[4] + 1) << " " << (face.vertexIndices[5] + 1)
                << " " << (face.vertexIndices[6] + 1) << " " << (face.vertexIndices[7] + 1) << "\n";
       obj_file << "f " << (face.vertexIndices[0] + 1) << " " << (face.vertexIndices[1] + 1)
                << " " << (face.vertexIndices[5] + 1) << " " << (face.vertexIndices[4] + 1) << "\n";
       obj_file << "f " << (face.vertexIndices[1] + 1) << " " << (face.vertexIndices[2] + 1)
                << " " << (face.vertexIndices[6] + 1) << " " << (face.vertexIndices[5] + 1) << "\n";
       obj_file << "f " << (face.vertexIndices[2] + 1) << " " << (face.vertexIndices[3] + 1)
                << " " << (face.vertexIndices[7] + 1) << " " << (face.vertexIndices[6] + 1) << "\n";
       obj_file << "f " << (face.vertexIndices[3] + 1) << " " << (face.vertexIndices[0] + 1)
                << " " << (face.vertexIndices[4] + 1) << " " << (face.vertexIndices[7] + 1) << "\n";
     }
   }
 
   obj_file.close();
   std::cout << "OBJ file saved: " << filename << std::endl;
 }
 
 }  // namespace
 
 int main() {
   auto start = std::chrono::high_resolution_clock::now();
 
   const std::string vtk_file = "../assets/robot/robot_2.vtk";
   const std::string obj_file = "output2.obj";
 
   std::vector<Vertex> points;
   std::vector<Face> faces;
   std::vector<double> normalized_z;
   std::vector<std::vector<std::string>> face_material_stacks;
 
   ParseVTK(vtk_file, points, faces, normalized_z, face_material_stacks);
 
   if (points.empty() || faces.empty()) {
     std::cerr << "Error: Failed to parse the VTK file." << std::endl;
     return 1;
   }
 
   WriteOBJ(obj_file, points, faces);
 
   for (size_t i = 0; i < face_material_stacks.size(); ++i) {
     std::cout << "Face " << i << " l/L = " << normalized_z[i] << " -> ";
     for (const auto& mat : face_material_stacks[i]) {
       std::cout << mat << " ";
     }
     std::cout << std::endl;
   }
 
   auto end = std::chrono::high_resolution_clock::now();
   auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
   std::cout << "Execution time: " << duration.count() << " microseconds" << std::endl;
 
   return 0;
 } 