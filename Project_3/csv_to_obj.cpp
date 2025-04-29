// =============================================================================
// Title: Humanoid Robot OBJ with Insulation Layer Generator
//
// Description:
//   - This program reads a 2D humanoid robot mesh from an OBJ file.
//   - It reads a thickness profile (height vs insulation thickness) from a CSV file.
//   - For each vertex in the humanoid OBJ, it interpolates the required insulation
//     thickness based on the X-coordinate (height).
//   - It then offsets each vertex along the +Y direction by the interpolated thickness.
//   - Finally, it writes a new OBJ file containing both the original and insulation
//     vertices.
//
// Assumptions:
//   - The first coordinate (X) of each OBJ vertex represents height along the robot.
//   - The thickness should be applied along the positive Y-axis.
//   - Linear interpolation is used between thickness points.
//
// Author: Parina Patel
// =============================================================================

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// Struct to represent a point from the thickness CSV.
struct ThicknessPoint {
  double height;     // X-position (height along the robot body).
  double thickness;  // Insulation thickness at this height.
};

// Struct to represent a vertex from the OBJ file.
struct Vertex {
  double x;
  double y;
  double z;
};

// Reads thickness data from a CSV file into a vector.
bool ReadCsv(const std::string& filename, std::vector<ThicknessPoint>* points) {
  std::ifstream infile(filename);
  if (!infile.is_open()) {
    std::cerr << "Failed to open CSV file: " << filename << std::endl;
    return false;
  }

  std::string line;
  // Skip the header line.
  std::getline(infile, line);

  while (std::getline(infile, line)) {
    std::istringstream iss(line);
    std::string token;
    ThicknessPoint point;

    if (!std::getline(iss, token, ',')) continue;
    point.height = std::stod(token);

    if (!std::getline(iss, token, ',')) continue;
    point.thickness = std::stod(token);

    points->push_back(point);
  }

  infile.close();
  return true;
}

// Reads vertices from an OBJ file, separating out non-vertex lines (e.g., faces).
bool ReadObj(const std::string& filename, std::vector<Vertex>* vertices,
             std::vector<std::string>* other_lines) {
  std::ifstream infile(filename);
  if (!infile.is_open()) {
    std::cerr << "Failed to open OBJ file: " << filename << std::endl;
    return false;
  }

  std::string line;
  while (std::getline(infile, line)) {
    if (line.empty()) continue;

    if (line.substr(0, 2) == "v ") {
      // Line describes a vertex.
      std::istringstream iss(line.substr(2));
      Vertex vertex;
      iss >> vertex.x >> vertex.y >> vertex.z;
      vertices->push_back(vertex);
    } else {
      // Store all other lines (faces, comments, etc.).
      other_lines->push_back(line);
    }
  }

  infile.close();
  return true;
}

// Interpolates insulation thickness at a given height (X-coordinate).
double InterpolateThickness(double height, const std::vector<ThicknessPoint>& points) {
  if (points.empty()) return 0.0;

  // Handle heights outside the given range.
  if (height <= points.front().height) {
    return points.front().thickness;
  }
  if (height >= points.back().height) {
    return points.back().thickness;
  }

  // Find two points between which 'height' lies and perform linear interpolation.
  for (size_t i = 0; i + 1 < points.size(); ++i) {
    if (height >= points[i].height && height <= points[i + 1].height) {
      double t = (height - points[i].height) /
                 (points[i + 1].height - points[i].height);
      return points[i].thickness +
             t * (points[i + 1].thickness - points[i].thickness);
    }
  }

  return 0.0;  // Default return; should not reach here if input is valid.
}

// Writes the combined original and offset vertices into a new OBJ file.
bool WriteCombinedObj(const std::string& filename,
                      const std::vector<Vertex>& original_vertices,
                      const std::vector<Vertex>& offset_vertices,
                      const std::vector<std::string>& other_lines) {
  std::ofstream outfile(filename);
  if (!outfile.is_open()) {
    std::cerr << "Failed to open output OBJ file: " << filename << std::endl;
    return false;
  }

  // Write original vertices.
  for (const auto& v : original_vertices) {
    outfile << "v " << v.x << " " << v.y << " " << v.z << "\n";
  }

  // Write offset (insulation) vertices.
  for (const auto& v : offset_vertices) {
    outfile << "v " << v.x << " " << v.y << " " << v.z << "\n";
  }

  // Write other OBJ elements (faces, comments, etc.).
  for (const auto& line : other_lines) {
    outfile << line << "\n";
  }

  outfile.close();
  return true;
}

// Main execution function.
int main() {
  // File paths.
  const std::string csv_filename = "Thickness_1_hr.csv";
  const std::string obj_filename = "humanoid_robot_2d.obj";
  const std::string output_filename = "humanoid_with_insulation.obj";

  // Read thickness profile from CSV.
  std::vector<ThicknessPoint> thickness_points;
  if (!ReadCsv(csv_filename, &thickness_points)) {
    return 1;
  }

  // Read humanoid vertices from OBJ.
  std::vector<Vertex> original_vertices;
  std::vector<std::string> other_lines;
  if (!ReadObj(obj_filename, &original_vertices, &other_lines)) {
    return 1;
  }

  // Generate offset vertices for insulation.
  std::vector<Vertex> offset_vertices;
  for (const auto& v : original_vertices) {
    double thickness = InterpolateThickness(v.x, thickness_points);
    Vertex v_offset = v;
    v_offset.y += thickness;  // Shift along +Y by the insulation thickness.
    offset_vertices.push_back(v_offset);
  }

  // Write new combined OBJ file.
  if (!WriteCombinedObj(output_filename, original_vertices, offset_vertices,
                        other_lines)) {
    return 1;
  }

  std::cout << "Successfully generated combined OBJ: " << output_filename
            << std::endl;
  return 0;
}
