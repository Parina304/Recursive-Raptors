#ifndef _MESH_MERGE_H
#define _MESH_MERGE_H

#include "mesh_lib.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/transform.h>
#include <CGAL/Aff_transformation_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> CGALMesh;
typedef Kernel::Aff_transformation_3 Transformation;

CGALMesh convertToCGALMesh(const Mesh& customMesh);
Mesh convertToNormalMesh(const CGALMesh& cgalMesh);
Mesh meshUnion(const Mesh& mesh1, const Mesh& mesh2);
bool caseInsensitiveCompare(const std::string& a, const std::string& b);
std::string getInputWithDefault(const std::string& prompt, const std::string& defaultValue);
#endif