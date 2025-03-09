#ifndef _MESH_IO
#define _MESH_IO

#include "mesh_struct.h"

bool read_obj(const std::string& filename, PointVector& points, TetrahedralMesh& mesh);
void write_mesh_as_obj(const PointVector& points, const TetrahedralMesh& mesh, const std::string& output_file);

#endif