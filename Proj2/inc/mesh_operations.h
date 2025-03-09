#ifndef _MESH_OPERARIONS
#define _MESH_OPERATIONS

#include "mesh_struct.h"

using namespace std;

void translate_mesh(PointVector &points, double x, double y, double z);
void rotate_mesh(PointVector &points, double rad, int axis);
void scale_mesh(PointVector& points, double scale_fac, double cx, double cy, double cz);

#endif