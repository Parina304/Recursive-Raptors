#ifndef _MESH_STRUCT
#define _MESH_STRUCT

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <set>
#include <tuple>

struct Point {
    double x, y, z;
};

typedef std::tuple<int, int, int, int> Tetrahedron;
typedef std::vector<Point> PointVector;
typedef std::set<Tetrahedron> TetrahedralMesh ;

#endif