#include "mesh_lib.h"

void Mesh::translate(float x, float y, float z){
    for(auto &v: vertices){
        v.x += x;
        v.y += y;
        v.z += z;
    }
    centroid.x += x;
    centroid.y += y;
    centroid.z += z;
}

void Mesh::rotate(float rad, int axis){
    float x, y, z, xx, yy, zz;
    for(auto &v: vertices){
        x = v.x - centroid.x;
        y = v.y - centroid.y;
        z = v.z - centroid.z;
        if(!axis){
            xx = cos(rad) * x + sin(rad) * y;
            yy = -sin(rad) * x + cos(rad) * y;
            v.x = xx + centroid.x;
            v.y = yy + centroid.y;
        }
        else if(axis == 1){
            yy = cos(rad) * y + sin(rad) * z;
            zz = -sin(rad) * y + cos(rad) * z;
            v.y = yy + centroid.y;
            v.z = zz + centroid.z;
        }
        else if(axis == 2){
            zz = -sin(rad) * z + cos(rad) * x;
            xx = cos(rad) * z + sin(rad) * x;
            v.z = zz + centroid.z;
            v.x = xx + centroid.x;
        }
    }
}

void Mesh::scale(float scale_fac){
    for(auto &v: vertices){
        v.x = scale_fac * (v.x - centroid.x) + centroid.x;
        v.y = scale_fac * (v.y - centroid.y) + centroid.y;
        v.z = scale_fac * (v.z - centroid.z) + centroid.z;
    }
}

void Mesh::to_origin(){
    Mesh::translate(-centroid.x, -centroid.y, -centroid.z);
}