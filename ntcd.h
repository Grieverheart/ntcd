//Interface
#ifndef __NTCD_H__
#define __NTCD_H__

//TODO: Make transform_t customizable
typedef struct{
    double pos[3];
    double rot[4];
    double size;
}transform_t;

int ntcd_gjk_boolean(const transform_t*, const void*, const transform_t*, const void*);
void ntcd_gjk_distance(const transform_t*, const void*, const transform_t*, const void*, double* dist_vec);
int ntcd_gjk_raycast(const transform_t*, const void*, const transform_t*, const void*, const double* ray_dir, double* distance, double* normal);

//Test
#if 0
#define NTCD_IMPLEMENTATION
#include "ntcd.h"

int main(int argc, char* argv[]){
    return 0;
}
#endif

#endif //__NTCD_H__

//Implementation
#ifdef NTCD_IMPLEMENTATION
#endif //NTCD_IMPLEMENTATION
