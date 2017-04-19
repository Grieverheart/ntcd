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
#if 1
#define NTCD_IMPLEMENTATION
#include "ntcd.h"
#undef NTCD_IMPLEMENTATION

int main(int argc, char* argv[]){
    return 0;
}
#endif

#endif //__NTCD_H__

//Implementation
#ifdef NTCD_IMPLEMENTATION

#define BARY_GEPP

static inline double ntcd__vec3_dot(const double* a, const double* b){
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static inline void ntcd__vec3_cross(const double* a, const double* b, double* c){
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

static inline void ntcd__vec3_triple_product(const double* a, const double* b, const double* c, double* d){
    d[0] = b[0] * ntcd__vec3_dot(a, c) - a[0] * ntcd__vec3_dot(b, c);
    d[1] = b[1] * ntcd__vec3_dot(a, c) - a[1] * ntcd__vec3_dot(b, c);
    d[2] = b[2] * ntcd__vec3_dot(a, c) - a[2] * ntcd__vec3_dot(b, c);
}

static inline void ntcd__vec3_add(const double* a, const double* b, double* c){
    c[0] = a[0] + b[0];
    c[1] = a[1] + b[1];
    c[2] = a[2] + b[2];
}

static inline void ntcd__vec3_sub(const double* a, const double* b, double* c){
    c[0] = a[0] - b[0];
    c[1] = a[1] - b[1];
    c[2] = a[2] - b[2];
}

static const unsigned char p_pos[16][3] =
{
    {0, 0, 0}, {0, 0, 0}, {1, 0, 0}, {0, 1, 0},
    {2, 0, 0}, {0, 2, 0}, {1, 2, 0}, {0, 1, 2},
    {3, 0, 0}, {0, 3, 0}, {1, 3, 0}, {0, 1, 3},
    {2, 3, 0}, {0, 2, 3}, {1, 2, 3}, {0, 0, 0}
};

static const unsigned char s_pos[] = {0, 0, 1, 0, 2, 0, 0, 0, 3}; //Lookup table for single enabled bit position
//________________________^__^_____^___________^

#if defined(BARY_ERICSON)
static inline double triangle_area_2D(double x1, double y1, double x2, double y2, double x3, double y3){
    return (x1 - x2) * (y2 - y3) - (x2 - x3) * (y1 - y2);
}

//Algorithm for calculating barycentric coordinates from
//'Real-time collision detection' by Christer Ericson.
static inline Vec3d barycentric_coordinates(const double* P, const double* A, const double* B, const double* C, double* R){
    double u, v, w;

    double m[3];
    double AB[3], AC[3];
    ntcd__vec3_sub(B, A, AB);
    ntcd__vec3_sub(C, A, AC);
    ntcd__vec3_cross(AB, AC, m);

    double nu, nv, ood;
    double x = fabs(m[0]), y = fabs(m[1]), z = fabs(m[2]);

    if(x >= y && x >= z){
        nu = triangle_area_2D(P[1], P[2], B[1], B[2], C[1], C[2]);
        nv = triangle_area_2D(P[1], P[2], C[1], C[2], A[1], A[2]);
        ood = 1.0 / m[0];
    }
    else if(y >= x && y >= z){
        nu = triangle_area_2D(P[0], P[2], B[0], B[2], C[0], C[2]);
        nv = triangle_area_2D(P[0], P[2], C[0], C[2], A[0], A[2]);
        ood = 1.0 / -m[1];
    }
    else{
        nu = triangle_area_2D(P[0], P[1], B[0], B[1], C[0], C[1]);
        nv = triangle_area_2D(P[0], P[1], C[0], C[1], A[0], A[1]);
        ood = 1.0 / m[2];
    }

    R[0] = nu * ood;
    R[1] = nv * ood;
    R[2] = 1.0 - R[0] - R[1];
}

#elif defined(BARY_CRAMER)
static inline void barycentric_coordinates(const double* P, const double* A, const double* B, const double* C, double* R){
    double v0[3], v1[3], v2[3];
    ntcd__vec3_sub(B, A, v0);
    ntcd__vec3_sub(C, A, v1);
    ntcd__vec3_sub(P, A, v2);

    double d00 = ntcd__vec3_dot(v0, v0);
    double d01 = ntcd__vec3_dot(v0, v1);
    double d02 = ntcd__vec3_dot(v0, v2);
    double d11 = ntcd__vec3_dot(v1, v1);
    double d12 = ntcd__vec3_dot(v1, v2);
    double denom = d00 * d11 - d01 * d01;

    R[0] = (d11 * d02 - d01 * d12) / denom;
    R[1] = (d00 * d12 - d01 * d02) / denom;
    R[2] = 1.0 - R[0] - R[1];
}
#elif defined(BARY_GEPP)
static inline void barycentric_coordinates(const double* P, const double* A, const double* B, const double* C, double* R){
    double v0[3], v1[3], v2[3];
    ntcd__vec3_sub(B, A, v0);
    ntcd__vec3_sub(C, A, v1);
    ntcd__vec3_sub(P, A, v2);

    double d00 = ntcd__vec3_dot(v0, v0);
    double d01 = ntcd__vec3_dot(v0, v1);
    double d02 = ntcd__vec3_dot(v0, v2);
    double d11 = ntcd__vec3_dot(v1, v1);
    double d12 = ntcd__vec3_dot(v1, v2);

    R[0] = (d00 * d12 - d01 * d02) / (d00 * d11 - d01 * d01);
    R[1] = (d00 >= d01)? (d02 - d01 * R[0]) / d00: (d12 - d11 * R[0]) / d01;
    //R[1] = (d00 >= d01)? d02 / d00 - (d01 / d00) * w: d12 / d01 - (d11 / d01) * w;
    R[2] = 1.0 - R[0] - R[1];
}
#endif

#endif //NTCD_IMPLEMENTATION
