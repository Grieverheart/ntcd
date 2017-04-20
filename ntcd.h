//Interface
#ifndef __NTCD_H__
#define __NTCD_H__

//TODO: Make ntcd_transform_t customizable
typedef struct{
    double pos[3];
    double rot[4];
    double size;
}ntcd_transform_t;

int ntcd_gjk_boolean(const ntcd_transform_t*, const void*, const ntcd_transform_t*, const void*);
void ntcd_gjk_distance(const ntcd_transform_t*, const void*, const ntcd_transform_t*, const void*, double* dist_vec);
int ntcd_gjk_raycast(const ntcd_transform_t*, const void*, const ntcd_transform_t*, const void*, const double* ray_dir, double* distance, double* normal);

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

#include <string.h>

#define BARY_GEPP

static inline double ntcd__vec3_dot(const double* a, const double* b){
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static inline void ntcd__vec3_cross(double* c, const double* a, const double* b){
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

static inline void ntcd__vec3_triple_product(double* d, const double* a, const double* b, const double* c){
    d[0] = b[0] * ntcd__vec3_dot(a, c) - a[0] * ntcd__vec3_dot(b, c);
    d[1] = b[1] * ntcd__vec3_dot(a, c) - a[1] * ntcd__vec3_dot(b, c);
    d[2] = b[2] * ntcd__vec3_dot(a, c) - a[2] * ntcd__vec3_dot(b, c);
}

static inline void ntcd__vec3_add(double*c, const double* a, const double* b){
    c[0] = a[0] + b[0];
    c[1] = a[1] + b[1];
    c[2] = a[2] + b[2];
}

static inline void ntcd__vec3_sub(double* c, const double* a, const double* b){
    c[0] = a[0] - b[0];
    c[1] = a[1] - b[1];
    c[2] = a[2] - b[2];
}

static inline double ntcd__vec3_length2(const double* vec){
    return vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
}

//TODO: Try memcmp performance
static inline int ntcd__vec3_equal(const double* a, const double* b){
    return (a[0] == b[0]) && (a[1] == b[1]) && (a[2] == b[2]);
}

//Lookup table which tells us at which positions in the simplex array
//to get our points a, b, c, d. I.e. if our bits are 0111 -> 7 -> {0, 1, 2}
static const unsigned char p_pos[16][3] =
{
    {0, 0, 0}, {0, 0, 0}, {1, 0, 0}, {0, 1, 0},
    {2, 0, 0}, {0, 2, 0}, {1, 2, 0}, {0, 1, 2},
    {3, 0, 0}, {0, 3, 0}, {1, 3, 0}, {0, 1, 3},
    {2, 3, 0}, {0, 2, 3}, {1, 2, 3}, {0, 0, 0}
};

#if defined(BARY_ERICSON)
static inline double triangle_area_2D(double x1, double y1, double x2, double y2, double x3, double y3){
    return (x1 - x2) * (y2 - y3) - (x2 - x3) * (y1 - y2);
}

//Algorithm for calculating barycentric coordinates from
//'Real-time collision detection' by Christer Ericson.
static inline Vec3d barycentric_coordinates(double* R, const double* P, const double* A, const double* B, const double* C){
    double u, v, w;

    double m[3];
    {
        double AB[3], AC[3];
        ntcd__vec3_sub(AB, B, A);
        ntcd__vec3_sub(AC, C, A);
        ntcd__vec3_cross(m, AB, AC);
    }

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
static inline void barycentric_coordinates(double* R, const double* P, const double* A, const double* B, const double* C){
    double v0[3], v1[3], v2[3];
    ntcd__vec3_sub(v0, B, A);
    ntcd__vec3_sub(v1, C, A);
    ntcd__vec3_sub(v2, P, A);

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
static inline void barycentric_coordinates(double* R, const double* P, const double* A, const double* B, const double* C){
    double v0[3], v1[3], v2[3];
    ntcd__vec3_sub(v0, B, A);
    ntcd__vec3_sub(v1, C, A);
    ntcd__vec3_sub(v2, P, A);

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

typedef struct{
    unsigned char bits_;
    unsigned char last_sb_;
    unsigned char size_;
    double p_[3 * 4]; //up to 4 points / 3-Simplex
    double a_[3 * 4]; //up to 4 points / 3-Simplex
    double b_[3 * 4]; //up to 4 points / 3-Simplex
    double max_vert2_;
}ntcd__simplex_t;

//void ntcd__simplex_print(const ntcd__simplex_t* simplex){
//    unsigned char bits = simplex->bits_;
//    for(int i = 0; i < 4; ++i, bits >>= 1){
//        if(bits & 1) printf("%d: %f, %f, %f\n", i, simplex->p_[3 * i + 0], simplex->p_[3 * i + 1], simplex->p_[3 * i + 2]);
//    }
//}

static inline void ntcd__simplex_add_point(ntcd__simplex_t* simplex, const double* point){
    unsigned char b = ~simplex->bits_; //Flip bits
    b &= -b; //Last set (available) bit
    unsigned char pos = 1 << b; //Get the bit position from the lookup table
    simplex->last_sb_ = pos;
    simplex->bits_ |= b; //Insert the new bit
    ++simplex->size_;
    memcpy(simplex->p_ + 3 * pos, point, 3 * sizeof(*point));
    double l2 = ntcd__vec3_length2(point);
    if(l2 > simplex->max_vert2_) simplex->max_vert2_ = l2;
}

static inline void ntcd__simplex_add_point_with_info(ntcd__simplex_t* simplex, const double* point, const double* pa, const double* pb){
    ntcd__simplex_add_point(simplex, point);
    memcpy(simplex->a_ + 3 * simplex->last_sb_, pa, 3 * sizeof(*pa));
    memcpy(simplex->b_ + 3 * simplex->last_sb_, pb, 3 * sizeof(*pb));
}

static inline void ntcd__simplex_remove_point(ntcd__simplex_t* simplex, int p){
    simplex->bits_ ^= (1 << p); //Erase the bit at position p
    --simplex->size_;
}

int ntcd__simplex_contains(const ntcd__simplex_t* simplex, const double* point){
    unsigned char bits = simplex->bits_;
    for(int i = 0; i < 4; ++i, bits >>= 1){
        if((bits & 1) && ntcd__vec3_equal(simplex->p_ + 3 * i, point)) return 1;
    }
    return 0;
}

void ntcd__simplex_translate(ntcd__simplex_t* simplex, const double* dr){
    //for(int k = 0; k < 4; ++k) p_[k] += dr;
    simplex->max_vert2_ = 0.0;
    unsigned char bits = simplex->bits_;
    for(int i = 0; i < 4; ++i, bits >>= 1){
        if(bits & 1){
            ntcd__vec3_add(simplex->p_ + 3 * i, dr, simplex->p_ + 3 * i);
            if(ntcd__vec3_length2(simplex->p_ + 3 * i) > simplex->max_vert2_) simplex->max_vert2_ = ntcd__vec3_length2(simplex->p_ + 3 * i);
        }
    }
}

void ntcd__simplex_compute_closest_points(ntcd__simplex_t* simplex, double* pa, double* pb, const double* P){
    switch(simplex->size_){
    //IMPORTANT: We are having accuracy problems with this projection.
    case 3:{
        const unsigned char* pos = p_pos[(simplex->bits_ ^ (1 << simplex->last_sb_))];
        const double* aA = simplex->a_ + 3 * simplex->last_sb_;
        const double* aB = simplex->a_ + 3 * pos[0];
        const double* aC = simplex->a_ + 3 * pos[1];
        const double* bA = simplex->b_ + 3 * simplex->last_sb_;
        const double* bB = simplex->b_ + 3 * pos[0];
        const double* bC = simplex->b_ + 3 * pos[1];

        const double* A = simplex->p_ + 3 * simplex->last_sb_;
        const double* B = simplex->p_ + 3 * pos[0];
        const double* C = simplex->p_ + 3 * pos[1];

        double bary[3];
        barycentric_coordinates(bary, P, A, B, C);

        for(int i = 0; i < 3; ++i){
            pa[i] = aA[i] * bary[0] + aB[i] * bary[1] + aC[i] * bary[2];
            pb[i] = bA[i] * bary[0] + bB[i] * bary[1] + bC[i] * bary[2];
        }

        break;
    }
    case 2:{
        const unsigned char* pos = p_pos[(simplex->bits_ ^ (1 << simplex->last_sb_))];
        const double* aA = simplex->a_ + 3 * simplex->last_sb_;
        const double* aB = simplex->a_ + 3 * pos[0];
        const double* bA = simplex->b_ + 3 * simplex->last_sb_;
        const double* bB = simplex->b_ + 3 * pos[0];
        double u, v;
        {
            const double* A = simplex->p_ + 3 * simplex->last_sb_;
            const double* B = simplex->p_ + 3 * pos[0];

            {
                double AB[3], AP[3];
                ntcd__vec3_sub(AB, B, A);
                ntcd__vec3_sub(AP, P, A);
                v = ntcd__vec3_dot(AB, AP) / ntcd__vec3_length2(AB);
                u = 1.0 - v;
            }
        }

        for(int i = 0; i < 3; ++i){
            pa[i] = aA[i] * u + aB[i] * v;
            pb[i] = bA[i] * u + bB[i] * v;
        }

        break;
    }
    case 1:{
        memcpy(pa, simplex->a_ + 3 * simplex->last_sb_, 3 * sizeof(*pa));
        memcpy(pb, simplex->b_ + 3 * simplex->last_sb_, 3 * sizeof(*pb));
        break;
    }
    default:
        break;
    }
}

#endif //NTCD_IMPLEMENTATION
