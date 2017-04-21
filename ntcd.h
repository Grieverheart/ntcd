//Interface
#ifndef __NTCD_H__
#define __NTCD_H__

//NOTE: best way to name structs?
//NOTE: Can I define ntcd_transform in the implementation and still be visible?

//TODO: Make ntcd_transform customizable
typedef struct{
    double pos[3];
    double rot[4];
    double size;
}ntcd_transform;

typedef void (*ntcd_support)(double*, const void*, const double*);

int ntcd_gjk_boolean(const ntcd_transform*, const void*, const ntcd_transform*, const void*);
void ntcd_gjk_distance(const ntcd_transform*, const void*, const ntcd_transform*, const void*, double* dist_vec);
int ntcd_gjk_raycast(const ntcd_transform*, const void*, const ntcd_transform*, const void*, const double* ray_dir, double* distance, double* normal);

typedef struct{
    ntcd_support support;
    double base_radius_, half_height_;
}ntcd_cylinder;
void ntcd_init_cylinder(ntcd_cylinder* cyl, double base_radius, double height);

//Test
#if 1
#define NTCD_IMPLEMENTATION
#include "ntcd.h"
#undef NTCD_IMPLEMENTATION
#include <stdio.h>

int main(int argc, char* argv[]){
    ntcd_cylinder cyl;
    ntcd_init_cylinder(&cyl, 1.0, 1.0);

    ntcd_transform ta = {
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0, 1.0},
        1.0
    };

    ntcd_transform tb = {
        {2.01, 0.0, 0.0},
        {1.0, 0.0, 0.0, 1.0},
        1.0
    };

    int overlap = ntcd_gjk_boolean(&ta, &cyl, &tb, &cyl);
    printf("%d\n", overlap);

    return 0;
}
#endif

#endif //__NTCD_H__

//Implementation
#ifdef NTCD_IMPLEMENTATION

#include <string.h>
#include <math.h>

#define BARY_GEPP

//Vectors
static inline double ntcd__vec3_dot(const double* a, const double* b){
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static inline void ntcd__vec3_cross(double* c, const double* a, const double* b){
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

//Calculates d = (a x b) x c
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

static inline void ntcd__vec3_smul(double* c, double f, const double* a){
    c[0] = f * a[0];
    c[1] = f * a[1];
    c[2] = f * a[2];
}

static inline void ntcd__vec3_fmadd(double* c, double a, const double* x, const double* y){
    c[0] = a * x[0] + y[0];
    c[1] = a * x[1] + y[1];
    c[2] = a * x[2] + y[2];
}

static inline double ntcd__vec3_length2(const double* vec){
    return vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
}

//TODO: Try memcmp performance
static inline int ntcd__vec3_equal(const double* a, const double* b){
    return (a[0] == b[0]) && (a[1] == b[1]) && (a[2] == b[2]);
}

// Quaternions
//TODO: Do we really need to calculate the length? Aren't all quaternions assumed of unit length?
static inline void ntcd__quat_inverse(double* r, const double* q){
    double ilength2 = 1.0  / (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
    for(int i = 0; i < 3; ++i) r[i] = -q[i] * ilength2;
    r[3] = q[3] * ilength2;
}

static inline void ntcd__quat_vec3_rotate(double* r, const double* q, const double* v){
    double u[3];
    {
        double a[3], b[3];
        ntcd__vec3_cross(a, q, v);
        ntcd__vec3_fmadd(b, q[3], v, a);
        ntcd__vec3_cross(u, q, b);
    }

    ntcd__vec3_fmadd(r, 2.0, u, v);

    //TODO: Check if this is faster
    //r[0] = v[0] + 2.0 * (q[1] * (q[0] * v[1] - q[1] * v[0] + q[3] * v[2]) - q[2] * (q[2] * v[0] - q[0] * v[2] + q[3] * v[1]));
    //r[1] = v[1] + 2.0 * (q[2] * (q[1] * v[2] - q[2] * v[1] + q[3] * v[0]) - q[0] * (q[0] * v[1] - q[1] * v[0] + q[3] * v[2]));
    //r[2] = v[2] + 2.0 * (q[0] * (q[2] * v[0] - q[0] * v[2] + q[3] * v[1]) - q[1] * (q[1] * v[2] - q[2] * v[1] + q[3] * v[0]));
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
}ntcd__simplex;

static inline void ntcd__simplex_init(ntcd__simplex* simplex){
    simplex->bits_      = 0;
    simplex->last_sb_   = 0;
    simplex->size_      = 0;
    simplex->max_vert2_ = 0.0;
}

//void ntcd__simplex_print(const ntcd__simplex* simplex){
//    unsigned char bits = simplex->bits_;
//    for(int i = 0; i < 4; ++i, bits >>= 1){
//        if(bits & 1) printf("%d: %f, %f, %f\n", i, simplex->p_[3 * i + 0], simplex->p_[3 * i + 1], simplex->p_[3 * i + 2]);
//    }
//}

static inline void ntcd__simplex_add_point(ntcd__simplex* simplex, const double* point){
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

static inline void ntcd__simplex_add_point_with_info(ntcd__simplex* simplex, const double* point, const double* pa, const double* pb){
    ntcd__simplex_add_point(simplex, point);
    memcpy(simplex->a_ + 3 * simplex->last_sb_, pa, 3 * sizeof(*pa));
    memcpy(simplex->b_ + 3 * simplex->last_sb_, pb, 3 * sizeof(*pb));
}

static inline void ntcd__simplex_remove_point(ntcd__simplex* simplex, int p){
    simplex->bits_ ^= (1 << p); //Erase the bit at position p
    --simplex->size_;
}

static inline int ntcd__simplex_contains(const ntcd__simplex* simplex, const double* point){
    unsigned char bits = simplex->bits_;
    for(int i = 0; i < 4; ++i, bits >>= 1){
        if((bits & 1) && ntcd__vec3_equal(simplex->p_ + 3 * i, point)) return 1;
    }
    return 0;
}

static inline void ntcd__simplex_translate(ntcd__simplex* simplex, const double* dr){
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

static void ntcd__simplex_compute_closest_points(ntcd__simplex* simplex, double* pa, double* pb, const double* P){
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

static void ntcd__simplex_closest(ntcd__simplex* simplex, double* dir){
    ///////////////////////////////////////////////
    //  Check if the origin is contained in the  //
    //  Minkowski sum.                           //
    ///////////////////////////////////////////////
    switch(simplex->size_){
    case 4:
    {
        const unsigned char* pos = p_pos[(simplex->bits_ ^ (1 << simplex->last_sb_))];

        const double* a = simplex->p_ + 3 * simplex->last_sb_;
        const double* b = simplex->p_ + 3 * pos[0];
        const double* c = simplex->p_ + 3 * pos[1];
        const double* d = simplex->p_ + 3 * pos[2];

        double ab[3], ac[3], ad[3];
        ntcd__vec3_sub(ab, b, a);
        ntcd__vec3_sub(ac, c, a);
        ntcd__vec3_sub(ad, d, a);

        ////////////////////* Vertex Case *///////////////////

        double dot_aba = ntcd__vec3_dot(ab, a);
        double dot_aca = ntcd__vec3_dot(ac, a);
        double dot_ada = ntcd__vec3_dot(ad, a);

        if(dot_aba >= 0.0 && dot_aca >= 0.0 && dot_ada >= 0.0){
            memcpy(dir, a, 3 * sizeof(*dir)); //Take direction passing through origin
            ntcd__simplex_remove_point(simplex, pos[0]);
            ntcd__simplex_remove_point(simplex, pos[1]);
            ntcd__simplex_remove_point(simplex, pos[2]);
            break;
        }

        ////////////////////* Edge Cases *///////////////////

        /* ab Edge case */
        double dot_abb = ntcd__vec3_dot(ab, b);
        double dot_abPerp1 = dot_aba * ntcd__vec3_dot(ac, b) - dot_abb * dot_aca;
        double dot_abPerp2 = dot_aba * ntcd__vec3_dot(ad, b) - dot_abb * dot_ada;
        // The origin must be inside the space defined by the intersection
        // of two half-space normal to the adjacent faces abc, abd
        if(dot_abPerp1 <= 0.0 && dot_abPerp2 <= 0.0 && dot_aba <= 0.0){
            double f = dot_aba / (dot_aba - dot_abb);
            ntcd__vec3_fmadd(dir, f, ab, a);
            ntcd__simplex_remove_point(simplex, pos[1]);
            ntcd__simplex_remove_point(simplex, pos[2]);
            break;
        }

        /* ac Edge case */
        double dot_acc = ntcd__vec3_dot(ac, c);
        double dot_acPerp1 = dot_aca * ntcd__vec3_dot(ad, c) - dot_acc * dot_ada;
        double dot_acPerp2 = dot_aca * ntcd__vec3_dot(ab, c) - dot_acc * dot_aba;
        // The origin must be inside the space defined by the intersection
        // of two half-space normal to the adjacent faces abc, acd
        if(dot_acPerp1 <= 0.0 && dot_acPerp2 <= 0.0 && dot_aca <= 0.0){
            double f = dot_aca / (dot_aca - dot_acc);
            ntcd__vec3_fmadd(dir, f, ac, a);
            ntcd__simplex_remove_point(simplex, pos[0]);
            ntcd__simplex_remove_point(simplex, pos[2]);
            break;
        }

        /* ad Edge case */
        double dot_add = ntcd__vec3_dot(ad, d);
        double dot_adPerp1 = dot_ada * ntcd__vec3_dot(ab, d) - dot_add * dot_aba;
        double dot_adPerp2 = dot_ada * ntcd__vec3_dot(ac, d) - dot_add * dot_aca;
        // The origin must be inside the space defined by the intersection
        // of two half-space normal to the adjacent faces acd, abd
        if(dot_adPerp1 <= 0.0 && dot_adPerp2 <= 0.0 && dot_ada <= 0.0){
            double f = dot_ada / (dot_ada - dot_add);
            ntcd__vec3_fmadd(dir, f, ad, a);
            ntcd__simplex_remove_point(simplex, pos[0]);
            ntcd__simplex_remove_point(simplex, pos[1]);
            break;
        }

        ////////////////////* Face Cases *///////////////////

        /* On abc side */
        // The origin should be on abc's side and between the half-spaces defined by ac and ab (normal to abc)
        {
            double abxac[3];
            ntcd__vec3_cross(abxac, ab, ac);
            if(ntcd__vec3_dot(ad, abxac) * ntcd__vec3_dot(a, abxac) > 0.0 && dot_abPerp1 >= 0.0 && dot_acPerp2 >= 0.0){
                /* Remove point d */
                ntcd__simplex_remove_point(simplex, pos[2]);
                double f = ntcd__vec3_dot(dir, a) / ntcd__vec3_length2(dir);
                if(ntcd__vec3_dot(ad, abxac) > 0.0) ntcd__vec3_smul(dir, -f, abxac);
                else ntcd__vec3_smul(dir, f, abxac);
                break;
            }
        }

        /* On abd side */
        // The origin should be on abd's side and between the half-spaces defined by ab and ad (normal to abd)
        {
            double abxad[3];
            ntcd__vec3_cross(abxad, ab, ad);
            if(ntcd__vec3_dot(ac, abxad) * ntcd__vec3_dot(a, abxad) > 0.0 && dot_abPerp2 >= 0.0 && dot_adPerp1 >= 0.0){
                /* Remove point c */
                ntcd__simplex_remove_point(simplex, pos[1]);
                double f = ntcd__vec3_dot(dir, a) / ntcd__vec3_length2(dir);
                if(ntcd__vec3_dot(ac, abxad) > 0.0) ntcd__vec3_smul(dir, -f, abxad);
                else ntcd__vec3_smul(dir, f, abxad);
                break;
            }
        }

        ///* On acd side */
        // The origin should be on acd's side and between the half-spaces defined by ac and ad (normal to acd)
        {
            double acxad[3];
            ntcd__vec3_cross(acxad, ac, ad);
            if(ntcd__vec3_dot(ab, acxad) * ntcd__vec3_dot(a, acxad) > 0.0 && dot_acPerp1 >= 0.0 && dot_adPerp2 >= 0.0){
                /* Remove point b */
                ntcd__simplex_remove_point(simplex, pos[0]);
                double f = ntcd__vec3_dot(dir, a) / ntcd__vec3_length2(dir);
                if(ntcd__vec3_dot(ab, acxad) > 0.0) ntcd__vec3_smul(dir, -f, acxad);
                else ntcd__vec3_smul(dir, f, acxad);
                break;
            }
        }

        ///* On bcd side */
        //// The origin should be on bcd's side
        {
            double bc[3], bd[3], bcxbd[3];
            ntcd__vec3_sub(bc, c, b);
            ntcd__vec3_sub(bd, d, b);
            ntcd__vec3_cross(bcxbd, bc, bd);
            if(ntcd__vec3_dot(bcxbd, ab) * ntcd__vec3_dot(bcxbd, b) < 0.0){
                /* Remove point a */
                ntcd__simplex_remove_point(simplex, simplex->last_sb_);
                simplex->last_sb_ = pos[0];
                double f = ntcd__vec3_dot(dir, b) / ntcd__vec3_length2(dir);
                if(ntcd__vec3_dot(ab, bcxbd) < 0.0) ntcd__vec3_smul(dir, -f, bcxbd);
                else ntcd__vec3_smul(dir, f, bcxbd);
                break;
            }
        }

        /* 'else' should only be when the origin is inside the tetrahedron */
        break;
    }
    case 3:
    {
        const unsigned char* pos = p_pos[(simplex->bits_ ^ (1 << simplex->last_sb_))];

        const double* a = simplex->p_ + 3 * simplex->last_sb_;
        const double* b = simplex->p_ + 3 * pos[0];
        const double* c = simplex->p_ + 3 * pos[1];

        double ab[3], ac[3];
        ntcd__vec3_sub(ab, b, a);
        ntcd__vec3_sub(ac, c, a);

        // Check if O in vertex region A
        double dot_aba = -ntcd__vec3_dot(ab, a);
        double dot_aca = -ntcd__vec3_dot(ac, a);
        if(dot_aba <= 0.0 && dot_aca <= 0.0){
            memcpy(dir, a, 3 * sizeof(*dir)); //Take direction passing through origin
            ntcd__simplex_remove_point(simplex, pos[0]);
            ntcd__simplex_remove_point(simplex, pos[1]);
            break;
        }

        // Check if O in edge region AB
        double dot_abb = -ntcd__vec3_dot(ab, b);
        double dot_acb = -ntcd__vec3_dot(ac, b);
        double vc = dot_aba * dot_acb - dot_abb * dot_aca;
        if(vc <= 0.0 && dot_aba >= 0.0 && dot_abb <= 0.0){
            double f = dot_aba / (dot_aba - dot_abb);
            ntcd__vec3_fmadd(dir, f, ab, a);
            /* Remove Point c */
            ntcd__simplex_remove_point(simplex, pos[1]);
            break;
        }

        // Check if O in edge region AC
        double dot_abc = -ntcd__vec3_dot(ab, c);
        double dot_acc = -ntcd__vec3_dot(ac, c);
        double vb = dot_abc * dot_aca - dot_aba * dot_acc;
        if(vb <= 0.0 && dot_aca >= 0.0 && dot_acc <= 0.0){
            double f = dot_aca / (dot_aca - dot_acc);
            ntcd__vec3_fmadd(dir, f, ac, a);
            /* Remove Point b */
            ntcd__simplex_remove_point(simplex, pos[0]);
            break;
        }

        double va = dot_abb * dot_acc - dot_abc * dot_acb;
        double w = 1.0 / (va + vb + vc);
        for(int i = 0; i < 3; ++i) dir[i] = a[i] + (ab[i] * vb + ac[i] * vc) * w;
        break;
    }
    case 2:
    {
        const unsigned char* pos = p_pos[(simplex->bits_ ^ (1 << simplex->last_sb_))];

        const double* a = simplex->p_ + 3 * simplex->last_sb_;
        const double* b = simplex->p_ + 3 * pos[0];

        double  ab[3];
        ntcd__vec3_add(ab, b, a);

        double t = -ntcd__vec3_dot(ab, a);
        if(t <= 0.0){
            memcpy(dir, a, 3 * sizeof(*dir)); //Take direction passing through origin
            ntcd__simplex_remove_point(simplex, pos[0]);
            break;
        }

        double denom = ntcd__vec3_length2(ab);
        if(t >= denom){
            ntcd__simplex_remove_point(simplex, simplex->last_sb_);
            simplex->last_sb_ = pos[0];
            memcpy(dir, b, 3 * sizeof(*dir));
            break;
        }

        ntcd__vec3_fmadd(dir, t / denom, ab, a);
        break;
    }
    case 1:
    {
        memcpy(dir, simplex->p_ + 3 * simplex->last_sb_, 3 * sizeof(*dir));
        break;
    }
    default: break;
    }
}

//TODO: dir -> -dir
static int ntcd__simplex_contains_origin(ntcd__simplex* simplex, double* dir){
    ///////////////////////////////////////////////
    //  Check if the origin is contained in the  //
    //  Minkowski sum.                           //
    ///////////////////////////////////////////////
    switch(simplex->size_){
    case 4:
    {
        const unsigned char* pos = p_pos[(simplex->bits_ ^ (1 << simplex->last_sb_))];

        const double* a = simplex->p_ + 3 * simplex->last_sb_;
        const double* b = simplex->p_ + 3 * pos[0];
        const double* c = simplex->p_ + 3 * pos[1];
        const double* d = simplex->p_ + 3 * pos[2];

        double ab[3], ac[3], ad[3];
        ntcd__vec3_sub(ab, b, a);
        ntcd__vec3_sub(ac, c, a);
        ntcd__vec3_sub(ad, d, a);

        ////////////////////* Face Cases *///////////////////

        /* On abc side */
        double abxac[3];
        ntcd__vec3_cross(abxac, ab, ac);
        int abPerp1Pos, acPerp2Pos;
        {
            double cross_temp_a[3], cross_temp_b[3];
            ntcd__vec3_cross(cross_temp_a, abxac, ab);
            ntcd__vec3_cross(cross_temp_b, ac, abxac);
            abPerp1Pos = (ntcd__vec3_dot(cross_temp_a, a) > 0.0);
            acPerp2Pos = (ntcd__vec3_dot(cross_temp_b, a) > 0.0);
        }
        // The origin should be on abc's side and between the half-spaces defined by ac and ab (normal to abc)
        {
            double abcPerp[3];
            if(ntcd__vec3_dot(abxac, ad) > 0.0) ntcd__vec3_smul(abcPerp, -1.0, abxac);
            else memcpy(abcPerp, abxac, 3 * sizeof(*abcPerp));

            if((ntcd__vec3_dot(abcPerp, a) < 0.0) && !abPerp1Pos && !acPerp2Pos){
                /* Remove point d */
                ntcd__simplex_remove_point(simplex, pos[2]);
                memcpy(dir, abcPerp, 3 * sizeof(*dir));
                break;
            }
        }

        /* On abd side */
        double abxad[3];
        ntcd__vec3_cross(abxad, ab, ad);
        int abPerp2Pos, adPerp1Pos;
        {
            double cross_temp_a[3], cross_temp_b[3];
            ntcd__vec3_cross(cross_temp_a, abxad, ab);
            ntcd__vec3_cross(cross_temp_b, ad, abxad);
            abPerp2Pos = (ntcd__vec3_dot(cross_temp_a, a) > 0.0);
            adPerp1Pos = (ntcd__vec3_dot(cross_temp_b, a) > 0.0);
        }
        // The origin should be on abd's side and between the half-spaces defined by ab and ad (normal to abd)
        {
            double abdPerp[3];
            if(ntcd__vec3_dot(abxad, ac) > 0.0) ntcd__vec3_smul(abdPerp, -1.0, abxad);
            else memcpy(abdPerp, abxad, 3 * sizeof(*abdPerp));

            if((ntcd__vec3_dot(abdPerp, a) < 0.0) && !abPerp2Pos && !adPerp1Pos){
                /* Remove point c */
                ntcd__simplex_remove_point(simplex, pos[1]);
                memcpy(dir, abdPerp, 3 * sizeof(*dir));
                break;
            }
        }

        /* On acd side */
        double acxad[3];
        ntcd__vec3_cross(acxad, ac, ad);
        int acPerp1Pos, adPerp2Pos;
        {
            double cross_temp_a[3], cross_temp_b[3];
            ntcd__vec3_cross(cross_temp_a, acxad, ac);
            ntcd__vec3_cross(cross_temp_b, ad, acxad);
            acPerp1Pos = (ntcd__vec3_dot(cross_temp_a, a) > 0.0);
            adPerp2Pos = (ntcd__vec3_dot(cross_temp_b, a) > 0.0);
        }
        // The origin should be on acd's side and between the half-spaces defined by ac and ad (normal to acd)
        {
            double acdPerp[3];
            if(ntcd__vec3_dot(acxad, ab) > 0.0) ntcd__vec3_smul(acdPerp, -1.0, acxad);
            else memcpy(acdPerp, acxad, 3 * sizeof(*acdPerp));

            if((ntcd__vec3_dot(acdPerp, a) < 0.0) && !acPerp1Pos && !adPerp2Pos){
                /* Remove point b */
                ntcd__simplex_remove_point(simplex, pos[0]);
                memcpy(dir, acdPerp, 3 * sizeof(*dir));
                break;
            }
        }

        ////////////////////* Edge Cases *///////////////////

        /* ab Edge case */
        // The origin must be inside the space defined by the intersection
        // of two half-space normal to the adjacent faces abc, abd
        if(abPerp1Pos && abPerp2Pos){
            ntcd__vec3_triple_product(dir, a, ab, ab);
            ntcd__simplex_remove_point(simplex, pos[1]);
            ntcd__simplex_remove_point(simplex, pos[2]);
            break;
        }

        /* ac Edge case */
        // The origin must be inside the space defined by the intersection
        // of two half-space normal to the adjacent faces abc, acd
        if(acPerp1Pos && acPerp2Pos){
            ntcd__vec3_triple_product(dir, a, ac, ac);
            ntcd__simplex_remove_point(simplex, pos[0]);
            ntcd__simplex_remove_point(simplex, pos[2]);
            break;
        }

        /* ad Edge case */
        // The origin must be inside the space defined by the intersection
        // of two half-space normal to the adjacent faces acd, abd
        if(adPerp1Pos && adPerp2Pos){
            ntcd__vec3_triple_product(dir, a, ad, ad);
            ntcd__simplex_remove_point(simplex, pos[0]);
            ntcd__simplex_remove_point(simplex, pos[1]);
            break;
        }

        /* 'else' should only be when the origin is inside the tetrahedron */
        return 1;
    }
    case 3:
    {
        const unsigned char* pos = p_pos[(simplex->bits_ ^ (1 << simplex->last_sb_))];

        const double* a = simplex->p_ + 3 * simplex->last_sb_;
        const double* b = simplex->p_ + 3 * pos[0];
        const double* c = simplex->p_ + 3 * pos[1];

        double ab[3], ac[3], abxac[3];
        ntcd__vec3_sub(ab, b, a);
        ntcd__vec3_sub(ac, c, a);
        ntcd__vec3_cross(abxac, ab, ac);

        ////////////////////* Edge Cases *///////////////////

        /* Origin on the outside of triangle and close to ab */
        double cross_temp[3];
        ntcd__vec3_cross(cross_temp, ab, abxac);
        if(ntcd__vec3_dot(cross_temp, a) < 0.0){
            ntcd__vec3_triple_product(dir, a, ab, ab);
            /* Remove Point c */
            ntcd__simplex_remove_point(simplex, pos[1]);
            break;
        }

        /* Origin on the outside of triangle and close to ac */
        ntcd__vec3_cross(cross_temp, abxac, ac);
        if(ntcd__vec3_dot(cross_temp, a) < 0.0){
            ntcd__vec3_triple_product(dir, a, ac, ac);
            /* Remove Point b */
            ntcd__simplex_remove_point(simplex, pos[0]);
            break;
        }

        /////////////////////* Face Case *///////////////////
        if(ntcd__vec3_dot(abxac, a) > 0.0) ntcd__vec3_smul(dir, -1.0, abxac);
        else memcpy(dir, abxac, 3 * sizeof(*dir));
        break;
    }
    case 2:
    {
        const unsigned char* pos = p_pos[(simplex->bits_ ^ (1 << simplex->last_sb_))];

        const double* a = simplex->p_ + 3 * simplex->last_sb_;
        const double* b = simplex->p_ + 3 * pos[0];

        double ab[3];
        ntcd__vec3_sub(ab, b, a);

        ntcd__vec3_triple_product(dir, a, ab, ab);
        break;
    }
    case 1:
    {
        const double* a = simplex->p_ + 3 * simplex->last_sb_;
        ntcd__vec3_smul(dir, -1.0, a);
        break;
    }
    default: break;
    }

    return 0;
}

int ntcd_gjk_boolean(
    const ntcd_transform* pa, const void* ca,
    const ntcd_transform* pb, const void* cb
){
    double dir[3];
    ntcd__vec3_sub(dir, pb->pos, pa->pos);
    ntcd__simplex simplex;
    ntcd__simplex_init(&simplex);

    unsigned int fail_safe = 0;

    double inv_rot_a[4], inv_rot_b[4];
    ntcd__quat_inverse(inv_rot_a, pa->rot);
    ntcd__quat_inverse(inv_rot_b, pb->rot);

    ntcd_support sa = *(const ntcd_support*)ca;
    ntcd_support sb = *(const ntcd_support*)cb;

    do{
        double vertex_a[3];
        {
            double inv_dir[3], support_point[3];
            ntcd__quat_vec3_rotate(inv_dir, inv_rot_a, dir);
            sa(support_point, ca, inv_dir);
            ntcd__quat_vec3_rotate(vertex_a, pa->rot, support_point);
            ntcd__vec3_fmadd(vertex_a, pa->size, vertex_a, pa->pos);
        }

        double vertex_b[3];
        {
            double inv_dir[3], support_point[3];
            ntcd__quat_vec3_rotate(inv_dir, inv_rot_b, dir);
            ntcd__vec3_smul(inv_dir, -1.0, inv_dir);
            sb(support_point, cb, inv_dir);
            ntcd__quat_vec3_rotate(vertex_b, pb->rot, support_point);
            ntcd__vec3_fmadd(vertex_b, pb->size, vertex_b, pb->pos);
        }

        double new_point[3];
        ntcd__vec3_sub(new_point, vertex_a, vertex_b);
        double dn = ntcd__vec3_dot(dir, new_point);
        if(dn < 0.0 || ntcd__simplex_contains(&simplex, new_point)) return 0;
        ntcd__simplex_add_point(&simplex, new_point);
        if(ntcd__simplex_contains_origin(&simplex, dir) || ntcd__vec3_length2(dir) == 0.0) return 1;
    }while(fail_safe++ < 100);

    //printf("Encountered error in GJK boolean: Infinite Loop.\n Direction (%f, %f, %f)\n", dir[0], dir[1], dir[2]);

    return 1;
}

//Shape implementations
static void ntcd__support_cylinder(double* support_point, const void* shape, const double* dir){
    ntcd_cylinder cyl = *(const ntcd_cylinder*)shape;

    double length = sqrt(dir[0] * dir[0] + dir[2] * dir[2]);
    if(length != 0.0){
        double d = cyl.base_radius_ / length;
        support_point[0] = d * dir[0];
        support_point[1] = copysign(cyl.half_height_, dir[1]);
        support_point[2] = d * dir[2];
    }
    else{
        support_point[0] = cyl.base_radius_;
        support_point[1] = copysign(cyl.half_height_, dir[1]);
        support_point[2] = 0.0;
    }
}

void ntcd_init_cylinder(ntcd_cylinder* cyl, double base_radius, double height){
    cyl->support = ntcd__support_cylinder;
    cyl->base_radius_ = base_radius;
    cyl->half_height_ = 0.5 * height;
}

#endif //NTCD_IMPLEMENTATION
