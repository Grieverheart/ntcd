//Interface
#ifndef __NTCD_H__
#define __NTCD_H__

//TODO: Define tolerances as macros.
//TODO: Implement more shapes.
//TODO: Investigate raycast false positives on y-axis.

#ifdef __cplusplus
extern "C"{
#endif

//TODO: Make ntcd_transform customizable.
typedef struct{
    double pos[3];
    double rot[4];
    double size;
}ntcd_transform;

typedef void (*ntcd_support)(double*, const void*, const double*);

int ntcd_gjk_boolean(const ntcd_transform*, const void*, const ntcd_transform*, const void*);
void ntcd_gjk_distance(double* dist_vec, const ntcd_transform*, const void*, const ntcd_transform*, const void*);
void ntcd_gjk_closest_points(double* point_on_a, double* point_on_b, const ntcd_transform*, const void*, const ntcd_transform*, const void*);
int ntcd_gjk_raycast(double* distance, double* normal, const ntcd_transform*, const void*, const ntcd_transform*, const void*, const double* ray_dir);

//Shapes
//NOTE: The long axis of a shape is defined to be along the y axis.
//TODO: Test all shapes

//Current supported shapes:
//Sphere, Point, Mesh, Cylinder, Box,
//Cone, Bicone, Leaf Cylinder, Hull

//Sphere
typedef struct{
    ntcd_support support;
}ntcd_sphere;
void ntcd_sphere_initialize(ntcd_sphere* sph);

//Point
typedef struct{
    ntcd_support support;
}ntcd_point;
void ntcd_point_initialize(ntcd_point* point);

//Disk
typedef struct{
    ntcd_support support;
}ntcd_disk;
void ntcd_disk_initialize(ntcd_disk* disk);

//Mesh
typedef struct{
    ntcd_support support;
    unsigned int n_vertices_;
    double* vertices_;
    unsigned int* n_vert_neighbours_;
    unsigned int** vert_neighbours_;
}ntcd_mesh;
void ntcd_mesh_initialize(ntcd_mesh* mesh, const unsigned int n_vertices, const double* vertices, const unsigned int n_faces, const unsigned int* face_start, const unsigned int* faces);
void ntcd_mesh_terminate(ntcd_mesh* mesh);

//Cylinder
typedef struct{
    ntcd_support support;
    double base_radius_, half_height_;
}ntcd_cylinder;
void ntcd_cylinder_initialize(ntcd_cylinder* cyl, double base_radius, double height);

//Box
typedef struct{
    ntcd_support support;
    double extent_[3];
}ntcd_box;
void ntcd_box_initialize(ntcd_box* box, const double* extent);

//Cone
typedef struct{
    ntcd_support support;
    double base_radius_, half_height_;
    double sintheta_;
}ntcd_cone;
void ntcd_cone_initialize(ntcd_cone* cone, double base_radius, double height);

//Bicone
typedef struct{
    ntcd_support support;
    double base_radius_, half_height_;
    double sintheta_;
}ntcd_bicone;
void ntcd_bicone_initialize(ntcd_bicone* bicone, double base_radius, double height);

//Leaf Cylinder
typedef struct{
    ntcd_support support;
    double half_width_, half_length_, half_height_;
    double circle_radius_, circle_distance_;
}ntcd_leaf_cylinder;
void ntcd_leaf_cylinder_initialize(ntcd_leaf_cylinder* leaf, double width, double length, double height);

//Hull
typedef struct{
    ntcd_support support;
    const void* ca_;
    const void* cb_;
    ntcd_transform ta_;
    ntcd_transform tb_;
    //Cache inverse rotations.
    double inv_rot_a_[4];
    double inv_rot_b_[4];
}ntcd_hull;
void ntcd_hull_initialize(ntcd_hull* hull, const ntcd_transform* ta, const void* ca, const ntcd_transform* tb, const void* cb);

#ifdef __cplusplus
}
#endif

//Test
#if 0
#define NTCD_IMPLEMENTATION
#include "ntcd.h"
#undef NTCD_IMPLEMENTATION
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif
#include <assert.h>

int main(int argc, char* argv[]){
    ntcd_cylinder cyl;
    ntcd_cylinder_initialize(&cyl, 1.0, 1.0);

    //Two cylinders of height 1, rotate by 45 degrees around the z-axis,
    //should overlap at x = 1 / sin(45).

    double d = 1.0 / sin(M_PI / 4.0);

    ntcd_transform ta = {
        {0.0, 0.0, 0.0},
        {0.0, 0.0, sin(M_PI / 8.0), cos(M_PI / 8.0)},
        1.0
    };

    ntcd_transform tb = {
        {d - 0.01, 0.0, 0.0},
        {0.0, 0.0, sin(M_PI / 8.0), cos(M_PI / 8.0)},
        1.0
    };

    ntcd_transform tc = {
        {d + 0.01, 0.0, 0.0},
        {0.0, 0.0, sin(M_PI / 8.0), cos(M_PI / 8.0)},
        1.0
    };

    int overlap_a = ntcd_gjk_boolean(&ta, &cyl, &tb, &cyl);
    int overlap_b = ntcd_gjk_boolean(&ta, &cyl, &tc, &cyl);
    assert(overlap_a != overlap_b);
    double dist[3];
    ntcd_gjk_distance(dist, &ta, &cyl, &tc, &cyl);
    double length = sqrt(dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2]);
    assert(length < 0.01);
    assert(length > 0.0);

    return 0;
}
#endif

#endif //__NTCD_H__

//Implementation
#ifdef NTCD_IMPLEMENTATION

#include <string.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>

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

static inline double ntcd__vec3_length(const double* vec){
    return sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}

//TODO: Try memcmp performance
static inline int ntcd__vec3_equal(const double* a, const double* b){
    return (a[0] == b[0]) && (a[1] == b[1]) && (a[2] == b[2]);
}

// Quaternions
//TODO: Do we really need to calculate the length? Aren't all quaternions assumed of unit length?
static inline void ntcd__quat_inverse(double* r, const double* q){
    double ilength2 = 1.0  / (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
    r[0] = -q[0] * ilength2;
    r[1] = -q[1] * ilength2;
    r[2] = -q[2] * ilength2;
    r[3] =  q[3] * ilength2;
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

//Lookup table for single enabled bit position -- Calculates the log2 of a 4-bit integer.
static const unsigned char s_pos[] = {0, 0, 1, 0, 2, 0, 0, 0, 3};
//_______________________________________^__^_____^___________^

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
static inline void barycentric_coordinates(double* R, const double* P, const double* A, const double* B, const double* C){
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

static inline void ntcd__simplex_initialize(ntcd__simplex* simplex){
    simplex->bits_      = 0;
    simplex->last_sb_   = 0;
    simplex->size_      = 0;
    simplex->max_vert2_ = 0.0;
}

//void ntcd__simplex_print(const ntcd__simplex* simplex){
//    printf("yo\n");
//    unsigned char bits = simplex->bits_;
//    for(int i = 0; i < 4; ++i, bits >>= 1){
//        if(bits & 1) printf("%d: %f, %f, %f\n", i, simplex->p_[3 * i + 0], simplex->p_[3 * i + 1], simplex->p_[3 * i + 2]);
//    }
//}

static inline void ntcd__simplex_add_point(ntcd__simplex* simplex, const double* point){
    unsigned char b = ~simplex->bits_; //Flip bits
    b &= -b; //Last set (available) bit
    unsigned char pos = s_pos[b]; //Get the bit position from the lookup table
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
                ntcd__vec3_smul(dir, ntcd__vec3_dot(abxac, a) / ntcd__vec3_length2(abxac), abxac);
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
                ntcd__vec3_smul(dir, ntcd__vec3_dot(abxad, a) / ntcd__vec3_length2(abxad), abxad);
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
                ntcd__vec3_smul(dir, ntcd__vec3_dot(acxad, a) / ntcd__vec3_length2(acxad), acxad);
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
                ntcd__vec3_smul(dir, ntcd__vec3_dot(bcxbd, b) / ntcd__vec3_length2(bcxbd), bcxbd);
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
        ntcd__vec3_sub(ab, b, a);

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
    ntcd__simplex_initialize(&simplex);

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
            ntcd__vec3_smul(inv_dir, -1.0, dir);
            ntcd__quat_vec3_rotate(inv_dir, inv_rot_b, inv_dir);
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
    }while(fail_safe++ < 2000);

    //printf("Encountered error in GJK boolean: Infinite Loop.\n Direction (%f, %f, %f)\n", dir[0], dir[1], dir[2]);

    return 1;
}

void ntcd_gjk_distance(
    double* dist,
    const ntcd_transform* pa, const void* ca,
    const ntcd_transform* pb, const void* cb
){
    double inv_rot_a[4], inv_rot_b[4];
    ntcd__quat_inverse(inv_rot_a, pa->rot);
    ntcd__quat_inverse(inv_rot_b, pb->rot);

    double dir[3] = {0.0};
    ntcd__simplex simplex;
    ntcd__simplex_initialize(&simplex);

    unsigned int fail_safe = 0;

    double dist2 = DBL_MAX;

    ntcd_support sa = *(const ntcd_support*)ca;
    ntcd_support sb = *(const ntcd_support*)cb;

    do{
        double vertex_a[3];
        {
            double inv_dir[3], support_point[3];
            ntcd__vec3_smul(inv_dir, -1.0, dir);
            ntcd__quat_vec3_rotate(inv_dir, inv_rot_a, inv_dir);
            sa(support_point, ca, inv_dir);
            ntcd__quat_vec3_rotate(vertex_a, pa->rot, support_point);
            ntcd__vec3_fmadd(vertex_a, pa->size, vertex_a, pa->pos);
        }

        double vertex_b[3];
        {
            double inv_dir[3], support_point[3];
            ntcd__quat_vec3_rotate(inv_dir, inv_rot_b, dir);
            sb(support_point, cb, inv_dir);
            ntcd__quat_vec3_rotate(vertex_b, pb->rot, support_point);
            ntcd__vec3_fmadd(vertex_b, pb->size, vertex_b, pb->pos);
        }
        double new_point[3];
        ntcd__vec3_sub(new_point, vertex_a, vertex_b);

        if(ntcd__simplex_contains(&simplex, new_point) || dist2 - ntcd__vec3_dot(dir, new_point) <= dist2 * 1.0e-8){
            memcpy(dist, dir, 3 * sizeof(*dist));
            return;
        }
        ntcd__simplex_add_point(&simplex, new_point);

        ntcd__simplex_closest(&simplex, dir);

        dist2 = ntcd__vec3_length2(dir);

        if(simplex.size_ == 4 || dist2 < 1.0e-12){
            memset(dist, 0, 3 * sizeof(*dist));
            return;
        }

    }while(++fail_safe < 200);

    //if(fail_safe == 2000) printf("Encountered error in GJK distance: Infinite Loop.\n Direction (%f, %f, %f)\n", dir[0], dir[1], dir[2]);

    memcpy(dist, dir, 3 * sizeof(*dist));
}

void ntcd_gjk_closest_points(
    double* point_on_a, double* point_on_b,
    const ntcd_transform* pa, const void* ca,
    const ntcd_transform* pb, const void* cb
){
    double inv_rot_a[4], inv_rot_b[4];
    ntcd__quat_inverse(inv_rot_a, pa->rot);
    ntcd__quat_inverse(inv_rot_b, pb->rot);

    double dir[3] = {0.0};
    ntcd__simplex simplex;
    ntcd__simplex_initialize(&simplex);

    unsigned int fail_safe = 0;

    double dist2 = DBL_MAX;

    ntcd_support sa = *(const ntcd_support*)ca;
    ntcd_support sb = *(const ntcd_support*)cb;

    do{
        double vertex_a[3];
        {
            double inv_dir[3], support_point[3];
            ntcd__vec3_smul(inv_dir, -1.0, dir);
            ntcd__quat_vec3_rotate(inv_dir, inv_rot_a, inv_dir);
            sa(support_point, ca, inv_dir);
            ntcd__quat_vec3_rotate(vertex_a, pa->rot, support_point);
            ntcd__vec3_fmadd(vertex_a, pa->size, vertex_a, pa->pos);
        }

        double vertex_b[3];
        {
            double inv_dir[3], support_point[3];
            ntcd__quat_vec3_rotate(inv_dir, inv_rot_b, dir);
            sb(support_point, cb, inv_dir);
            ntcd__quat_vec3_rotate(vertex_b, pb->rot, support_point);
            ntcd__vec3_fmadd(vertex_b, pb->size, vertex_b, pb->pos);
        }
        double new_point[3];
        ntcd__vec3_sub(new_point, vertex_a, vertex_b);

        if(ntcd__simplex_contains(&simplex, new_point) || dist2 - ntcd__vec3_dot(dir, new_point) <= dist2 * 1.0e-8){
            ntcd__simplex_compute_closest_points(&simplex, point_on_a, point_on_b, dir);
            return;
        }

        ntcd__simplex_add_point_with_info(&simplex, new_point, vertex_a, vertex_b);

        ntcd__simplex_closest(&simplex, dir);

        dist2 = ntcd__vec3_length2(dir);

        if(simplex.size_ == 4 || dist2 < 1.0e-12){
            memset(point_on_a, 0, 3 * sizeof(*point_on_a));
            memset(point_on_b, 0, 3 * sizeof(*point_on_b));
            return;
        }

    }while(++fail_safe < 2000);

    //if(fail_safe == 2000) printf("Encountered error in GJK closest points: Infinite Loop.\n Direction (%f, %f, %f)\n", dir[0], dir[1], dir[2]);

    ntcd__simplex_compute_closest_points(&simplex, point_on_a, point_on_b, dir);
}

//TODO: Do we need some way to check if a specific point has been already added to the simplex.
int ntcd_gjk_raycast(
    double* distance, double* normal,
    const ntcd_transform* pa, const void* ca,
    const ntcd_transform* pb, const void* cb,
    const double* ray_dir
)
{
    static const double etol = 10.0 * DBL_EPSILON;

    double inv_rot_a[4], inv_rot_b[4];
    ntcd__quat_inverse(inv_rot_a, pa->rot);
    ntcd__quat_inverse(inv_rot_b, pb->rot);

    double dir[3];
    ntcd__vec3_sub(dir, pb->pos, pa->pos);
    ntcd__simplex simplex;
    ntcd__simplex_initialize(&simplex);

    unsigned int fail_safe = 0;

    double x[3] = {0.0};
    double lambda = 0.0;

    double dist2 = DBL_MAX;

    ntcd_support sa = *(const ntcd_support*)ca;
    ntcd_support sb = *(const ntcd_support*)cb;

    do{
        double vertex_a[3];
        {
            double inv_dir[3], support_point[3];
            ntcd__vec3_smul(inv_dir, -1.0, dir);
            ntcd__quat_vec3_rotate(inv_dir, inv_rot_a, inv_dir);
            sa(support_point, ca, inv_dir);
            ntcd__quat_vec3_rotate(vertex_a, pa->rot, support_point);
            ntcd__vec3_fmadd(vertex_a, pa->size, vertex_a, pa->pos);
        }

        double vertex_b[3];
        {
            double inv_dir[3], support_point[3];
            ntcd__quat_vec3_rotate(inv_dir, inv_rot_b, dir);
            sb(support_point, cb, inv_dir);
            ntcd__quat_vec3_rotate(vertex_b, pb->rot, support_point);
            ntcd__vec3_fmadd(vertex_b, pb->size, vertex_b, pb->pos);
        }
        double new_point[3];
        ntcd__vec3_sub(new_point, vertex_a, vertex_b);

        double new_point_trans[3];
        ntcd__vec3_sub(new_point_trans, new_point, x);

        if(ntcd__vec3_dot(dir, new_point_trans) > 0.0){
            if(ntcd__vec3_dot(dir, ray_dir) >= 0.0) return 0;

            double delta = ntcd__vec3_dot(dir, new_point_trans) / ntcd__vec3_dot(dir, ray_dir);
            lambda -= delta;
            if(lambda > *distance) return 0;
            ntcd__vec3_smul(x, -lambda, ray_dir);
            if(ntcd__vec3_length2(x) > 100.0) return 0;
            ntcd__vec3_smul(normal, -1.0 / ntcd__vec3_length(dir), dir);
            double dr[3];
            ntcd__vec3_smul(dr, -delta, ray_dir);
            ntcd__simplex_translate(&simplex, dr);
        }

        ntcd__simplex_add_point(&simplex, new_point_trans);
        ntcd__simplex_closest(&simplex, dir);

        dist2 = ntcd__vec3_length2(dir);

        if(simplex.size_ == 4 || dist2 < etol * simplex.max_vert2_){
            *distance = lambda;
            return 1;
        }
    }while(fail_safe++ < 1000);
    *distance = lambda;
    //printf("Encountered error in GJK raycast: Infinite Loop.\n Direction (%f, %f, %f)\n", dir[0], dir[1], dir[2]);
    return 0;
}

//Shape implementations

//Cylinder
static void ntcd__support_cylinder(double* support_point, const void* shape, const double* dir){
    ntcd_cylinder cyl = *(const ntcd_cylinder*)shape;

    double length = sqrt(dir[0] * dir[0] + dir[2] * dir[2]);
    //TODO: Make length > small_number?
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

void ntcd_cylinder_initialize(ntcd_cylinder* cyl, double base_radius, double height){
    cyl->support      = ntcd__support_cylinder;
    cyl->base_radius_ = base_radius;
    cyl->half_height_ = 0.5 * height;
}

//Box
static void ntcd__support_box(double* support_point, const void* shape, const double* dir){
    const double* extent = ((const ntcd_box*)shape)->extent_;
    support_point[0] = copysign(extent[0], dir[0]);
    support_point[1] = copysign(extent[1], dir[1]);
    support_point[2] = copysign(extent[2], dir[2]);
}

void ntcd_box_initialize(ntcd_box* box, const double* extent){
    box->support = ntcd__support_box;
    memcpy(box->extent_, extent, 3 * sizeof(*extent));
}

//Cone
static void ntcd__support_cone(double* support_point, const void* shape, const double* dir){
    ntcd_cone cone = *(const ntcd_cone*)shape;

    double test = dir[1] / ntcd__vec3_length(dir);
    if(test >= cone.sintheta_){
        support_point[0] = 0.0;
        support_point[1] = cone.half_height_;
        support_point[2] = 0.0;
    }
    else if(test < cone.sintheta_ && (dir[0] != 0.0 || dir[2] != 0.0)){
        double length = sqrt(dir[0] * dir[0] + dir[2] * dir[2]);
        support_point[0] = cone.base_radius_ * dir[0] / length;
        support_point[1] = -cone.half_height_;
        support_point[2] = cone.base_radius_ * dir[2] / length;
    }
    else{
        support_point[0] = 0.0;
        support_point[1] = -cone.half_height_;
        support_point[2] = 0.0;
    }
}

void ntcd_cone_initialize(ntcd_cone* cone, double base_radius, double height){
    cone->support      = ntcd__support_cone;
    cone->base_radius_ = base_radius;
    cone->half_height_ = 0.5 * height;
    cone->sintheta_    = base_radius / sqrt(base_radius * base_radius + height * height);
}

//Bicone
static void ntcd__support_bicone(double* support_point, const void* shape, const double* dir){
    ntcd_bicone bicone = *(const ntcd_bicone*)shape;

    double test = dir[1] / ntcd__vec3_length(dir);
    if(test >= bicone.sintheta_){
        support_point[0] = 0.0;
        support_point[1] = bicone.half_height_;
        support_point[2] = 0.0;
    }
    else if(test <= -bicone.sintheta_){
        support_point[0] = 0.0;
        support_point[1] = -bicone.half_height_;
        support_point[2] = 0.0;
    }
    else{
        double factor = bicone.base_radius_ / sqrt(dir[0] * dir[0] + dir[2] * dir[2]);
        support_point[0] = factor * dir[0];
        support_point[1] = 0.0;
        support_point[2] = factor * dir[2];
    }
}

void ntcd_bicone_initialize(ntcd_bicone* bicone, double base_radius, double height){
    bicone->support      = ntcd__support_bicone;
    bicone->base_radius_ = base_radius;
    bicone->half_height_ = 0.5 * height;
    bicone->sintheta_    = base_radius / sqrt(base_radius * base_radius + bicone->half_height_ * bicone->half_height_);
}

//Leaf Cylinder
static void ntcd__support_leaf_cylinder(double* support_point, const void* shape, const double* dir){
    ntcd_leaf_cylinder leaf = *(const ntcd_leaf_cylinder*)shape;

    double x = 0.0, z = 0.0;
    if(dir[0] != 0.0 || dir[2] != 0.0){
        double l = sqrt(dir[0] * dir[0] + dir[2] * dir[2]);
        double test = dir[2] / l;
        if(test >= leaf.half_length_ / leaf.circle_radius_) z = leaf.half_length_;
        else if(test <= -leaf.half_length_ / leaf.circle_radius_) z = -leaf.half_length_;
        else{
            x = leaf.circle_radius_ * dir[0] / l - copysign(leaf.circle_distance_, dir[0]);
            z = leaf.circle_radius_ * dir[2] / l;
        }
    }

    support_point[0] = x;
    support_point[1] = copysign(leaf.half_height_, dir[1]);
    support_point[2] = z;
}

void ntcd_leaf_cylinder_initialize(ntcd_leaf_cylinder* leaf, double width, double length, double height){
    leaf->support          = ntcd__support_leaf_cylinder;
    leaf->half_width_      = 0.5 * width;
    leaf->half_length_     = 0.5 * length;
    leaf->half_height_     = 0.5 * height;
    leaf->circle_radius_   = 0.25 * (length * length + width * width) / width;
    leaf->circle_distance_ = 0.25 * (length * length - width * width) / width;
}

//Sphere
static void ntcd__support_sphere(double* support_point, const void* shape, const double* dir){
    double norm = ntcd__vec3_length(dir);
    if(norm > 0.0){
        support_point[0] = dir[0] / norm;
        support_point[1] = dir[1] / norm;
        support_point[2] = dir[2] / norm;
    }
    else{
        support_point[0] = 1.0;
        support_point[1] = 0.0;
        support_point[2] = 0.0;
    }
}

void ntcd_sphere_initialize(ntcd_sphere* sph){
    sph->support = ntcd__support_sphere;
}

//Point
static void ntcd__support_point(double* support_point, const void* shape, const double* dir){
    memset(support_point, 0, 3 * sizeof(*dir));
}

void ntcd_point_initialize(ntcd_point* point){
    point->support = ntcd__support_point;
}

//Disk
static void ntcd__support_disk(double* support_point, const void* shape, const double* dir){
    double length2 = dir[0] * dir[0] + dir[2] * dir[2];
    if(length2 != 0.0){
        double length = 1.0 / sqrt(length2);
        support_point[0] = dir[0] * length;
        support_point[1] = 0.0;
        support_point[2] = dir[2] * length;
    }
    else{
        support_point[0] = 1.0;
        support_point[1] = 0.0;
        support_point[2] = 0.0;
    }
}

void ntcd_disk_initialize(ntcd_disk* disk){
    disk->support = ntcd__support_disk;
}

//Mesh

//TODO: Perhaps choose one of the two code paths depending on the number of vertices.
static void ntcd__support_mesh(double* support_point, const void* shape, const double* dir){
    ntcd_mesh mesh = *(const ntcd_mesh*)shape;
    //unsigned int curr = 0;
    //double p = 0.0;
    //double max = ntcd__vec3_dot(mesh.vertices_, dir);
    //for(unsigned int i = 1; i < mesh.n_vertices_; ++i){
    //    p = ntcd__vec3_dot(mesh.vertices_ + 3 * i, dir);
    //    if(p > max){
    //        curr = i;
    //        max = p;
    //    }
    //}
    //memcpy(support_point, mesh.vertices_ + 3 * curr, 3 * sizeof(support_point));
    unsigned int next = 0, last = 0, curr = 0;
    double p = 0.0;
    double max = ntcd__vec3_dot(mesh.vertices_, dir);
    for(;;){
        for(unsigned int vid = 0; vid < mesh.n_vert_neighbours_[curr]; ++vid){
            next = mesh.vert_neighbours_[curr][vid];
            if(next != last){
                p = ntcd__vec3_dot(mesh.vertices_ + 3 * next, dir);
                if(p > max){
                    max = p;
                    last = curr;
                    curr = next;
                    break;
                }
            }
            if(vid == mesh.n_vert_neighbours_[curr] - 1){
                memcpy(support_point, mesh.vertices_ + 3 * curr, 3 * sizeof(*support_point));
                return;
            }
        }
    }
}

void ntcd_mesh_initialize(ntcd_mesh* mesh, const unsigned int n_vertices, const double* vertices, const unsigned int n_faces, const unsigned int* face_start, const unsigned int* faces){
    mesh->support     = ntcd__support_mesh;
    mesh->n_vertices_ = n_vertices;
    mesh->vertices_   = (double*)malloc(3 * n_vertices * sizeof(*vertices));
    memcpy(mesh->vertices_, vertices, 3 * n_vertices * sizeof(*vertices));

    /* Find Edges */
    unsigned int n_edges = n_vertices + n_faces - 2; //Euler's formula.
    unsigned int* edges = (unsigned int*)malloc(n_edges * 2 * sizeof(*edges)); //Each edge is represented by its two vertex indices.
    //Iterate over each face and for each next face in the list, check if they
    //share two vertices, this defines an edge.
    for(unsigned int fi = 0, eid = 0; fi < n_faces; ++fi){
        const unsigned int* face = faces + face_start[fi];
        unsigned int face_size = ((fi < n_faces - 1)? face_start[fi + 1]: n_vertices) - face_start[fi];
        double normal[3];
        {
            double ab[3], ac[3];
            ntcd__vec3_sub(ab, vertices + 3 * face[1], vertices + 3 * face[0]);
            ntcd__vec3_sub(ac, vertices + 3 * face[face_size - 1], vertices + 3 * face[0]);
            ntcd__vec3_cross(normal, ab, ac);
        }
        ntcd__vec3_smul(normal, ntcd__vec3_length(normal), normal);

        for(unsigned int fj = fi + 1; fj < n_faces; ++fj){
            unsigned int fcount = 0;
            unsigned int edge[3];
            for(unsigned int i = 0; i < face_size; ++i){
                unsigned int vid_fi = face[i];
                for(unsigned int j = 0; j < face_size; ++j){
                    if(vid_fi == faces[face_start[fj] + j]){
                        edge[fcount++] = vid_fi;
                    }
                }
                if(fcount == 2){
                    memcpy(edges + 2 * eid, edge, 2 * sizeof(*edge));
                    ++eid;
                    fcount = 0;
                }
            }
        }
    }

    /* Find Vertex Neighbours */
    //For all vertices, check if two edges share this vertex. If they do and it
    //isn't vertex 0, append the other vertices of these edge to the neighbor list
    mesh->vert_neighbours_   = (unsigned int**)malloc(n_vertices * sizeof(*mesh->vert_neighbours_));
    mesh->n_vert_neighbours_ = (unsigned int*)calloc(n_vertices, sizeof(*mesh->n_vert_neighbours_));
    for(unsigned int vid = 0; vid < n_vertices; ++vid){
        unsigned int capacity = 5;
        unsigned int n_neighbours = 0;
        mesh->vert_neighbours_[vid] = (unsigned int *)malloc(capacity * sizeof(**mesh->vert_neighbours_));
        unsigned int* neighbours = mesh->vert_neighbours_[vid];
        for(unsigned int ei = 0; ei < n_edges; ++ei){
            for(unsigned int i = 0; i < 2; ++i){
                if(edges[2 * ei + i] == vid && edges[2 * ei + (i + 1) % 2] != 0){
                    neighbours[n_neighbours++] = edges[2 * ei + (i + 1) % 2];
                    if(n_neighbours == capacity){
                        capacity *= 2;
                        unsigned int* temp = (unsigned int*)realloc(mesh->vert_neighbours_[vid], capacity * sizeof(**mesh->vert_neighbours_));
                        //assert(temp != NULL);
                        mesh->vert_neighbours_[vid] = temp;
                        neighbours = temp;
                    }
                }
            }
        }
        if(n_neighbours) mesh->n_vert_neighbours_[vid] = n_neighbours;
    }

    free(edges);
}

void ntcd_mesh_terminate(ntcd_mesh* mesh){
    free(mesh->vertices_);
    free(mesh->n_vert_neighbours_);
    for(unsigned int vid = 0; vid < mesh->n_vertices_; ++vid){
        free(mesh->vert_neighbours_[vid]);
    }
    free(mesh->vert_neighbours_);
}

//Hull
static void ntcd__support_hull(double* support_point, const void* shape, const double* dir){
    ntcd_hull hull = *(const ntcd_hull*)shape;
    ntcd_support sa = *(const ntcd_support*)hull.ca_;
    ntcd_support sb = *(const ntcd_support*)hull.cb_;

    double inv_dir_a[3];
    ntcd__quat_vec3_rotate(inv_dir_a, hull.inv_rot_a_, dir);
    sa(support_point, hull.ca_, inv_dir_a);
    ntcd__quat_vec3_rotate(support_point, hull.ta_.rot, support_point);
    ntcd__vec3_fmadd(support_point, hull.ta_.size, support_point, hull.ta_.pos);

    double inv_dir_b[3], temp_point[3];
    ntcd__quat_vec3_rotate(inv_dir_b, hull.inv_rot_b_, dir);
    sb(temp_point, hull.cb_, inv_dir_b);
    ntcd__quat_vec3_rotate(temp_point, hull.tb_.rot, temp_point);
    ntcd__vec3_fmadd(temp_point, hull.tb_.size, temp_point, hull.tb_.pos);

    if(ntcd__vec3_dot(temp_point, dir) > ntcd__vec3_dot(support_point, dir)){
        memcpy(support_point, temp_point, 3 * sizeof(*support_point));
    }
}

void ntcd_hull_initialize(ntcd_hull* hull, const ntcd_transform* ta, const void* ca, const ntcd_transform* tb, const void* cb){
    hull->support = ntcd__support_hull;
    hull->ca_     = ca;
    hull->cb_     = cb;
    hull->ta_     = *ta;
    hull->tb_     = *tb;

    ntcd__quat_inverse(hull->inv_rot_a_, hull->ta_.rot);
    ntcd__quat_inverse(hull->inv_rot_b_, hull->tb_.rot);
}

#endif //NTCD_IMPLEMENTATION
