#ifndef SEM_H
#define SEM_H

#include "common.h"
#include "assert.h"

#define FACE_NUM 6

typedef struct Vector_s
{
    STG_float x, y, z;
} Vector;


typedef struct Plane_s
{
    STG_float a, b, c, d;
} Plane;


typedef struct Limits_s
{
    STG_float x_min, x_max;
    STG_float y_min, y_max;
    STG_float z_min, z_max;
} Limits;


void compute_line_plane_intersection(Vector p0, Vector vec, Plane pl, Vector * p_inter);

STG_float get_scalar_prod(Vector v1, Vector v2);

Limits * get_in_planes_lims(Limits vol_lim, Vector *eddies_pos, STG_int num_eddies, Vector vel);


void compute_init_eddies_params_fort(
        STG_float * x_e, STG_float * y_e, STG_float * z_e,
        STG_float * eps1, STG_float * eps2, STG_float * eps3,
        STG_int num_eddies,
        STG_float x_min, STG_float x_max,
        STG_float y_min, STG_float y_max,
        STG_float z_min, STG_float z_max
);


void compute_new_eddies_params_fort(
        STG_float * x_e_cur, STG_float * y_e_cur, STG_float * z_e_cur,
        STG_float * eps_x_cur, STG_float * eps_y_cur, STG_float * eps_z_cur,
        STG_int num_eddies,
        STG_float x_min, STG_float x_max,
        STG_float y_min, STG_float y_max,
        STG_float z_min, STG_float z_max,
        STG_float u0, STG_float v0, STG_float w0,
        STG_float * x_min_in, STG_float * x_max_in,
        STG_float * y_min_in, STG_float * y_max_in,
        STG_float * z_min_in, STG_float * z_max_in,
        STG_float * x_e_new, STG_float * y_e_new, STG_float * z_e_new,
        STG_float * eps_x_new, STG_float * eps_y_new, STG_float * eps_z_new
);

void compute_init_eddies_params(Vector * eddies_pos, Vector * eddies_int, STG_int num_eddies, Limits vol_lims);

void compute_new_eddies_params(
        Vector * eddies_pos_cur, Vector * eddies_int_cur, STG_int num_eddies, Limits vol_lims, Vector vel,
        Limits * in_planes_lims, Vector * eddies_pos_new, Vector * eddies_int_new
);




#endif // !SEM_H
