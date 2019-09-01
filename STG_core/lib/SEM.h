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


typedef struct STG_SEMData_Stationary_s
{

    STG_int i_cnt;
    STG_int j_cnt;
    STG_int k_cnt;

    STG_float * a11;
    STG_float * a12;
    STG_float * a13;
    STG_float * a21;
    STG_float * a22;
    STG_float * a23;
    STG_float * a31;
    STG_float * a32;
    STG_float * a33;

    STG_int num_eddies;
    Limits * in_planes_lims;
    // Вопрос: что задавать - скорость переноса вихрей или временной масштаб
    Vector eddies_vel;
} STG_SEMData_Stationary;

typedef struct STG_SEMData_Transient_s
{
    STG_int num_ts;
    STG_int num_eddies;
    Vector * eddies_int;
    Vector * eddies_pos;

} STG_SEMData_Transient;

void compute_line_plane_intersection(Vector p0, Vector vec, Plane pl, Vector * p_inter);

void compute_limits(STG_InitData init_data, Limits * lims);

void STG_compute_SEM_matrix_data(
        STG_float re_uu, STG_float re_vv, STG_float re_ww,
        STG_float re_uv, STG_float re_uw, STG_float re_vw,
        STG_float * a11, STG_float * a12, STG_float * a13,
        STG_float * a21, STG_float * a22, STG_float * a23,
        STG_float * a31, STG_float * a32, STG_float * a33
);

STG_float get_scalar_prod(Vector v1, Vector v2);

Limits * get_in_planes_lims(Limits vol_lim, Vector *eddies_pos, STG_int num_eddies, Vector vel);


void STG_compute_SEM_init_eddies_params_fort(
        STG_float * x_e, STG_float * y_e, STG_float * z_e,
        STG_float * eps1, STG_float * eps2, STG_float * eps3,
        STG_int num_eddies,
        STG_float x_min, STG_float x_max,
        STG_float y_min, STG_float y_max,
        STG_float z_min, STG_float z_max
);


void STG_compute_SEM_new_eddies_params_fort(
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

void STG_compute_SEM_init_eddies_params(Vector * eddies_pos, Vector * eddies_int, STG_int num_eddies, Limits vol_lims);

void STG_compute_SEM_new_eddies_params(
        Vector * eddies_pos_cur, Vector * eddies_int_cur, STG_int num_eddies, Limits vol_lims, Vector vel,
        STG_float dt, Limits * in_planes_lims, Vector * eddies_pos_new, Vector * eddies_int_new
);




#endif // !SEM_H
