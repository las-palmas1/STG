#include "precompiled.h"
#include "SEM.h"
#include <stdlib.h>
#include <stdio.h>


void compute_line_plane_intersection(Vector p0, Vector vec, Plane pl, Vector * p_inter)
{
    STG_float k = pl.a * vec.x + pl.b * vec.y + pl.c * vec.z;
    STG_float x_res = (p0.x * (pl.b * vec.y + pl.c * vec.z) - vec.x * (pl.b * p0.y + pl.c * p0.z + pl.d));
    STG_float y_res = (p0.y * (pl.a * vec.x + pl.c * vec.z) - vec.y * (pl.a * p0.x + pl.c * p0.z + pl.d));
    STG_float z_res = (p0.z * (pl.b * vec.y + pl.a * vec.x) - vec.z * (pl.b * p0.y + pl.a * p0.x + pl.d));
    p_inter->x = x_res / k;
    p_inter->y = y_res / k;
    p_inter->z = z_res / k;
}


void STG_compute_SEM_matrix_data(
        STG_float re_uu, STG_float re_vv, STG_float re_ww,
        STG_float re_uv, STG_float re_uw, STG_float re_vw,
        STG_float * a11, STG_float * a12, STG_float * a13,
        STG_float * a21, STG_float * a22, STG_float * a23,
        STG_float * a31, STG_float * a32, STG_float * a33
)
{
    *a11 = sqrt(re_uu);
    *a12 = 0.;
    *a13 = 0.;
    *a21 = re_uv / (*a11);
    *a22 = sqrt(re_vv - (*a21) * (*a21));
    *a23 = 0.;
    *a31 = re_uw / (*a11);
    *a32 = (re_vw - (*a21) * (*a31)) / (*a22);
    *a33 = sqrt(re_ww - (*a31)*(*a31) - (*a32)*(*a32));
}

STG_float get_scalar_prod(Vector v1, Vector v2)
{
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}


static STG_float find_max(STG_float * arr, STG_int num)
{
    STG_float res = arr[0];
    for (STG_int i = 1; i < num; i++)
    {
        if (arr[i] > arr[i - 1])
        {
            res = arr[i];
        }
    }
    return  res;
}

static STG_float find_min(STG_float * arr, STG_int num)
{
    STG_float res = arr[0];
    for (STG_int i = 1; i < num; i++)
    {
        if (arr[i] < arr[i - 1])
        {
            res = arr[i];
        }
    }
    return  res;
}

void compute_limits(STG_InitData init_data, Limits * lims)
{
    STG_float ls_x_max, ls_y_max, ls_z_max;
    STG_Scales scales = init_data.scales;
    STG_int mesh_size = init_data.i_cnt * init_data.j_cnt * init_data.k_cnt;

    STG_float ls_ux_max = find_max(scales.ls_ux, mesh_size);
    STG_float ls_uy_max = find_max(scales.ls_uy, mesh_size);
    STG_float ls_uz_max = find_max(scales.ls_uz, mesh_size);
    STG_float ls_vx_max = find_max(scales.ls_vx, mesh_size);
    STG_float ls_vy_max = find_max(scales.ls_vy, mesh_size);
    STG_float ls_vz_max = find_max(scales.ls_vz, mesh_size);
    STG_float ls_wx_max = find_max(scales.ls_wx, mesh_size);
    STG_float ls_wy_max = find_max(scales.ls_wy, mesh_size);
    STG_float ls_wz_max = find_max(scales.ls_wz, mesh_size);

    STG_float ls_x_max_arr[3] = { ls_ux_max, ls_vx_max, ls_wx_max };
    ls_x_max = find_max(ls_x_max_arr, 3);
    STG_float ls_y_max_arr[3] = { ls_uy_max, ls_vy_max, ls_wy_max };
    ls_y_max = find_max(ls_y_max_arr, 3);
    STG_float ls_z_max_arr[3] = { ls_uz_max, ls_vz_max, ls_wz_max };
    ls_z_max = find_max(ls_z_max_arr, 3);

    STG_float x_min = find_min(init_data.mesh.x, mesh_size) - ls_x_max;
    STG_float x_max = find_max(init_data.mesh.x, mesh_size) + ls_x_max;
    STG_float y_min = find_min(init_data.mesh.y, mesh_size) - ls_y_max;
    STG_float y_max = find_max(init_data.mesh.y, mesh_size) + ls_y_max;
    STG_float z_min = find_min(init_data.mesh.z, mesh_size) - ls_z_max;
    STG_float z_max = find_max(init_data.mesh.z, mesh_size) + ls_z_max;

    lims->x_max = x_max;
    lims->x_min = x_min;
    lims->y_max = y_max;
    lims->y_min = y_min;
    lims->z_max = z_max;
    lims->z_min = z_min;
}

Limits * get_in_planes_lims(Limits vol_lim, Vector *eddies_pos, STG_int num_eddies, Vector vel)
{
    Plane volume_faces[FACE_NUM] = {
        {1., 0., 0., -vol_lim.x_min}, {1., 0., 0., -vol_lim.x_max},
        {0., 1., 0., -vol_lim.y_min}, {0., 1., 0., -vol_lim.y_max},
        {0., 0., 1., -vol_lim.z_min}, {0., 0., 1., -vol_lim.z_max},
    };
    Limits face_limits[FACE_NUM] = {
        {.x_min = vol_lim.x_min, .x_max = vol_lim.x_min, .y_min = vol_lim.y_min, .y_max = vol_lim.y_max,
        .z_min = vol_lim.z_min, .z_max = vol_lim.z_max},

        {.x_min = vol_lim.x_max, .x_max = vol_lim.x_max, .y_min = vol_lim.y_min, .y_max = vol_lim.y_max,
        .z_min = vol_lim.z_min, .z_max = vol_lim.z_max},

        {.x_min = vol_lim.x_min, .x_max = vol_lim.x_max, .y_min = vol_lim.y_min, .y_max = vol_lim.y_min,
        .z_min = vol_lim.z_min, .z_max = vol_lim.z_max},

        {.x_min = vol_lim.x_min, .x_max = vol_lim.x_max, .y_min = vol_lim.y_max, .y_max = vol_lim.y_max,
        .z_min = vol_lim.z_min, .z_max = vol_lim.z_max},

        {.x_min = vol_lim.x_min, .x_max = vol_lim.x_max, .y_min = vol_lim.y_min, .y_max = vol_lim.y_max,
        .z_min = vol_lim.z_min, .z_max = vol_lim.z_min},

        {.x_min = vol_lim.x_min, .x_max = vol_lim.x_max, .y_min = vol_lim.y_min, .y_max = vol_lim.y_max,
        .z_min = vol_lim.z_max, .z_max = vol_lim.z_max}
    };

    Limits * in_planes_lims = (Limits * )malloc(sizeof (Limits) * num_eddies);
    for (STG_int i = 0; i < num_eddies; i++)
    {
        for (STG_int j = 0; j < FACE_NUM; j++)
        {
            Vector intersec;
            compute_line_plane_intersection(eddies_pos[i], vel, volume_faces[j], &intersec);
            if (isfinite(intersec.x) && isfinite(intersec.y) && isfinite(intersec.z))
            {
                if ((intersec.x >= vol_lim.x_min) && (intersec.x <= vol_lim.x_max) &&
                        (intersec.y >= vol_lim.y_min) && (intersec.y <= vol_lim.y_max) &&
                        (intersec.z >= vol_lim.z_min) && (intersec.z <= vol_lim.z_max))
                {
                    Vector inters_vec = {
                        .x = intersec.x - eddies_pos[i].x, .y = intersec.y - eddies_pos[i].y, .z = intersec.z - eddies_pos[i].z
                    };
                    STG_float s_prod = get_scalar_prod(inters_vec, vel);
                    if (s_prod < 0)
                    {
                        in_planes_lims[i] = face_limits[j];
                    }
                }
            }
        }
    }
    return in_planes_lims;
}


void STG_compute_SEM_init_eddies_params(Vector * eddies_pos, Vector * eddies_int, STG_int num_eddies, Limits vol_lims)
{
    for (STG_int i = 0; i < num_eddies; i++)
    {
        get_uniform_ref(vol_lims.x_min, vol_lims.x_max, 1, &(eddies_pos[i].x));
        get_uniform_ref(vol_lims.y_min, vol_lims.y_max, 1, &(eddies_pos[i].y));
        get_uniform_ref(vol_lims.z_min, vol_lims.z_max, 1, &(eddies_pos[i].z));
        get_normal_ref(0, 1, 1, &(eddies_int[i].x));
        get_normal_ref(0, 1, 1, &(eddies_int[i].y));
        get_normal_ref(0, 1, 1, &(eddies_int[i].z));
    }
}


void STG_compute_SEM_new_eddies_params(
        Vector * eddies_pos_cur, Vector * eddies_int_cur, STG_int num_eddies, Limits vol_lims, Vector vel,
        STG_float dt, Limits * in_planes_lims, Vector * eddies_pos_new, Vector * eddies_int_new
)
{
    for (STG_int i = 0; i < num_eddies; i++)
    {
        eddies_pos_new[i].x = eddies_pos_cur[i].x + dt * vel.x;
        eddies_pos_new[i].y = eddies_pos_cur[i].y + dt * vel.y;
        eddies_pos_new[i].z = eddies_pos_cur[i].z + dt * vel.z;

        if (eddies_pos_new[i].x < vol_lims.x_min || eddies_int_new[i].x > vol_lims.x_max ||
                eddies_pos_new[i].y < vol_lims.y_min || eddies_int_new[i].y > vol_lims.y_max ||
                eddies_pos_new[i].z < vol_lims.z_min || eddies_int_new[i].z > vol_lims.z_max)
        {
            get_uniform_ref(in_planes_lims[i].x_min, in_planes_lims[i].x_max, 1, &(eddies_pos_new[i].x));
            get_uniform_ref(in_planes_lims[i].y_min, in_planes_lims[i].y_max, 1, &(eddies_pos_new[i].y));
            get_uniform_ref(in_planes_lims[i].z_min, in_planes_lims[i].z_max, 1, &(eddies_pos_new[i].z));
            get_normal_ref(0, 1, 1, &(eddies_int_new[i].x));
            get_normal_ref(0, 1, 1, &(eddies_int_new[i].y));
            get_normal_ref(0, 1, 1, &(eddies_int_new[i].z));
        }
        else {
            eddies_int_new[i].x = eddies_int_cur[i].x;
            eddies_int_new[i].y = eddies_int_cur[i].y;
            eddies_int_new[i].z = eddies_int_cur[i].z;
        }
    }

}









