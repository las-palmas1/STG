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

STG_float get_scalar_prod(Vector v1, Vector v2)
{
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
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
}


void compute_init_eddies_params(Vector * eddies_pos, Vector * eddies_int, STG_int num_eddies, Limits vol_lims)
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








