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
        STG_float ts, Limits * in_planes_lims, Vector * eddies_pos_new, Vector * eddies_int_new
)
{
    for (STG_int i = 0; i < num_eddies; i++)
    {
        eddies_pos_new[i].x = eddies_pos_cur[i].x + ts * vel.x;
        eddies_pos_new[i].y = eddies_pos_cur[i].y + ts * vel.y;
        eddies_pos_new[i].z = eddies_pos_cur[i].z + ts * vel.z;

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


void STG_compute_SEM_stat_data(
	STG_InitData init_data, STG_int num_eddies, Vector eddies_vel, STG_SEMData_Stationary * stat_data
)
{
	stat_data->i_cnt = init_data.i_cnt;
	stat_data->j_cnt = init_data.j_cnt;
	stat_data->k_cnt = init_data.k_cnt;
	STG_int mesh_size = init_data.i_cnt * init_data.j_cnt * init_data.k_cnt;

	stat_data->a11 = (STG_float *)malloc(sizeof(STG_float) * mesh_size);
	stat_data->a12 = (STG_float *)malloc(sizeof(STG_float) * mesh_size);
	stat_data->a13 = (STG_float *)malloc(sizeof(STG_float) * mesh_size);
	stat_data->a21 = (STG_float *)malloc(sizeof(STG_float) * mesh_size);
	stat_data->a22 = (STG_float *)malloc(sizeof(STG_float) * mesh_size);
	stat_data->a23 = (STG_float *)malloc(sizeof(STG_float) * mesh_size);
	stat_data->a31 = (STG_float *)malloc(sizeof(STG_float) * mesh_size);
	stat_data->a32 = (STG_float *)malloc(sizeof(STG_float) * mesh_size);
	stat_data->a33 = (STG_float *)malloc(sizeof(STG_float) * mesh_size);

	for (STG_int i = 0; i < mesh_size; i++)
	{
		STG_compute_SEM_matrix_data(
			init_data.re.re_uu[i], init_data.re.re_vv[i], init_data.re.re_ww[i],
			init_data.re.re_uv[i], init_data.re.re_uw[i], init_data.re.re_vw[i],
			&(stat_data->a11[i]), &(stat_data->a12[i]), &(stat_data->a13[i]),
			&(stat_data->a21[i]), &(stat_data->a21[i]), &(stat_data->a21[i]),
			&(stat_data->a31[i]), &(stat_data->a31[i]), &(stat_data->a31[i])
		);
	}

	stat_data->num_eddies = num_eddies;
	stat_data->eddies_vel = eddies_vel;
	stat_data->eddies_pos_init = (Vector*)malloc(sizeof(Vector) * num_eddies);
	stat_data->eddies_int_init = (Vector*)malloc(sizeof(Vector) * num_eddies);
	compute_limits(init_data, &(stat_data->vol_lims));

	STG_compute_SEM_init_eddies_params(stat_data->eddies_pos_init, stat_data->eddies_int_init, 
		num_eddies, stat_data->vol_lims);

	stat_data->in_planes_lims = get_in_planes_lims(stat_data->vol_lims, stat_data->eddies_pos_init, num_eddies, eddies_vel);
}

void STG_free_SEM_stat_data(STG_SEMData_Stationary * stat_data)
{
	free(stat_data->a11);
	free(stat_data->a12);
	free(stat_data->a13);
	free(stat_data->a21);
	free(stat_data->a22);
	free(stat_data->a23);
	free(stat_data->a31);
	free(stat_data->a32);
	free(stat_data->a33);

	free(stat_data->eddies_int_init);
	free(stat_data->eddies_pos_init);
	free(stat_data->in_planes_lims);
}


void STG_compute_SEM_trans_data(
	STG_SEMData_Stationary stat_data, STG_float ts, STG_int num_ts, STG_SEMData_Transient * trans_data
)
{
	trans_data->num_eddies = stat_data.num_eddies;
	trans_data->num_ts = num_ts;
	trans_data->ts = ts;
	trans_data->eddies_int = (Vector *)malloc(sizeof(Vector) * num_ts * stat_data.num_eddies);
	trans_data->eddies_int = (Vector *)malloc(sizeof(Vector) * num_ts * stat_data.num_eddies);

	for (STG_int i = 0; i < stat_data.num_eddies; i++)
	{
		trans_data->eddies_int[i] = stat_data.eddies_int_init[i];
		trans_data->eddies_pos[i] = stat_data.eddies_pos_init[i];
	}
	for (STG_int it = 0; it < num_ts - 1; it++)
	{
		for (STG_int i = 0; i < stat_data.num_eddies; i++)
		{
			STG_int i_cur = it * stat_data.num_eddies + i;
			STG_int i_next = (it + 1) * stat_data.num_eddies + i;
			STG_compute_SEM_new_eddies_params(
				&(trans_data->eddies_pos[i_cur]), &(trans_data->eddies_int[i_cur]), stat_data.num_eddies,
				stat_data.vol_lims, stat_data.eddies_vel, ts, stat_data.in_planes_lims,
				&(trans_data->eddies_pos[i_next]), &(trans_data->eddies_int[i_next])
			);
		}
	}
}

void STG_free_SEM_trans_data(STG_SEMData_Transient * trans_data)
{
	free(trans_data->eddies_int);
	free(trans_data->eddies_pos);
}


static STG_float form_func(STG_float x, STG_float x_e, STG_float ls)
{
	STG_float res;
	if ((abs(x - x_e) / ls) < 1)
	{
		res = sqrt(1.5) * (1 - abs((x - x_e) / ls));
	}
	else
	{
		res = 0.;
	}
	return res;
}


void STG_compute_SEM_pulsations(
        STG_float * x_e, STG_float * y_e, STG_float * z_e,
        STG_float * eps_x, STG_float * eps_y, STG_float * eps_z,
        STG_int num_eddies, STG_float x, STG_float y, STG_float z, STG_float volume,
        STG_float ls_ux, STG_float ls_uy, STG_float ls_uz,
        STG_float ls_vx, STG_float ls_vy, STG_float ls_vz,
        STG_float ls_wx, STG_float ls_wy, STG_float ls_wz,
        STG_float a11, STG_float a12, STG_float a13,
        STG_float a21, STG_float a22, STG_float a23,
        STG_float a31, STG_float a32, STG_float a33,
        STG_float * u, STG_float * v, STG_float * w
)
{
    STG_float f_u = 0.;
    STG_float f_v = 0.;
    STG_float f_w = 0.;

	for (STG_int i = 0; i < num_eddies; i++)
	{
        STG_float f_ux = form_func(x, x_e[i], ls_ux);
        STG_float f_uy = form_func(y, y_e[i], ls_uy);
        STG_float f_uz = form_func(z, z_e[i], ls_uz);
        STG_float f_vx = form_func(x, x_e[i], ls_vx);
        STG_float f_vy = form_func(y, y_e[i], ls_vy);
        STG_float f_vz = form_func(z, z_e[i], ls_vz);
        STG_float f_wx = form_func(x, x_e[i], ls_wx);
        STG_float f_wy = form_func(y, y_e[i], ls_wy);
        STG_float f_wz = form_func(z, z_e[i], ls_wz);
        f_u += eps_x[i] * f_ux * f_uy * f_uz;
        f_v += eps_y[i] * f_vx * f_vy * f_vz;
        f_w += eps_z[i] * f_wx * f_wy * f_wz;
	}
    // note: check if sqrt in the last case is necessary
    STG_float u1 = 1 / sqrt(num_eddies) * sqrt(volume) / sqrt(ls_ux * ls_uy * ls_uz) * f_u;
    STG_float v1 = 1 / sqrt(num_eddies) * sqrt(volume) / sqrt(ls_vx * ls_vy * ls_vz) * f_v;
    STG_float w1 = 1 / sqrt(num_eddies) * sqrt(volume) / sqrt(ls_wx * ls_wy * ls_wz) * f_w;
    *u = a11 * u1 + a12 * v1 + a13 * w1;
    *v = a21 * u1 + a22 * v1 + a23 * w1;
    *w = a31 * u1 + a32 * v1 + a33 * w1;
}


static void extract_arrays(Vector * vec, STG_int num, STG_float * x, STG_float * y, STG_float * z)
{
    x = (STG_float *)malloc(sizeof (STG_float) * num);
    y = (STG_float *)malloc(sizeof (STG_float) * num);
    z = (STG_float *)malloc(sizeof (STG_float) * num);
    for (STG_int i = 0; i < num; i++)
    {
        x[i] = vec[i].x;
        y[i] = vec[i].y;
        z[i] = vec[i].z;
    }
}

void STG_compute_SEM_moment_field(
        STG_SEMData_Stationary stat_data, STG_SEMData_Transient trans_data, STG_float ts,
        STG_int num_ts, STG_VelMomField * mom_field
)
{
    mom_field->i_cnt = stat_data.i_cnt;
    mom_field->j_cnt = stat_data.j_cnt;
    mom_field->k_cnt = stat_data.k_cnt;
    STG_int mesh_size = stat_data.i_cnt * stat_data.j_cnt * stat_data.k_cnt;
    mom_field->time = num_ts * ts;

    mom_field->u_p = (STG_float *)malloc(sizeof (STG_float) * mesh_size);
    mom_field->v_p = (STG_float *)malloc(sizeof (STG_float) * mesh_size);
    mom_field->w_p = (STG_float *)malloc(sizeof (STG_float) * mesh_size);

    STG_float * x_e, * y_e, * z_e;
    STG_float * eps_x, * eps_y, * eps_z;
    STG_int i_ts_start = num_ts * stat_data.num_eddies;

    extract_arrays(&(trans_data.eddies_pos[i_ts_start]), stat_data.num_eddies, x_e, y_e, z_e);
    extract_arrays(&(trans_data.eddies_int[i_ts_start]), stat_data.num_eddies, eps_x, eps_x, eps_z);

    for (STG_int i = 0; i < mesh_size; i++)
    {
//        STG_compute_SEM_pulsations()
    }

}










