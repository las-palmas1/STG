#include "precompiled.h"
#include "Smirnov.h"
#include <stdlib.h>
#include <stdio.h>


void STG_compute_Smirnov_matrix_data_homo(
	STG_float re_uu, STG_float re_vv, STG_float re_ww,
	STG_float re_uv, STG_float re_uw, STG_float re_vw,
	STG_float * c1, STG_float * c2, STG_float * c3,
	STG_float * a11, STG_float * a12, STG_float * a13,
	STG_float * a21, STG_float * a22, STG_float * a23,
	STG_float * a31, STG_float * a32, STG_float * a33
)
{
	STG_float eig_vals[3], eig_vec1[3], eig_vec2[3], eig_vec3[3];
	compute_eig(
		re_uu, re_vv, re_ww, re_uv, re_uw, re_vw,
		eig_vals, eig_vec1, eig_vec2, eig_vec3
	);
	*c1 = eig_vals[0];
	*c2 = eig_vals[1];
	*c3 = eig_vals[2];
	*a11 = eig_vec1[0];
	*a21 = eig_vec1[1];
	*a31 = eig_vec1[2];
	*a12 = eig_vec2[0];
	*a22 = eig_vec2[1];
	*a32 = eig_vec2[2];
	*a13 = eig_vec3[0];
	*a23 = eig_vec3[1];
	*a33 = eig_vec3[2];
}

void STG_compute_Smirnov_random_data(
	STG_int num_modes, STG_float * omega,
	STG_float * k1, STG_float * k2, STG_float * k3,
	STG_float * zeta1, STG_float * zeta2, STG_float * zeta3,
	STG_float * xi1, STG_float * xi2, STG_float * xi3,
	STG_float * p1, STG_float * p2, STG_float * p3,
	STG_float * q1, STG_float * q2, STG_float * q3
)
{
	get_normal_ref(0, 0.5, num_modes, k1);
	get_normal_ref(0, 0.5, num_modes, k2);
	get_normal_ref(0, 0.5, num_modes, k3);
	get_normal_ref(0, 1, num_modes, zeta1);
	get_normal_ref(0, 1, num_modes, zeta2);
	get_normal_ref(0, 1, num_modes, zeta3);
	get_normal_ref(0, 1, num_modes, xi1);
	get_normal_ref(0, 1, num_modes, xi2);
	get_normal_ref(0, 1, num_modes, xi3);
	get_normal_ref(0, 1, num_modes, omega);
	for (STG_int i = 0; i < num_modes; i++)
	{
		p1[i] = zeta2[i] * k3[i] - zeta3[i] * k2[i];
		p2[i] = zeta3[i] * k1[i] - zeta1[i] * k3[i];
		p3[i] = zeta1[i] * k2[i] - zeta2[i] * k1[i];

		q1[i] = xi2[i] * k3[i] - xi3[i] * k2[i];
		q2[i] = xi3[i] * k1[i] - xi1[i] * k3[i];
		q3[i] = xi1[i] * k2[i] - xi2[i] * k1[i];
	}
}

void STG_compute_Smirnov_data_time_indep(STG_InitData init_data, STG_int num_modes, STG_SmirnovData_TimeIndep * data)
{
	data->num_modes = num_modes;
	STG_int is = init_data.i_cnt;
	STG_int js = init_data.j_cnt;
	STG_int ks = init_data.k_cnt;

	data->i_cnt = is;
	data->j_cnt = js;
	data->k_cnt = ks;

	data->c1 = (STG_float*)malloc(sizeof(STG_float) * is * js * ks);
	data->c2 = (STG_float*)malloc(sizeof(STG_float) * is * js * ks);
	data->c3 = (STG_float*)malloc(sizeof(STG_float) * is * js * ks);
	data->a11 = (STG_float*)malloc(sizeof(STG_float) * is *	 js * ks);
	data->a12 = (STG_float*)malloc(sizeof(STG_float) * is * js * ks);
	data->a13 = (STG_float*)malloc(sizeof(STG_float) * is * js * ks);
	data->a21 = (STG_float*)malloc(sizeof(STG_float) * is * js * ks);
	data->a22 = (STG_float*)malloc(sizeof(STG_float) * is * js * ks);
	data->a23 = (STG_float*)malloc(sizeof(STG_float) * is * js * ks);
	data->a31 = (STG_float*)malloc(sizeof(STG_float) * is * js * ks);
	data->a32 = (STG_float*)malloc(sizeof(STG_float) * is * js * ks);
	data->a33 = (STG_float*)malloc(sizeof(STG_float) * is * js * ks);

	data->k1 = get_normal(0, 0.5, num_modes);
	data->k2 = get_normal(0, 0.5, num_modes);
	data->k3 = get_normal(0, 0.5, num_modes);
	data->zeta1 = get_normal(0, 1, num_modes);
	data->zeta2 = get_normal(0, 1, num_modes);
	data->zeta3 = get_normal(0, 1, num_modes);
	data->xi1 = get_normal(0, 1, num_modes);
	data->xi2 = get_normal(0, 1, num_modes);
	data->xi3 = get_normal(0, 1, num_modes);
	data->omega = get_normal(0, 1, num_modes);
	data->p1 = (STG_float*)malloc(sizeof(STG_float) * num_modes);
	data->p2 = (STG_float*)malloc(sizeof(STG_float) * num_modes);
	data->p3 = (STG_float*)malloc(sizeof(STG_float) * num_modes);
	data->q1 = (STG_float*)malloc(sizeof(STG_float) * num_modes);
	data->q2 = (STG_float*)malloc(sizeof(STG_float) * num_modes);
	data->q3 = (STG_float*)malloc(sizeof(STG_float) * num_modes);

	for (STG_int i = 0; i < num_modes; i++)
	{
		data->p1[i] = data->zeta2[i] * data->k3[i] - data->zeta3[i] * data->k2[i];
		data->p2[i] = data->zeta3[i] * data->k1[i] - data->zeta1[i] * data->k3[i];
		data->p3[i] = data->zeta1[i] * data->k2[i] - data->zeta2[i] * data->k1[i];

		data->q1[i] = data->xi2[i] * data->k3[i] - data->xi3[i] * data->k2[i];
		data->q2[i] = data->xi3[i] * data->k1[i] - data->xi1[i] * data->k3[i];
		data->q3[i] = data->xi1[i] * data->k2[i] - data->xi2[i] * data->k1[i];
	}


	for (STG_int i = 0; i < init_data.i_cnt; i++)
	{
		for (STG_int j = 0; j < init_data.j_cnt; j++)
		{
			for (STG_int k = 0; k < init_data.k_cnt; k++)
			{
				STG_float eig_vals[3], eig_vec1[3], eig_vec2[3], eig_vec3[3];
				STG_int index = GET_INDEX(i, j, k, is, js, ks);
				compute_eig(
					init_data.re.re_uu[index], init_data.re.re_vv[index], init_data.re.re_ww[index],
					init_data.re.re_uv[index], init_data.re.re_uw[index], init_data.re.re_vw[index],
					eig_vals, eig_vec1, eig_vec2, eig_vec3
				);
				data->c1[index] = eig_vals[0];
				data->c2[index] = eig_vals[1];
				data->c3[index] = eig_vals[2];
				data->a11[index] = eig_vec1[0];
				data->a21[index] = eig_vec1[1];
				data->a31[index] = eig_vec1[2];
				data->a12[index] = eig_vec2[0];
				data->a22[index] = eig_vec2[1];
				data->a32[index] = eig_vec2[2];
				data->a13[index] = eig_vec3[0];
				data->a23[index] = eig_vec3[1];
				data->a33[index] = eig_vec3[2];
			}
		}
	}
}


void STG_compute_Smirnov_data_time_dep(STG_float ts, STG_int num_ts, STG_SmirnovData_TimeDep * data)
{
	data->num_ts = num_ts;
	data->ts = ts;
}


void STG_free_Smirnov_data_time_indep(STG_SmirnovData_TimeIndep * data)
{
	free(data->c1);
	free(data->c2);
	free(data->c3);
	free(data->a11);
	free(data->a12);
	free(data->a13);
	free(data->a21);
	free(data->a22);
	free(data->a23);
	free(data->a31);
	free(data->a32);
	free(data->a33);

	free(data->k1);
	free(data->k2);
	free(data->k3);
	free(data->zeta1);
	free(data->zeta2);
	free(data->zeta3);
	free(data->xi1);
	free(data->xi2);
	free(data->xi3);
	free(data->p1);
	free(data->p2);
	free(data->p3);
	free(data->q1);
	free(data->q2);
	free(data->q3);
	free(data->omega);
}


void STG_compute_Smirnov_pulsations(
	STG_float * k1, STG_float * k2, STG_float * k3, 
	STG_float * p1, STG_float * p2, STG_float * p3,
	STG_float * q1, STG_float * q2, STG_float * q3, STG_float * omega, 
	STG_float c1, STG_float c2, STG_float c3,
	STG_float a11, STG_float a12, STG_float a13,
	STG_float a21, STG_float a22, STG_float a23,
	STG_float a31, STG_float a32, STG_float a33,
	STG_float x, STG_float y, STG_float z,
	STG_float length_scale, STG_float time_scale, STG_int num_modes, STG_float time,
	STG_float * u, STG_float * v, STG_float * w
)
{	
	STG_float mode1_sum = 0;
	STG_float mode2_sum = 0;
	STG_float mode3_sum = 0;
	STG_float k1_p, k2_p, k3_p;

	for (STG_int i = 0; i < num_modes; i++)
	{
		k1_p = k1[i] * (length_scale / time_scale) / c1;
		k2_p = k2[i] * (length_scale / time_scale) / c2;
		k3_p = k3[i] * (length_scale / time_scale) / c3;
		mode1_sum += MODE(1, i);
		mode2_sum += MODE(2, i);
		mode3_sum += MODE(3, i);
	}
 	STG_float v1 = sqrt((STG_float)2 / num_modes) * mode1_sum;
	STG_float v2 = sqrt((STG_float)2 / num_modes) * mode2_sum;
	STG_float v3 = sqrt((STG_float)2 / num_modes) * mode3_sum;
	STG_float w1 = c1 * v1;
	STG_float w2 = c2 * v2;
	STG_float w3 = c3 * v3;
	*u = a11 * w1 + a12 * w2 + a13 * w3;
	*v = a21 * w1 + a22 * w2 + a23 * w3;
	*w = a31 * w1 + a32 * w2 + a33 * w3;
}


void STG_compute_Smirnov(
	STG_InitData init_data, STG_SmirnovData_TimeIndep data_tind, 
	STG_SmirnovData_TimeDep data_tdep, STG_OutData * out_data
)
{
	STG_int is = data_tind.i_cnt;
	STG_int js = data_tind.j_cnt;
	STG_int ks = data_tind.k_cnt;
	out_data->i_cnt = is;
	out_data->j_cnt = js;
	out_data->k_cnt = ks;
	out_data->u_p = (STG_float**)malloc(sizeof(STG_float*) * (data_tdep.num_ts + 1));
	out_data->v_p = (STG_float**)malloc(sizeof(STG_float*) * (data_tdep.num_ts + 1));
	out_data->w_p = (STG_float**)malloc(sizeof(STG_float*) * (data_tdep.num_ts + 1));
	out_data->time = (STG_float*)malloc(sizeof(STG_float) * (data_tdep.num_ts + 1));
	out_data->num_ts = data_tdep.num_ts;

	STG_int num = is * js * ks;
	for (STG_int it = 0; it < data_tdep.num_ts + 1; it++)
	{
		STG_float time = data_tdep.ts * it;
		out_data->time[it] = time;
		out_data->u_p[it] = (STG_float*)malloc(sizeof(STG_float) * num);
		out_data->v_p[it] = (STG_float*)malloc(sizeof(STG_float) * num);
		out_data->w_p[it] = (STG_float*)malloc(sizeof(STG_float) * num);
		for (STG_int i = 0; i < num; i++)
		{
			STG_compute_Smirnov_pulsations(
				data_tind.k1, data_tind.k2, data_tind.k3, 
				data_tind.p1, data_tind.p2, data_tind.p3, 
				data_tind.q1, data_tind.q2, data_tind.q3, data_tind.omega,
				data_tind.c1[i], data_tind.c2[i], data_tind.c3[i], 
				data_tind.a11[i], data_tind.a12[i], data_tind.a13[i],
				data_tind.a21[i], data_tind.a22[i], data_tind.a23[i],
				data_tind.a31[i], data_tind.a32[i], data_tind.a33[i], 
				init_data.mesh.x[i], init_data.mesh.y[i], init_data.mesh.z[i],
				init_data.scales.length_scale[i], init_data.scales.time_scale[i], 
				data_tind.num_modes, time,  
				&(out_data->u_p[it][i]), &(out_data->v_p[it][i]), &(out_data->w_p[it][i])
			);
		}
	}
}


void STG_compute_Smirnov_node(
	STG_InitData init_data, STG_SmirnovData_TimeIndep data_tind, 
	STG_SmirnovData_TimeDep data_tdep, STG_OutDataNode * out_data, STG_int i, STG_int j, STG_int k
)
{
	out_data->i = i;
	out_data->j = j;
	out_data->k = k;
	STG_int num = GET_INDEX(i, j, k, init_data.i_cnt, init_data.j_cnt, init_data.k_cnt);
	out_data->time = (STG_float*)malloc(sizeof(STG_float) * (data_tdep.num_ts + 1));
	out_data->num_ts = data_tdep.num_ts;
	out_data->u_p = (STG_float*)malloc(sizeof(STG_float) * (data_tdep.num_ts + 1));
	out_data->v_p = (STG_float*)malloc(sizeof(STG_float) * (data_tdep.num_ts + 1));
	out_data->w_p = (STG_float*)malloc(sizeof(STG_float) * (data_tdep.num_ts + 1));
	for (STG_int it = 0; it < data_tdep.num_ts + 1; it++)
	{
		STG_float time = data_tdep.ts * it;
		out_data->time[it] = time;
		STG_compute_Smirnov_pulsations(
			data_tind.k1, data_tind.k2, data_tind.k3, 
			data_tind.p1, data_tind.p2, data_tind.p3, 
			data_tind.q1, data_tind.q2, data_tind.q3, data_tind.omega,
			data_tind.c1[num], data_tind.c2[num], data_tind.c3[num], 
			data_tind.a11[num], data_tind.a12[num], data_tind.a13[num],
			data_tind.a21[num], data_tind.a22[num], data_tind.a23[num],
			data_tind.a31[num], data_tind.a32[num], data_tind.a33[num], 
			init_data.mesh.x[num], init_data.mesh.y[num], init_data.mesh.z[num],
			init_data.scales.length_scale[num], init_data.scales.time_scale[num], 
			data_tind.num_modes, out_data->time[it],
			&(out_data->u_p[it]), &(out_data->v_p[it]), &(out_data->w_p[it])
		);
	}
}