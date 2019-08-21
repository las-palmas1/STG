#include "precompiled.h"
#include "Davidson.h"
#include <stdlib.h>
#include <stdio.h>


void STG_free_Davidson_stat_data(STG_DavidsonData_Stationary * data)
{
	free(data->a11);
	free(data->a12);
	free(data->a13);
	free(data->a21);
	free(data->a22);
	free(data->a23);
	free(data->a31);
	free(data->a32);
	free(data->a33);
	free(data->c1);
	free(data->c2);
	free(data->c3);
	free(data->k_arr);
	free(data->energy);
	free(data->u_abs);
    free(data->a);
    free(data->b);
}

void STG_free_Davidson_trans_data(STG_DavidsonData_Transient * data)
{
	free(data->phi);
	free(data->psi);
	free(data->alpha);
	free(data->theta);
	free(data->k1);
	free(data->k2);
	free(data->k3);
	free(data->sigma1);
	free(data->sigma2);
	free(data->sigma3);
	free(data->u_p_prev);
	free(data->v_p_prev);
	free(data->w_p_prev);
}

void STG_compute_Davidson_matrix_data(
	STG_float re_uu, STG_float re_vv, STG_float re_ww,
	STG_float re_uv, STG_float re_uw, STG_float re_vw,
	STG_float * c1, STG_float * c2, STG_float * c3,
	STG_float * a11, STG_float * a12, STG_float * a13,
	STG_float * a21, STG_float * a22, STG_float * a23,
	STG_float * a31, STG_float * a32, STG_float * a33
) 
{
	STG_float trace = re_uu + re_vv + re_ww;
	STG_float eig_vals[3], eig_vec1[3], eig_vec2[3], eig_vec3[3];
	compute_eig(
		re_uu, re_vv, re_ww, re_uv, re_uw, re_vw,
		eig_vals, eig_vec1, eig_vec2, eig_vec3
	);
	*c1 = sqrt(eig_vals[0] * 3 / trace);
	*c2 = sqrt(eig_vals[1] * 3 / trace);
	*c3 = sqrt(eig_vals[2] * 3 / trace);
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

void STG_compute_Davidson_spectrum(
	STG_float delta_min, STG_int num_modes, STG_float re_uu, STG_float re_vv, STG_float re_ww, STG_float ls_i,
	STG_float dissip_rate, STG_float visc, STG_float * energy, STG_float * k_arr, STG_float * u_abs
)
{
	STG_float alpha = 1.483;
    STG_float pi = 2 * asin(1);
	STG_float k_t = 0.5 * (re_uu + re_vv + re_ww);
	STG_float k_e = alpha * 9 * pi / (55 * ls_i);
	STG_float k_min = k_e / 2;
	STG_float k_max = 2 * pi / (2 * delta_min);
	STG_float delta_k = (k_max - k_min) / (num_modes - 1);
	
	for (STG_int i = 0; i < num_modes; i++) 
	{
		k_arr[i] = k_min + delta_k * i;
	}
    STG_float k_eta = powf(dissip_rate, 0.25) / powf(visc, 0.75);
    STG_float u_rms = sqrt(2. / 3. * k_t);

	for (STG_int i = 0; i < num_modes; i++)
	{
        energy[i] = alpha * u_rms*u_rms / k_e * pow(k_arr[i] / k_e, 4) * exp(-2 * pow(k_arr[i] / k_eta, 2)) /
            pow(1 + k_arr[i] * k_arr[i] / (k_e * k_e), 17. / 6.);
		u_abs[i] = sqrt(energy[i] * delta_k);
	}
}

void STG_alloc_Davidson_trans_data(STG_InitData init_data, STG_int num_modes, STG_int num_ts_tot,
	STG_DavidsonData_Transient * data)
{
	data->phi = (STG_float *)malloc(sizeof(STG_float) * num_modes * (num_ts_tot + 1));
	data->psi = (STG_float *)malloc(sizeof(STG_float) * num_modes * (num_ts_tot + 1));
	data->alpha = (STG_float *)malloc(sizeof(STG_float) * num_modes * (num_ts_tot + 1));
	data->theta = (STG_float *)malloc(sizeof(STG_float) * num_modes * (num_ts_tot + 1));
	data->k1 = (STG_float *)malloc(sizeof(STG_float) * num_modes * (num_ts_tot + 1));
	data->k2 = (STG_float *)malloc(sizeof(STG_float) * num_modes * (num_ts_tot + 1));
	data->k3 = (STG_float *)malloc(sizeof(STG_float) * num_modes * (num_ts_tot + 1));
	data->sigma1 = (STG_float *)malloc(sizeof(STG_float) * num_modes * (num_ts_tot + 1));
	data->sigma2 = (STG_float *)malloc(sizeof(STG_float) * num_modes * (num_ts_tot + 1));
	data->sigma3 = (STG_float *)malloc(sizeof(STG_float) * num_modes * (num_ts_tot + 1));

	STG_int mesh_size = init_data.i_cnt * init_data.j_cnt * init_data.k_cnt;
	data->u_p_prev = (STG_float *)malloc(sizeof(STG_float) * mesh_size);
	data->v_p_prev = (STG_float *)malloc(sizeof(STG_float) * mesh_size);
	data->w_p_prev = (STG_float *)malloc(sizeof(STG_float) * mesh_size);

	for (STG_int i = 0; i < mesh_size; i++)
	{
		data->u_p_prev[i] = 0.;
		data->v_p_prev[i] = 0.;
		data->w_p_prev[i] = 0.;
	}

}

static void STG_compute_delta_min(STG_InitData init_data, STG_float * delta_min)
{
	// Считаем, что сетка равномерная по всем трем координатным направлениям
	STG_int i0 = GET_INDEX(0, 0, 0, init_data.i_cnt, init_data.j_cnt, init_data.k_cnt);
	STG_int ii_next = GET_INDEX(1, 0, 0, init_data.i_cnt, init_data.j_cnt, init_data.k_cnt);
	STG_int ij_next = GET_INDEX(0, 1, 0, init_data.i_cnt, init_data.j_cnt, init_data.k_cnt);
	STG_int ik_next = GET_INDEX(0, 0, 1, init_data.i_cnt, init_data.j_cnt, init_data.k_cnt);
	
    STG_float dx = fabs(init_data.mesh.x[ii_next] - init_data.mesh.x[i0]);
	STG_float dy = fabs(init_data.mesh.y[ij_next] - init_data.mesh.y[i0]);
	STG_float dz = fabs(init_data.mesh.z[ik_next] - init_data.mesh.z[i0]);
	*delta_min = min(dx, min(dy, dz));
}


void STG_compute_Davidson_auto_coef(STG_float ts, STG_float ts_i, STG_float *a, STG_float *b)
{
    *a = exp(-ts / ts_i);
    *b = sqrt(1 - (*a)*(*a));
}


void STG_compute_Davidson_stat_data(
    STG_InitData init_data, STG_int num_modes, STG_float dissip_rate, STG_float visc, STG_float ts, 
	STG_DavidsonData_Stationary * data
)
{
	data->i_cnt = init_data.i_cnt;
	data->j_cnt = init_data.j_cnt;
	data->k_cnt = init_data.k_cnt;
	data->num_modes = num_modes;
	STG_int mesh_size = init_data.i_cnt * init_data.j_cnt * init_data.k_cnt;

    data->energy = (STG_float *)malloc(sizeof(STG_float) * num_modes * mesh_size);
    data->k_arr = (STG_float *)malloc(sizeof(STG_float) * num_modes * mesh_size);
    data->u_abs = (STG_float *)malloc(sizeof(STG_float) * num_modes * mesh_size);

	data->a11 = (STG_float *)malloc(sizeof(STG_float) * mesh_size);
	data->a12 = (STG_float *)malloc(sizeof(STG_float) * mesh_size);
	data->a13 = (STG_float *)malloc(sizeof(STG_float) * mesh_size);
	data->a21 = (STG_float *)malloc(sizeof(STG_float) * mesh_size);
	data->a22 = (STG_float *)malloc(sizeof(STG_float) * mesh_size);
	data->a23 = (STG_float *)malloc(sizeof(STG_float) * mesh_size);
	data->a31 = (STG_float *)malloc(sizeof(STG_float) * mesh_size);
	data->a32 = (STG_float *)malloc(sizeof(STG_float) * mesh_size);
	data->a33 = (STG_float *)malloc(sizeof(STG_float) * mesh_size);
	data->c1 = (STG_float *)malloc(sizeof(STG_float) * mesh_size);
	data->c2 = (STG_float *)malloc(sizeof(STG_float) * mesh_size);
	data->c3 = (STG_float *)malloc(sizeof(STG_float) * mesh_size);

    data->a = (STG_float *)malloc(sizeof(STG_float) * mesh_size);
    data->b = (STG_float *)malloc(sizeof(STG_float) * mesh_size);

    STG_float delta_min;
    STG_compute_delta_min(init_data, &delta_min);

    for (STG_int i = 0; i < mesh_size; i++)
    {
        STG_compute_Davidson_matrix_data(
                    init_data.re.re_uu[i], init_data.re.re_vv[i], init_data.re.re_ww[i],
                    init_data.re.re_uv[i], init_data.re.re_uw[i], init_data.re.re_vw[i],
                    &(data->c1[i]), &(data->c2[i]), &(data->c3[i]),
                    &(data->a11[i]), &(data->a12[i]), &(data->a13[i]),
                    &(data->a21[i]), &(data->a22[i]), &(data->a23[i]),
                    &(data->a31[i]), &(data->a32[i]), &(data->a33[i])
                    );
        STG_compute_Davidson_auto_coef((STG_float)ts, (STG_float)init_data.scales.ts_i[i], &(data->a[i]), &(data->b[i]));
        STG_int i_spec = i * num_modes;
        STG_compute_Davidson_spectrum(
                    delta_min, num_modes, init_data.re.re_uu[i], init_data.re.re_vv[i], init_data.re.re_ww[i],
                    (STG_float)init_data.scales.ls_i[i], dissip_rate, visc,
                    &(data->energy[i_spec]), &(data->k_arr[i_spec]), &(data->u_abs[i_spec])
                    );

    }
}


void STG_compute_Davidson_random_data(
        STG_int num_modes, STG_float * k_arr, STG_float * phi, STG_float * psi, STG_float * alpha, STG_float * theta,
        STG_float * k1, STG_float * k2, STG_float * k3,
        STG_float * sigma1, STG_float * sigma2, STG_float * sigma3)
{
    STG_float pi = 2 * asin(1.);
    get_uniform_ref(0., 2 * pi, num_modes, phi);
    get_uniform_ref(0., 2 * pi, num_modes, psi);
    get_uniform_ref(0., 2 * pi, num_modes, alpha);
    get_trigon_ref(num_modes, theta);
    for (STG_int i = 0; i < num_modes; i++)
    {
        k1[i] = sin(theta[i]) * cos(phi[i]) * k_arr[i];
        k2[i] = sin(theta[i]) * sin(phi[i]) * k_arr[i];
        k3[i] = cos(theta[i]) * k_arr[i];
        sigma1[i] = cos(phi[i]) * cos(theta[i]) * cos(alpha[i]) - sin(phi[i]) * sin(alpha[i]);
        sigma2[i] = sin(phi[i]) * cos(theta[i]) * cos(alpha[i]) - cos(phi[i]) * sin(alpha[i]);
        sigma3[i] = -sin(theta[i]) * cos(alpha[i]);
    }

}

void STG_compute_Davidson_trans_data(
        STG_DavidsonData_Stationary stat_data, STG_int num_ts_tot, STG_DavidsonData_Transient * data)
{
	data->num_modes = stat_data.num_modes;
	data->num_ts = num_ts_tot;
	for (STG_int i = 0; i < num_ts_tot + 1; i++)
	{
		STG_int i_mode = i * stat_data.num_modes;
		STG_compute_Davidson_random_data(
			stat_data.num_modes, stat_data.k_arr, &(data->phi[i_mode]), &(data->psi[i_mode]), 
			&(data->alpha[i_mode]), &(data->theta[i_mode]),
			&(data->k1[i_mode]), &(data->k2[i_mode]), &(data->k3[i_mode]), 
			&(data->sigma1[i_mode]), &(data->sigma2[i_mode]), &(data->sigma3[i_mode]));
	}
}


void STG_compute_Davidson_pulsations(
	STG_float * k1, STG_float * k2, STG_float * k3,
	STG_float * sigma1, STG_float * sigma2, STG_float * sigma3, STG_float * psi, STG_float * u_abs,
	STG_float c1, STG_float c2, STG_float c3,
	STG_float a11, STG_float a12, STG_float a13,
	STG_float a21, STG_float a22, STG_float a23,
	STG_float a31, STG_float a32, STG_float a33,
	STG_float x, STG_float y, STG_float z,
    STG_float a, STG_float b, STG_int num_modes,
	STG_float u_p_prev, STG_float v_p_prev, STG_float w_p_prev,
	STG_float * u_p, STG_float * v_p, STG_float * w_p
)
{
	STG_float v1 = 0.;
	STG_float v2 = 0.;
	STG_float v3 = 0.;
	for (STG_int i = 0; i < num_modes; i++)
	{
		v1 += 2 * u_abs[i] * cos(k1[i] * x + k2[i] * y + k3[i] * z + psi[i]) * sigma1[i];
		v2 += 2 * u_abs[i] * cos(k1[i] * x + k2[i] * y + k3[i] * z + psi[i]) * sigma2[i];
		v3 += 2 * u_abs[i] * cos(k1[i] * x + k2[i] * y + k3[i] * z + psi[i]) * sigma3[i];
	}
	STG_float w1 = c1 * v1;
	STG_float w2 = c2 * v2;
	STG_float w3 = c3 * v3;
	STG_float v1_prime = a11 * w1 + a12 * w2 + a13 * w3;
	STG_float v2_prime = a21 * w1 + a22 * w2 + a23 * w3;
	STG_float v3_prime = a31 * w1 + a32 * w2 + a33 * w3;
	*u_p = a * u_p_prev + b * v1_prime;
	*v_p = a * v_p_prev + b * v2_prime;
	*w_p = a * w_p_prev + b * v3_prime;
}


void STG_compute_Davidson_moment_field(STG_InitData init_data,
	STG_DavidsonData_Stationary stat_data, STG_DavidsonData_Transient * trans_data, STG_float ts, STG_int num_ts,
	STG_VelMomField * mom_field)
{
	mom_field->time = ts * num_ts;
	mom_field->i_cnt = init_data.i_cnt;
	mom_field->j_cnt = init_data.j_cnt;
	mom_field->k_cnt = init_data.k_cnt;
	STG_int mesh_size = init_data.i_cnt * init_data.j_cnt * init_data.k_cnt;
	mom_field->u_p = (STG_float*)malloc(sizeof(STG_float) * mesh_size);
	mom_field->v_p = (STG_float*)malloc(sizeof(STG_float) * mesh_size);
	mom_field->w_p = (STG_float*)malloc(sizeof(STG_float) * mesh_size);

	STG_float a;
	STG_float b;
	for (STG_int i = 0; i < mesh_size; i++)
	{
		if (num_ts == 0)
		{
			a = 0.;
			b = 1.;
		}
		else
		{
			a = stat_data.a[i];
			b = stat_data.b[i];
		}
		trans_data->u_p_prev[i] = 0.;
		trans_data->v_p_prev[i] = 0.;
		trans_data->w_p_prev[i] = 0.;
		STG_int i_mode = num_ts * stat_data.num_modes;
		STG_int i_spec = i * stat_data.num_modes;
		STG_compute_Davidson_pulsations(
			&(trans_data->k1[i_mode]), &(trans_data->k2[i_mode]), &(trans_data->k3[i_mode]),
			&(trans_data->sigma1[i_mode]), &(trans_data->sigma2[i_mode]), &(trans_data->sigma3[i_mode]), 
			&(trans_data->psi[i_mode]), &(stat_data.u_abs[i_spec]),
			stat_data.c1[i], stat_data.c2[i], stat_data.c3[i],
			stat_data.a11[i], stat_data.a12[i], stat_data.a13[i],
			stat_data.a21[i], stat_data.a22[i], stat_data.a23[i],
			stat_data.a31[i], stat_data.a32[i], stat_data.a33[i],
			init_data.mesh.x[i], init_data.mesh.y[i], init_data.mesh.z[i],
			a, b, stat_data.num_modes,
			trans_data->u_p_prev[i], trans_data->v_p_prev[i], trans_data->w_p_prev[i],
			&(mom_field->u_p[i]), &(mom_field->v_p[i]), &(mom_field->w_p[i])
		);
		trans_data->u_p_prev[i] = mom_field->u_p[i];
		trans_data->v_p_prev[i] = mom_field->v_p[i];
		trans_data->w_p_prev[i] = mom_field->w_p[i];
	}
}


void STG_compute_Davidson_node_hist(STG_InitData init_data,
	STG_DavidsonData_Stationary stat_data, STG_float ts, STG_int num_ts_tot,
	STG_DavidsonData_Transient * trans_data, STG_VelNodeHist * node_hist, STG_int i, STG_int j, STG_int k)
{
	node_hist->i = i;
	node_hist->j = j;
	node_hist->k = k;
	STG_int num = GET_INDEX(i, j, k, init_data.i_cnt, init_data.j_cnt, init_data.k_cnt);
	node_hist->num_ts = num_ts_tot;
	node_hist->time = (STG_float *)malloc(sizeof(STG_float) * (num_ts_tot + 1));
	node_hist->u_p = (STG_float *)malloc(sizeof(STG_float) * (num_ts_tot + 1));
	node_hist->v_p = (STG_float *)malloc(sizeof(STG_float) * (num_ts_tot + 1));
	node_hist->w_p = (STG_float *)malloc(sizeof(STG_float) * (num_ts_tot + 1));

	STG_float time;
	STG_float a, b;
	for (STG_int it = 0; it < num_ts_tot + 1; it++)
	{
		if (it == 0) 
		{
			trans_data->u_p_prev[num] = 0.;
			trans_data->v_p_prev[num] = 0.;
			trans_data->w_p_prev[num] = 0.;
			a = 0.;
			b = 1.;
		}
		else
		{
			trans_data->u_p_prev[num] = node_hist->u_p[it - 1];
			trans_data->v_p_prev[num] = node_hist->v_p[it - 1];
			trans_data->w_p_prev[num] = node_hist->w_p[it - 1];
			a = stat_data.a[num];
			b = stat_data.b[num];
		}
		time = ts * it;
		node_hist->time[it] = time;
		STG_int i_mode = it * stat_data.num_modes;
		STG_int num_spec = num * stat_data.num_modes;
		STG_compute_Davidson_pulsations(
			&(trans_data->k1[i_mode]), &(trans_data->k2[i_mode]), &(trans_data->k3[i_mode]),
			&(trans_data->sigma1[i_mode]), &(trans_data->sigma2[i_mode]), &(trans_data->sigma3[i_mode]), 
			&(trans_data->psi[i_mode]), &(stat_data.u_abs[num_spec]), 
			stat_data.c1[num], stat_data.c2[num], stat_data.c3[num],
			stat_data.a11[num], stat_data.a12[num], stat_data.a13[num],
			stat_data.a21[num], stat_data.a22[num], stat_data.a23[num],
			stat_data.a31[num], stat_data.a32[num], stat_data.a33[num],
			init_data.mesh.x[num], init_data.mesh.y[num], init_data.mesh.z[num],
			stat_data.a[num], stat_data.b[num], stat_data.num_modes,
			trans_data->u_p_prev[num], trans_data->v_p_prev[num], trans_data->w_p_prev[num],
			&(node_hist->u_p[it]), &(node_hist->v_p[it]), &(node_hist->w_p[it])
		);
	}
}