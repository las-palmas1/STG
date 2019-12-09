#include "precompiled.h"
#include "Spectral.h"
#include <stdlib.h>
#include <stdio.h>


void STG_free_Spectral_data(STG_SpectralData * data)
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

    free(data->u_abs);

    free(data->phi);
    free(data->psi);
    free(data->alpha);
    free(data->theta);
    free(data->omega);
}



void STG_compute_Spectral_matrix_data(
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
    if ((re_uu == 0) && (re_vv == 0) && (re_ww == 0))
    {
        trace = eig_vals[0] + eig_vals[1] + eig_vals[2];
    }
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



void STG_compute_Spectral_random_angles_and_phase(
        STG_int num_modes, STG_float * phi, STG_float * psi,
        STG_float * alpha, STG_float * theta, STG_float * omega
)
{
    STG_float pi = 2 * asin(1.);

    get_uniform_ref(0., 2 * pi, num_modes, phi);
    get_uniform_ref(0., 2 * pi, num_modes, psi);
    get_uniform_ref(0., 2 * pi, num_modes, alpha);
    get_trigon_ref(num_modes, theta);
    get_normal_ref(2 * pi, 2 * pi, num_modes, omega);
}



void STG_compute_Spectral_modes_params(
        STG_int num_modes, STG_float * k_arr, STG_float * phi,
        STG_float * alpha, STG_float * theta,
        STG_float * k1, STG_float * k2, STG_float * k3,
        STG_float * sigma1, STG_float * sigma2, STG_float * sigma3
)
{
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



void STG_compute_Spectral_data(
    STG_InitData init_data, STG_int num_modes, STG_SpectralData * data
)
{
    data->i_cnt = init_data.i_cnt;
    data->j_cnt = init_data.j_cnt;
    data->k_cnt = init_data.k_cnt;

    data->num_modes = num_modes;
    STG_int mesh_size = init_data.i_cnt * init_data.j_cnt * init_data.k_cnt;

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

    data->phi = (STG_float *)malloc(sizeof(STG_float) * num_modes);
    data->psi = (STG_float *)malloc(sizeof(STG_float) * num_modes);
    data->alpha = (STG_float *)malloc(sizeof(STG_float) * num_modes);
    data->theta = (STG_float *)malloc(sizeof(STG_float) * num_modes);
    data->omega = (STG_float *)malloc(sizeof(STG_float) * num_modes);

    for (STG_int i_node = 0; i_node < mesh_size; i_node++)
    {
        for (STG_int i_mode = 0; i_mode < num_modes; i_mode++)
        {
            STG_int i = i_node * num_modes + i_mode;
            STG_float delta_k;

            if (i_mode < num_modes - 1){
                delta_k = init_data.spectrum.k_arr[i + 1] - init_data.spectrum.k_arr[i];
            } else {
                delta_k = init_data.spectrum.k_arr[i] - init_data.spectrum.k_arr[i - 1];
            }
            data->u_abs[i] = sqrt(init_data.spectrum.energy[i] * delta_k);
        }
    }

    for (STG_int i_node = 0; i_node < mesh_size; i_node++)
    {
        STG_compute_Spectral_matrix_data(
            init_data.re.re_uu[i_node], init_data.re.re_vv[i_node], init_data.re.re_ww[i_node],
            init_data.re.re_uv[i_node], init_data.re.re_uw[i_node], init_data.re.re_vw[i_node],
            &(data->c1[i_node]), &(data->c2[i_node]), &(data->c3[i_node]),
            &(data->a11[i_node]), &(data->a12[i_node]), &(data->a13[i_node]),
            &(data->a21[i_node]), &(data->a22[i_node]), &(data->a23[i_node]),
            &(data->a31[i_node]), &(data->a32[i_node]), &(data->a33[i_node])
        );
    }
    STG_compute_Spectral_random_angles_and_phase(
        num_modes, data->phi, data->psi, data->alpha, data->theta, data->omega
    );
}



void STG_compute_Spectral_pulsations(
    STG_float * k1, STG_float * k2, STG_float * k3,
    STG_float * sigma1, STG_float * sigma2, STG_float * sigma3,
    STG_float * psi, STG_float * omega, STG_float * u_abs,
    STG_float c1, STG_float c2, STG_float c3,
    STG_float a11, STG_float a12, STG_float a13,
    STG_float a21, STG_float a22, STG_float a23,
    STG_float a31, STG_float a32, STG_float a33,
    STG_float x, STG_float y, STG_float z,
    STG_int num_modes, STG_float time_scale, STG_float time,
    STG_float * u_p, STG_float * v_p, STG_float * w_p
)
{
    STG_float v1 = 0.;
    STG_float v2 = 0.;
    STG_float v3 = 0.;
    for (STG_int i = 0; i < num_modes; i++)
    {
        v1 += 2 * u_abs[i] * cos(k1[i] * x + k2[i] * y + k3[i] * z + psi[i] + omega[i] * time / time_scale) * sigma1[i];
        v2 += 2 * u_abs[i] * cos(k1[i] * x + k2[i] * y + k3[i] * z + psi[i] + omega[i] * time / time_scale) * sigma2[i];
        v3 += 2 * u_abs[i] * cos(k1[i] * x + k2[i] * y + k3[i] * z + psi[i] + omega[i] * time / time_scale) * sigma3[i];
    }
    STG_float w1 = c1 * v1;
    STG_float w2 = c2 * v2;
    STG_float w3 = c3 * v3;
    *u_p = a11 * w1 + a12 * w2 + a13 * w3;
    *v_p = a21 * w1 + a22 * w2 + a23 * w3;
    *w_p = a31 * w1 + a32 * w2 + a33 * w3;
}



void STG_compute_Spectral_moment_field(
    STG_InitData init_data, STG_SpectralData data, STG_float ts, STG_int num_ts, STG_VelMomField * mom_field
)
{
    mom_field->time = ts * num_ts;
    mom_field->i_cnt = init_data.i_cnt;
    mom_field->j_cnt = init_data.j_cnt;
    mom_field->k_cnt = init_data.k_cnt;
    STG_int mesh_size = init_data.i_cnt * init_data.j_cnt * init_data.k_cnt;
    mom_field->u_p = (STG_float*)malloc(sizeof(STG_float) * mesh_size);
    mom_field->v_p = (STG_float*)malloc(sizeof(STG_float) * mesh_size);
    mom_field->w_p = (STG_float*)malloc(sizeof(STG_float) * mesh_size);

    STG_float * k1 = (STG_float *)malloc(sizeof(STG_float) * data.num_modes);
    STG_float * k2 = (STG_float *)malloc(sizeof(STG_float) * data.num_modes);
    STG_float * k3 = (STG_float *)malloc(sizeof(STG_float) * data.num_modes);
    STG_float * sigma1 = (STG_float *)malloc(sizeof(STG_float) * data.num_modes);
    STG_float * sigma2 = (STG_float *)malloc(sizeof(STG_float) * data.num_modes);
    STG_float * sigma3 = (STG_float *)malloc(sizeof(STG_float) * data.num_modes);

    for (STG_int i_node = 0; i_node < mesh_size; i_node++)
    {
        STG_int i_spec = data.num_modes * i_node;
        STG_compute_Spectral_modes_params(
            data.num_modes, &(init_data.spectrum.k_arr[i_spec]), data.phi, data.alpha,
            data.theta, k1, k2, k3, sigma1, sigma2, sigma3
        );

        STG_compute_Spectral_pulsations(
            k1, k2, k3, sigma1, sigma2, sigma3,
            data.psi, data.omega, &(data.u_abs[i_spec]),
            data.c1[i_node], data.c2[i_node], data.c3[i_node],
            data.a11[i_node], data.a12[i_node], data.a13[i_node],
            data.a21[i_node], data.a22[i_node], data.a23[i_node],
            data.a31[i_node], data.a32[i_node], data.a33[i_node],
            init_data.mesh.x[i_node], init_data.mesh.y[i_node], init_data.mesh.z[i_node],
            data.num_modes, init_data.scales.ts_i[i_node], mom_field->time,
            &(mom_field->u_p[i_node]), &(mom_field->v_p[i_node]), &(mom_field->w_p[i_node])
        );
    }
    free(k1);
    free(k2);
    free(k3);
    free(sigma1);
    free(sigma2);
    free(sigma3);
}



void STG_compute_Spectral_node_hist(STG_InitData init_data,
    STG_SpectralData data, STG_float ts, STG_int num_ts_tot,
    STG_VelNodeHist * node_hist, STG_int i, STG_int j, STG_int k
)
{
    node_hist->i = i;
    node_hist->j = j;
    node_hist->k = k;
    STG_int i_node = GET_INDEX(i, j, k, init_data.i_cnt, init_data.j_cnt, init_data.k_cnt);
    STG_int i_spec = i_node * data.num_modes;

    node_hist->num_ts = num_ts_tot;
    node_hist->time = (STG_float *)malloc(sizeof(STG_float) * (num_ts_tot + 1));
    node_hist->u_p = (STG_float *)malloc(sizeof(STG_float) * (num_ts_tot + 1));
    node_hist->v_p = (STG_float *)malloc(sizeof(STG_float) * (num_ts_tot + 1));
    node_hist->w_p = (STG_float *)malloc(sizeof(STG_float) * (num_ts_tot + 1));


    STG_float * k1 = (STG_float *)malloc(sizeof(STG_float) * data.num_modes);
    STG_float * k2 = (STG_float *)malloc(sizeof(STG_float) * data.num_modes);
    STG_float * k3 = (STG_float *)malloc(sizeof(STG_float) * data.num_modes);
    STG_float * sigma1 = (STG_float *)malloc(sizeof(STG_float) * data.num_modes);
    STG_float * sigma2 = (STG_float *)malloc(sizeof(STG_float) * data.num_modes);
    STG_float * sigma3 = (STG_float *)malloc(sizeof(STG_float) * data.num_modes);

    STG_compute_Spectral_modes_params(
        data.num_modes, &(init_data.spectrum.k_arr[i_spec]), data.phi, data.alpha,
        data.theta, k1, k2, k3, sigma1, sigma2, sigma3
    );

    STG_float time;
    for (STG_int it = 0; it < num_ts_tot + 1; it++)
    {
        time = ts * it;
        node_hist->time[it] = time;
        STG_compute_Spectral_pulsations(
            k1, k2, k3, sigma1, sigma2, sigma3,
            data.psi, data.omega, &(data.u_abs[i_spec]),
            data.c1[i_node], data.c2[i_node], data.c3[i_node],
            data.a11[i_node], data.a12[i_node], data.a13[i_node],
            data.a21[i_node], data.a22[i_node], data.a23[i_node],
            data.a31[i_node], data.a32[i_node], data.a33[i_node],
            init_data.mesh.x[i_node], init_data.mesh.y[i_node], init_data.mesh.z[i_node],
            data.num_modes, init_data.scales.ts_i[i_node], time,
            &(node_hist->u_p[it]), &(node_hist->v_p[it]), &(node_hist->w_p[it])
        );
    }
    free(k1);
    free(k2);
    free(k3);
    free(sigma1);
    free(sigma2);
    free(sigma3);
}





