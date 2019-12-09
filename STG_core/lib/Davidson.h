#ifndef DAVIDSON_H
#define DAVIDSON_H

#include "common.h"

#define MODE(iaxe, imode)  p##iaxe[imode] * cos((k1_p * x + k2_p * y + k3_p * z) / length_scale + omega[imode] * time / time_scale) + \
                           q##iaxe[imode] * sin((k1_p * x + k2_p * y + k3_p * z) / length_scale + omega[imode] * time / time_scale)


typedef struct STG_DavidsonData_Stationary_s
{
    STG_int i_cnt;
    STG_int j_cnt;
    STG_int k_cnt;
    STG_float *a;
    STG_float *b;

    STG_float * c1;
    STG_float * c2;
    STG_float * c3;
    STG_float * a11;
    STG_float * a12;
    STG_float * a13;
    STG_float * a21;
    STG_float * a22;
    STG_float * a23;
    STG_float * a31;
    STG_float * a32;
    STG_float * a33;

    STG_int num_modes;
	// ????????? ???????
	STG_float * energy;
	STG_float * k_arr;
	STG_float * u_abs;

} STG_DavidsonData_Stationary;


typedef struct STG_DavidsonData_Transient_s
{
	STG_int num_ts;
	STG_int num_modes;

	STG_float * phi;
	STG_float * psi;
	STG_float * alpha;
	STG_float * theta;

	STG_float * u_p_prev;
	STG_float * v_p_prev;
	STG_float * w_p_prev;

} STG_DavidsonData_Transient;


STG_SHARED_LIB_API void STG_compute_Davidson_stat_data(
        STG_InitData init_data, STG_int num_modes, STG_float dissip_rate,
        STG_float visc, STG_float ts, STG_DavidsonData_Stationary * data
);

STG_SHARED_LIB_API void STG_free_Davidson_stat_data(STG_DavidsonData_Stationary * data);

// num_ts_tot - total number of time steps
STG_SHARED_LIB_API void STG_alloc_Davidson_trans_data(
        STG_InitData init_data, STG_int num_modes, STG_int num_ts_tot,
        STG_DavidsonData_Transient * data
);


STG_SHARED_LIB_API void STG_compute_Davidson_trans_data(
	STG_DavidsonData_Stationary stat_data, STG_int num_ts_tot, STG_DavidsonData_Transient * data
);


STG_SHARED_LIB_API void STG_free_Davidson_trans_data(STG_DavidsonData_Transient * data);

// num_ts - number of current time step (start with 0)
STG_SHARED_LIB_API void STG_compute_Davidson_moment_field(
    STG_InitData init_data, STG_DavidsonData_Stationary stat_data, STG_DavidsonData_Transient * trans_data,
    STG_float ts, STG_int num_ts, STG_VelMomField * mom_field
);


// num_ts_tot - total number of time steps
STG_SHARED_LIB_API void STG_compute_Davidson_node_hist(
    STG_InitData init_data, STG_DavidsonData_Stationary stat_data,
    STG_float ts, STG_int num_ts_tot, STG_DavidsonData_Transient * trans_data,
    STG_VelNodeHist * node_hist, STG_int i, STG_int j, STG_int k
);


STG_SHARED_LIB_API void STG_compute_Davidson_pulsations(
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
);

STG_SHARED_LIB_API void STG_compute_Davidson_auto_coef(STG_float ts, STG_float ts_i, STG_float *a, STG_float *b);


STG_SHARED_LIB_API void STG_compute_Davidson_matrix_data(
	STG_float re_uu, STG_float re_vv, STG_float re_ww,
	STG_float re_uv, STG_float re_uw, STG_float re_vw,
	STG_float * c1, STG_float * c2, STG_float * c3,
	STG_float * a11, STG_float * a12, STG_float * a13,
	STG_float * a21, STG_float * a22, STG_float * a23,
	STG_float * a31, STG_float * a32, STG_float * a33
);


STG_SHARED_LIB_API void STG_compute_Davidson_spectrum(
	STG_float delta_min, STG_int num_modes, STG_float re_uu, STG_float re_vv, STG_float re_ww, STG_float ls_i,
    STG_float dissip_rate, STG_float visc, STG_float * energy, STG_float * k_arr, STG_float * u_abs
);


void STG_compute_Davidson_random_angles_and_phase(
	STG_int num_modes, STG_float * phi, STG_float * psi, STG_float * alpha, STG_float * theta
);


void STG_compute_Davidson_modes_params(
	STG_int num_modes, STG_float * k_arr, STG_float * phi, STG_float * psi, STG_float * alpha, STG_float * theta,
	STG_float * k1, STG_float * k2, STG_float * k3,
	STG_float * sigma1, STG_float * sigma2, STG_float * sigma3
);

STG_SHARED_LIB_API void STG_compute_Davidson_random_data(
        STG_int num_modes, STG_float * k_arr, STG_float * phi, STG_float * psi,
        STG_float * alpha, STG_float * theta,
        STG_float * k1, STG_float * k2, STG_float * k3,
        STG_float * sigma1, STG_float * sigma2, STG_float * sigma3
);
#endif // !DAVIDSON_H
