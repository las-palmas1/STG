#ifndef SPECTRAL_H
#define SPECTRAL_H

#include "common.h"


typedef struct STG_SpectralData_s
{
    STG_int i_cnt;
    STG_int j_cnt;
    STG_int k_cnt;

    STG_int num_modes;
    STG_float * u_abs;

    STG_float * phi;
    STG_float * psi;
    STG_float * alpha;
    STG_float * theta;
    // multiplier before time term in expression for velocity
    STG_float * omega;

} STG_SpectralData;


STG_SHARED_LIB_API void STG_compute_Spectral_data(
    STG_InitData init_data, STG_int num_modes, STG_SpectralData * data
);

STG_SHARED_LIB_API void STG_free_Spectral_data(STG_SpectralData * data);


// num_ts - number of current time step (start with 0)
STG_SHARED_LIB_API void STG_compute_Spectral_moment_field(
        STG_InitData init_data, STG_SpectralData data, STG_float ts, STG_int num_ts, STG_VelMomField * mom_field
);


// num_ts_tot - total number of time steps
STG_SHARED_LIB_API void STG_compute_Spectral_node_hist(
	STG_InitData init_data, STG_SpectralData data, STG_float ts, STG_int num_ts_tot,
    STG_VelNodeHist * node_hist, STG_int i, STG_int j, STG_int k
);


STG_SHARED_LIB_API void STG_compute_Spectral_pulsations(
	STG_float * k1, STG_float * k2, STG_float * k3,
    STG_float * sigma1, STG_float * sigma2, STG_float * sigma3,
    STG_float * psi, STG_float * omega, STG_float * u_abs,
	STG_float x, STG_float y, STG_float z,
    STG_int num_modes, STG_float time_scale, STG_float time,
	STG_float * u_p, STG_float * v_p, STG_float * w_p
);


// These params (phi, psi, alpha, theta and omega) are shared for all nodes
STG_SHARED_LIB_API void STG_compute_Spectral_random_angles_and_phase(
        STG_int num_modes, STG_float * phi, STG_float * psi,
        STG_float * alpha, STG_float * theta, STG_float * omega
);

STG_SHARED_LIB_API void STG_compute_Spectral_modes_params(
        STG_int num_modes, STG_float * k_arr, STG_float * phi,
        STG_float * alpha, STG_float * theta,
        STG_float * k1, STG_float * k2, STG_float * k3,
        STG_float * sigma1, STG_float * sigma2, STG_float * sigma3
);



#endif // !SPECTRAL_H
