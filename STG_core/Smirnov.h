#ifndef SMIRNOV_H
#define SMIRNOV_H

#include "common.h"

#define MODE(iaxe, imode)  p##iaxe[imode] * cos((k1_p * x + k2_p * y + k3_p * z) / length_scale + omega[imode] * time / time_scale) + \
						   q##iaxe[imode] * sin((k1_p * x + k2_p * y + k3_p * z) / length_scale + omega[imode] * time / time_scale)


// —труктура хран€ща€ вспомогательные данные дл€ метода —мирнова дл€ случа€ неоднородной турбулентности.
typedef struct STG_SmirnovData_TimeIndep_s
{
	STG_int i_cnt;
	STG_int j_cnt;
	STG_int k_cnt;
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
	STG_float * k1;
	STG_float * k2;
	STG_float * k3;
	STG_float * zeta1;
	STG_float * zeta2;
	STG_float * zeta3;
	STG_float * xi1;
	STG_float * xi2;
	STG_float * xi3;
	STG_float * omega;
	STG_float * p1;
	STG_float * p2;
	STG_float * p3;
	STG_float * q1;
	STG_float * q2;
	STG_float * q3;

} STG_SmirnovData_TimeIndep;


typedef struct STG_SmirnovData_TimeDep_s
{
	STG_float ts;
	STG_int num_ts;

} STG_SmirnovData_TimeDep;

// Ёта функци€ будет использоватьс€ в решателе дл€ расчета параметров однородной турбулентности
void STG_compute_Smirnov_matrix_data_homo(
	STG_float re_uu, STG_float re_vv, STG_float re_ww,
	STG_float re_uv, STG_float re_uw, STG_float re_vw,
	STG_float * c1, STG_float * c2, STG_float * c3,
	STG_float * a11, STG_float * a12, STG_float * a13,
	STG_float * a21, STG_float * a22, STG_float * a23,
	STG_float * a31, STG_float * a32, STG_float * a33
);

// Ѕудет использоватьс€ в решателе дл€ расчета случайных величин
void STG_compute_Smirnov_random_data(
	STG_int num_modes, STG_float * omega,
	STG_float * k1, STG_float * k2, STG_float * k3,
	STG_float * zeta1, STG_float * zeta2, STG_float * zeta3,
	STG_float * xi1, STG_float * xi2, STG_float * xi3,
	STG_float * p1, STG_float * p2, STG_float * p3,
	STG_float * q1, STG_float * q2, STG_float * q3
);

STG_SHARED_LIB_API void STG_compute_Smirnov_data_time_indep(
	STG_InitData init_data, STG_int num_modes, STG_SmirnovData_TimeIndep * data
);


STG_SHARED_LIB_API void STG_compute_Smirnov_data_time_dep(STG_float ts, STG_int num_ts, STG_SmirnovData_TimeDep * data);

STG_SHARED_LIB_API void STG_free_Smirnov_data_time_indep(STG_SmirnovData_TimeIndep * data);

// ƒанна€ функци€ будет использовать в решателе дл€ расчета пульсаций
STG_SHARED_LIB_API void STG_compute_Smirnov_pulsations(
	STG_float * k1, STG_float * k2, STG_float * k3, 
	STG_float * p1, STG_float * p2, STG_float * p3,
	STG_float * q1, STG_float * q2, STG_float * q3, 
	STG_float * omega, STG_float c1, STG_float c2, STG_float c3,
	STG_float a11, STG_float a12, STG_float a13,
	STG_float a21, STG_float a22, STG_float a23,
	STG_float a31, STG_float a32, STG_float a33,
	STG_float x, STG_float y, STG_float z,
	STG_float length_scale, STG_float time_scale, 
	STG_int num_modes, STG_float time, 
	STG_float * u, STG_float * v, STG_float * w
);

// Ќужна дл€ априорных тестов в питоне, рассчитывает трехмерное поле
STG_SHARED_LIB_API void STG_compute_Smirnov(
	STG_InitData init_data, STG_SmirnovData_TimeIndep data_tind,
	STG_SmirnovData_TimeDep data_tdep, STG_OutData * out_data
);

// Ќужна дл€ априорных тестов в питоне, рассчитывает полусации в точке
STG_SHARED_LIB_API void STG_compute_Smirnov_node(
	STG_InitData init_data, STG_SmirnovData_TimeIndep data_tind, STG_SmirnovData_TimeDep data_tdep, 
	STG_OutDataNode * out_data, STG_int i, STG_int j, STG_int k
);

#endif // !SMIRNOV_H
