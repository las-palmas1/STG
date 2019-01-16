#ifndef SMIRNOV_H
#define SMIRNOV_H

#include "precompiled.h"
#include "common.h"


#define MODE(iaxe, imode)  p##iaxe[imode] * cos((k1_p * x + k2_p * y + k3_p * z) / length_scale + omega[imode] * time / time_scale) + \
						   q##iaxe[imode] * sin((k1_p * x + k2_p * y + k3_p * z) / length_scale + omega[imode] * time / time_scale)

typedef struct SmirnovData_s
{
	STG_float ts;
	STG_int num_ts;
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

} SmirnovData;

STG_SHARED_LIB_API void compute_Smirnov_data(InitData init_data, STG_int num_modes, STG_float ts, STG_int num_ts, SmirnovData * data);

STG_SHARED_LIB_API void free_Smirnov_data(SmirnovData * data);

STG_SHARED_LIB_API void compute_Smirnov_pulsations(
	STG_float * k1, STG_float * k2, STG_float * k3, STG_float * p1, STG_float * p2, STG_float * p3,
	STG_float * q1, STG_float * q2, STG_float * q3, STG_float * omega, STG_float c1, STG_float c2, STG_float c3,
	STG_float a11, STG_float a12, STG_float a13,
	STG_float a21, STG_float a22, STG_float a23,
	STG_float a31, STG_float a32, STG_float a33,
	STG_float x, STG_float y, STG_float z,
	STG_float length_scale, STG_float time_scale, STG_int num_modes, STG_float time, 
	STG_float * u, STG_float * v, STG_float * w
);

STG_SHARED_LIB_API void compute_Smirnov_field(InitData init_data, SmirnovData data, OutData * out_data);

#endif // !SMIRNOV_H
