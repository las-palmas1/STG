#ifndef COMMON_H
#define COMMON_H

#include "precompiled.h"

#ifdef OS_WIN
	#ifdef STG_SHARED_LIB_EXPORTS
		#define STG_SHARED_LIB_API __declspec(dllexport)
	#else
		#define STG_SHARED_LIB_API __declspec(dllimport)
	#endif
#else
    #define STG_SHARED_LIB_API 
#endif // OS_WIN

#define GET_INDEX(i, j, k, is, js, ks) ks * js * i + ks * j + k  

typedef float STG_float;
typedef unsigned long int STG_int;


typedef struct Mesh_s
{
	STG_float * x;
	STG_float * y;
	STG_float * z;

} Mesh;

typedef struct ReStress_s
{
	STG_float * re_uu;
	STG_float * re_vv;
	STG_float * re_ww;
	STG_float * re_uv;
    STG_float * re_uw;
	STG_float * re_vw;

} ReStress;

typedef struct Scales_s
{
	STG_float * length_scale;
	STG_float * time_scale;

} Scales;

typedef struct InitData_s
{
	STG_int i_cnt;
	STG_int j_cnt;
	STG_int k_cnt;
	Mesh mesh;
    ReStress re;
    Scales scales;

} InitData;

typedef struct OutData_s
{
	STG_float * time;
	STG_int num_ts;
	STG_int i_cnt;
	STG_int j_cnt;
	STG_int k_cnt;
	STG_float ** u_p;
	STG_float ** v_p;
	STG_float ** w_p;

} OutData;

void free_InitData(InitData * init_data);

void free_OutData(OutData * out_data);

STG_SHARED_LIB_API void test_func(STG_float * arr, STG_int num);

STG_SHARED_LIB_API void init_rand();

STG_SHARED_LIB_API STG_float * get_uniform(STG_float min, STG_float max, STG_int num);

STG_SHARED_LIB_API STG_float bisection(STG_float x1, STG_float x2, STG_float(*func)(STG_float), STG_float a, STG_float eps, STG_int max_iter);

STG_float box_muller();

STG_SHARED_LIB_API STG_float * get_normal(STG_float mu, STG_float sigma, STG_int num);

//STG_float inverse_transform(STG_float u, STG_float(*cdf)(STG_float));

STG_SHARED_LIB_API STG_float * get_trigon(STG_int num);

void eig_3x3_sym(STG_float m11, STG_float m22, STG_float m33, STG_float m12, STG_float m13, STG_float m23, STG_float * eig_vals);

void compute_eig_vector1(
	STG_float m11, STG_float m22, STG_float m33, STG_float m12, STG_float m13, STG_float m23,
	STG_float eig_val1, STG_float * eig_vec1
);

void compute_ort_component(STG_float W[3], STG_float * U, STG_float * V);

void compute_eig_vector2(
	STG_float m11, STG_float m22, STG_float m33, STG_float m12, STG_float m13, STG_float m23,
	STG_float eig_vec1[3], STG_float eig_val2, STG_float U[3], STG_float V[3], STG_float * eig_vec2
);

void cross(STG_float v1[3], STG_float v2[3], STG_float *cross);

STG_float dot(STG_float v1[3], STG_float v2[3]);

void compute_eig(
	STG_float m11, STG_float m22, STG_float m33, STG_float m12, STG_float m13, STG_float m23,
	STG_float * eig_vals, STG_float * eig_vec1, STG_float * eig_vec2, STG_float * eig_vec3
);

#endif // !COMMON_H
