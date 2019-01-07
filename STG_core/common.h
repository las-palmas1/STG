#ifndef COMMON_H
#define COMMON_H

#define GET_INDEX(i, j, k, is, js, ks) ks*js*i + ks*j + k  

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
	STG_int i_cnt;
	STG_int j_cnt;
	STG_int k_cnt;
	STG_float * u_p;
	STG_float * v_p;
	STG_float * w_p;

} OutData;


STG_float * get_uniform(STG_float min, STG_float max, STG_int num);

STG_float bisection(STG_float x1, STG_float x2, STG_float (*mono_func)(STG_float), STG_float eps, STG_int max_iter);

//STG_float * get_normal(STG_float mu, STG_float sigma, STG_int num);





#endif // !COMMON_H
