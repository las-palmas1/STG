#include "precompiled.h"
#include "common.h"
#include <stdlib.h>
#include <stdio.h>


void STG_free_InitData(STG_InitData * init_data)
{
	free(init_data->mesh.x);
	free(init_data->mesh.y);
	free(init_data->mesh.z);
	free(init_data->re.re_uu);
	free(init_data->re.re_vv);
	free(init_data->re.re_ww);
	free(init_data->re.re_uv);
	free(init_data->re.re_uw);
	free(init_data->re.re_vw);
	free(init_data->scales.length_scale);
	free(init_data->scales.time_scale);
}

void STG_free_OutData(STG_OutData * out_data)
{
	free(out_data->time);
	free(out_data->u_p);
	free(out_data->v_p);
	free(out_data->w_p);
}

void STG_free_OutDataNode(STG_OutDataNode * out_data)
{
	free(out_data->time);
	free(out_data->u_p);
	free(out_data->v_p);
	free(out_data->w_p);
}

void STG_free_OutDataTS(STG_OutDataTS * out_data)
{
	free(out_data->u_p);
	free(out_data->v_p);
	free(out_data->w_p);
}

void STG_init_rand()
{
	srand(time(NULL));
}

void test_func(STG_float * arr, STG_int num)
{
	for (STG_int i = 0; i < num; i++)
	{
		printf("%.2f \n", arr[i]);
	}
}

// Функция возвращает массив заданного числа точек равномерно распределенных на заданном интервале
STG_float * get_uniform(STG_float min, STG_float max, STG_int num)
{
    STG_float * res = (STG_float*)malloc(sizeof (STG_float) * num);
	for (int i = 0; i < num; i++)
    {
        res[i] = min + (STG_float)rand() / (STG_float)RAND_MAX * (max - min);
    }
    return res;
}

void get_uniform_ref(STG_float min, STG_float max, STG_int num, STG_float * res)
{
	for (int i = 0; i < num; i++)
	{
		res[i] = min + (STG_float)rand() / (STG_float)RAND_MAX * (max - min);
	}
}

static STG_int sign(STG_float val)
{
    if (val > 0.){
        return 1;
    }
    if (val == 0.){
        return 0;
    }
    if (val < 0){
        return -1;
    }
}

// Функция для нахождения с заданной точностью корня уравнения func(x) - a = 0 на заданном промежутке по алгоритму бисекции.
// func(x) - монотонная на промежутке (x1, x2) функция.
// Задается также максимальное число итераций max_iter и точность eps.
STG_float bisection(STG_float x1, STG_float x2, STG_float (*func)(STG_float), STG_float a, STG_float eps, STG_int max_iter)
{
    STG_float res;
    if (func(x1) - a == 0.){
        res = func(x1);
    }
    if (func(x2) - a == 0.){
        res = func(x2);
    }
    else {

        STG_float dx;
        STG_float x_start = x1;
        STG_float x_end = x2;
        STG_float x_inter = (x2 + x1) / 2;
        STG_int count = 0;
        STG_float r = x_end - x_start;

        while (r > eps) {
            dx = (x_end - x_start) / 2;
            x_inter = x_start + dx;
            if (sign(func(x_inter) - a) != sign(func(x_start) - a)){
                x_end = x_inter;
            }
            else {
                x_start = x_inter;
            }
            count += 1;
            r = x_end - x_start;
            if (count > max_iter){
                break;
            }
        }
        res = x_inter;
    }
    return res;
}

// Реализация преобразования Бокса-Мюллера для моделирования стандартного нормального распределения
STG_float box_muller()
{
	STG_float pi = 2 * asin(1);
	STG_float phi = (STG_float)rand() / (STG_float)RAND_MAX;
	STG_float r = (STG_float)rand() / (STG_float)RAND_MAX;
	STG_float res = cos(2 * pi * phi) * sqrt(-2 * log(r));
	return res;
}

STG_float * get_normal(STG_float mu, STG_float sigma, STG_int num)
{
	STG_float * res = (STG_float*)malloc(sizeof(STG_float) * num);
	for (STG_int i = 0; i < num; i++)
	{
		res[i] = mu + sigma * box_muller();
	}
	return res;
}

void get_normal_ref(STG_float mu, STG_float sigma, STG_int num, STG_float * res)
{
	for (STG_int i = 0; i < num; i++)
	{
		res[i] = mu + sigma * box_muller();
	}
}

// Рассчитывает массив случайных чисел с распределением вероятностей p(x) = 1/2 * sin(x) по методу обратного преобразования 
STG_float * get_trigon(STG_int num)
{
	STG_float * res = (STG_float*)malloc(sizeof(STG_float) * num);
	for (STG_int i = 0; i < num; i++)
	{
		STG_float u = (STG_float)rand() / (STG_float)RAND_MAX;
		res[i] = acos(1 - 2 * u);
	}
	return res;
}

void get_trigon_ref(STG_int num, STG_float * res)
{
	for (STG_int i = 0; i < num; i++)
	{
		STG_float u = (STG_float)rand() / (STG_float)RAND_MAX;
		res[i] = acos(1 - 2 * u);
	}
}

// Расчет собственных значений симметричной матрицы размером 3 x 3.
void eig_3x3_sym(STG_float m11, STG_float m22, STG_float m33, STG_float m12, STG_float m13, STG_float m23, STG_float * eig_vals)
{
	if (m11 == 0.)
	{
		m11 = 1e-7;
	}
	if (m22 == 0.)
	{
		m22 = 1e-7;
	}
	if (m33 == 0.)
	{
		m33 = 1e-7;
	}
	STG_float p1 = m12 * m12 + m13 * m13 + m23 * m23;
	if (p1 == 0)
	{
		eig_vals[0] = m11;
		eig_vals[1] = m22;
		eig_vals[2] = m33;
	}
	else
	{
		STG_float q = (m11 + m22 + m33) / (STG_float)3;
		STG_float p2 = (m11 - q)*(m11 - q) + (m22 - q)*(m22 - q) + (m33 - q)*(m33 - q) + 2 * p1;
		STG_float p = sqrt(p2 / (STG_float)6);
		STG_float B[3][3] = {
			{(1 / p) * (m11 - q), (1 / p) * (m12 - 0), (1 / p) * (m13 - 0) },
			{ (1 / p) * (m12 - 0), (1 / p) * (m22 - q), (1 / p) * (m23 - 0) },
			{ (1 / p) * (m13 - 0), (1 / p) * (m23 - 0), (1 / p) * (m33 - q) },
		};
		STG_float r = 0.5 * (
			B[0][0] * B[1][1] * B[2][2] - B[0][0] * B[1][2] * B[2][1] - B[0][1] * B[1][0] * B[2][2] + 
			B[0][1] * B[1][2] * B[2][0] + B[0][2] * B[1][0] * B[2][1] - B[0][2] * B[1][1] * B[2][0]
		);
		
		STG_float pi = 2 * asin(1);
		STG_float phi;
		if (r <= -1) {
			phi = pi / (STG_float)3;
		}
		if (r >= 1) {
			phi = (STG_float)0;
		}
		if (r > -1 && r < 1) {
			phi = (STG_float)acos(r) / (STG_float)3;
		}
		
		eig_vals[0] = q + 2 * p * cos(phi);
		eig_vals[2] = q + 2 * p * cos(phi + (STG_float)2 * pi / (STG_float)3);
		eig_vals[1] = 3 * q - eig_vals[0] - eig_vals[2];
	} 
}

// Векторное произведение векторов в 3-х мерном пространстве
void cross(STG_float v1[3], STG_float v2[3], STG_float *cross)
{
	cross[0] = v1[1] * v2[2] - v2[1] * v1[2];
	cross[1] = v2[0] * v1[2] - v1[0] * v2[2];
	cross[2] = v1[0] * v2[1] - v2[0] * v1[1];
}

// Скалярное произведение 3-х мерных векторов.
STG_float dot(STG_float v1[3], STG_float v2[3])
{
	return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

// Расчет собственного вектора для наибольшего собственного значения
void compute_eig_vector1(
	STG_float m11, STG_float m22, STG_float m33, STG_float m12, STG_float m13, STG_float m23,
	STG_float eig_val1, STG_float * eig_vec1
)
{
	STG_float row1[3] = { m11 - eig_val1, m12, m13 };
	STG_float row2[3] = { m12, m22 - eig_val1, m23 };
	STG_float row3[3] = { m13, m23, m33 - eig_val1 };
	STG_float r1xr2[3];
	STG_float r1xr3[3];
	STG_float r2xr3[3];
	cross(row1, row2, r1xr2);
	cross(row1, row3, r1xr3);
	cross(row2, row3, r2xr3);
	STG_float d1 = dot(r1xr2, r1xr2);
	STG_float d2 = dot(r1xr3, r1xr3);
	STG_float d3 = dot(r2xr3, r2xr3);
	STG_float dmax = d1;
	STG_int imax = 1;
	if (d2 > dmax) { dmax = d2; imax = 2; }
	if (d3 > dmax) { imax = 3; }
	if (imax == 1)
	{
		eig_vec1[0] = r1xr2[0] / sqrt(d1);
		eig_vec1[1] = r1xr2[1] / sqrt(d1);
		eig_vec1[2] = r1xr2[2] / sqrt(d1);
	}
	else if (imax == 2)
	{
		eig_vec1[0] = r1xr3[0] / sqrt(d2);
		eig_vec1[1] = r1xr3[1] / sqrt(d2);
		eig_vec1[2] = r1xr3[2] / sqrt(d2);
	}
	else
	{
		eig_vec1[0] = r2xr3[0] / sqrt(d3);
		eig_vec1[1] = r2xr3[1] / sqrt(d3);
		eig_vec1[2] = r2xr3[2] / sqrt(d3);
	}
}

// Расчет пары единичных взаимноперпендикулярных векторов, перпендикулярных первому собственному вектору W.
void compute_ort_component(STG_float W[3], STG_float * U, STG_float * V)
{
	STG_float inv_length;
	if (fabs(W[0]) > fabs(W[1]))
	{
		inv_length = 1 / sqrt(W[0] * W[0] + W[2] * W[2]);
		U[0] = -W[2] * inv_length;
		U[1] = 0;
		U[2] = W[0] * inv_length;
	}
	else
	{
		inv_length = 1 / sqrt(W[1] * W[1] + W[2] * W[2]);
		U[0] = 0;
		U[1] = W[2] * inv_length;
		U[2] = -W[1] * inv_length;
	}
	cross(W, U, V);
}

// Расчет собственного вектора, соотвествующего второму по величине собственному значению.
void compute_eig_vector2(
	STG_float m11, STG_float m22, STG_float m33, STG_float m12, STG_float m13, STG_float m23, 
	STG_float eig_vec1[3], STG_float eig_val2, STG_float U[3], STG_float V[3], STG_float * eig_vec2
)
{
	STG_float MU[3] = {
		m11 * U[0] + m12 * U[1] + m13 * U[2],
		m12 * U[0] + m22 * U[1] + m23 * U[2],
		m13 * U[0] + m23 * U[1] + m33 * U[2],
	};
	STG_float MV[3] = {
		m11 * V[0] + m12 * V[1] + m13 * V[2],
		m12 * V[0] + m22 * V[1] + m23 * V[2],
		m13 * V[0] + m23 * V[1] + m33 * V[2],
	};
	STG_float k11 = dot(U, MU) - eig_val2;
	STG_float k12 = dot(U, MV);
	STG_float k22 = dot(V, MV) - eig_val2;

	STG_float absK11 = fabs(k11), absK12 = fabs(k12), absK22 = fabs(k22);
	if (absK11 >= absK22)
	{
		STG_float maxAbsComp = max(absK11, absK12);
		if (maxAbsComp > 0)
		{
			if (absK11 >= absK12)
			{
				k12 /= k11;   k11 = (STG_float)1 / sqrt(1 + k12 * k12);   k12 *= k11;
			}
			else
			{
				k11 /= k12;   k12 = (STG_float)1 / sqrt(1 + k11 * k11);   k11 *= k12;
			}
			eig_vec2[0] = k12 * U[0] - k11 * V[0];
			eig_vec2[1] = k12 * U[1] - k11 * V[1];
			eig_vec2[2] = k12 * U[2] - k11 * V[2];
		}
		else
		{
			eig_vec2[0] = U[0];
			eig_vec2[1] = U[1];
			eig_vec2[2] = U[2];
		}
	}
	else
	{
		STG_float maxAbsComp = max(absK22, absK12);
		if (maxAbsComp > 0)
		{
			if (absK22 >= absK12)
			{
				k12 /= k22;   k22 = (STG_float)1 / sqrt(1 + k12 * k12);   k12 *= k22;
			}
			else
			{
				k22 /= k12;   k12 = (STG_float)1 / sqrt(1 + k22 * k22);   k22 *= k12;
			}
			eig_vec2[0] = k22 * U[0] - k12 * V[0];
			eig_vec2[1] = k22 * U[1] - k12 * V[1];
			eig_vec2[2] = k22 * U[2] - k12 * V[2];
		}
		else
		{
			eig_vec2[0] = U[0];
			eig_vec2[1] = U[1];
			eig_vec2[2] = U[2];
		}
	}
}

void compute_eig(
	STG_float m11, STG_float m22, STG_float m33, STG_float m12, STG_float m13, STG_float m23,
	STG_float * eig_vals, STG_float * eig_vec1, STG_float * eig_vec2, STG_float * eig_vec3
)
{
	if (m12 == 0 && m13 == 0 && m23 == 0)
	{
		eig_vals[0] = m11;
		eig_vals[1] = m22;
		eig_vals[2] = m33;
		eig_vec1[0] = 1; eig_vec1[1] = 0; eig_vec1[2] = 0;
		eig_vec2[0] = 0; eig_vec2[1] = 1; eig_vec2[2] = 0;
		eig_vec3[0] = 0; eig_vec3[1] = 0; eig_vec3[2] = 1;
	}
	else
	{
		STG_float U[3], V[3];
		eig_3x3_sym(m11, m22, m33, m12, m13, m23, eig_vals);
		compute_eig_vector1(m11, m22, m33, m12, m13, m23, eig_vals[0], eig_vec1);
		compute_ort_component(eig_vec1, U, V);
		compute_eig_vector2(m11, m22, m33, m12, m13, m23, eig_vec1, eig_vals[1], U, V, eig_vec2);
		cross(eig_vec1, eig_vec2, eig_vec3);
	}
}


