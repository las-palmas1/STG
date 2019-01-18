#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include "Smirnov.h"


static STG_float test_func_linear(STG_float x)
{
    return x - 1;
}


// Экстремум в точке x = -1
static STG_float test_func_pow2(STG_float x)
{
    return (x - 2) * (x + 4);
}

// Экстремум в x=-3 и x=5
static STG_float test_func_pow3(STG_float x)
{
    return (x - 1) * (x - 8) * (x + 6);
}

static void bisection_wrap(
        STG_float x1, STG_float x2, STG_float (*func)(STG_float), STG_float a,
        STG_float eps, STG_int max_iter, char * func_name, STG_float expect_res
    )
{
    printf("Test bisection. Function:     %s\n", func_name);
    printf("Interval:   x1 = %.3f,   x2 = %.3f\n", x1, x2);
    printf("eps = %.4f,  max_iter = %d\n", eps, max_iter);

    STG_float res = bisection(x1, x2, func, a, eps, max_iter);

    printf("Result:   %.4f\n", res);
    printf("Expected:   %.4f\n", expect_res);
    printf("\n");
}

static void test_eig_values_and_vectors_computing(STG_float m11, STG_float m22, STG_float m33, STG_float m12, STG_float m13, STG_float m23)
{
	printf("Eigvalues test.\n");
	printf("Matrix: \n");
	printf("%.3f  %.3f  %.3f \n", m11, m12, m13);
	printf("%.3f  %.3f  %.3f \n", m12, m22, m23);
	printf("%.3f  %.3f  %.3f \n", m13, m23, m33);
	STG_float eig_vals[3];
	STG_float eig_vec1[3], eig_vec2[3], eig_vec3[3];
	compute_eig(m11, m22, m33, m12, m13, m23, eig_vals, eig_vec1, eig_vec2, eig_vec3);

	STG_float M1[3][3] = {
		{ m11 - eig_vals[0], m12, m13 },
		{ m12, m22 - eig_vals[1], m23 },
		{ m13, m23, m33 - eig_vals[2] },
	};
	STG_float det = (
		M1[0][0] * M1[1][1] * M1[2][2] - M1[0][0] * M1[1][2] * M1[2][1] - M1[0][1] * M1[1][0] * M1[2][2] +
		M1[0][1] * M1[1][2] * M1[2][0] + M1[0][2] * M1[1][0] * M1[2][1] - M1[0][2] * M1[1][1] * M1[2][0]
	);
	STG_float M_x_eig_vec1[3] = {
		m11 * eig_vec1[0] + m12 * eig_vec1[1] + m13 * eig_vec1[2],
		m12 * eig_vec1[0] + m22 * eig_vec1[1] + m23 * eig_vec1[2],
		m13 * eig_vec1[0] + m23 * eig_vec1[1] + m33 * eig_vec1[2],
	};
	STG_float eig_val1_x_eig_vec1[3] = {
		eig_vals[0] * eig_vec1[0], eig_vals[0] * eig_vec1[1], eig_vals[0] * eig_vec1[2]
	};
	STG_float M_x_eig_vec2[3] = {
		m11 * eig_vec2[0] + m12 * eig_vec2[1] + m13 * eig_vec2[2],
		m12 * eig_vec2[0] + m22 * eig_vec2[1] + m23 * eig_vec2[2],
		m13 * eig_vec2[0] + m23 * eig_vec2[1] + m33 * eig_vec2[2],
	};
	STG_float eig_val2_x_eig_vec2[3] = {
		eig_vals[1] * eig_vec2[0], eig_vals[1] * eig_vec2[1], eig_vals[1] * eig_vec2[2]
	};
	STG_float M_x_eig_vec3[3] = {
		m11 * eig_vec3[0] + m12 * eig_vec3[1] + m13 * eig_vec3[2],
		m12 * eig_vec3[0] + m22 * eig_vec3[1] + m23 * eig_vec3[2],
		m13 * eig_vec3[0] + m23 * eig_vec3[1] + m33 * eig_vec3[2],
	};
	STG_float eig_val3_x_eig_vec3[3] = {
		eig_vals[2] * eig_vec3[0], eig_vals[2] * eig_vec3[1], eig_vals[2] * eig_vec3[2]
	};

	printf("\n");
	printf("eig1 = %.3f\n", eig_vals[0]);
	printf("eig2 = %.3f\n", eig_vals[1]);
	printf("eig3 = %.3f\n", eig_vals[2]);
	printf("det(M - eig*I) = %.4f\n", det);
	printf("\n");
	printf("Eigvectors:\n");
	printf("%.3f  %.3f  %.3f\n", eig_vec1[0], eig_vec2[0], eig_vec3[0]);
	printf("%.3f  %.3f  %.3f\n", eig_vec1[1], eig_vec2[1], eig_vec3[1]);
	printf("%.3f  %.3f  %.3f\n", eig_vec1[2], eig_vec2[2], eig_vec3[2]);
	printf("(M - eig1*I)*V1 = \n");
	printf("%.3f \n", M_x_eig_vec1[0] - eig_val1_x_eig_vec1[0]);
	printf("%.3f \n", M_x_eig_vec1[1] - eig_val1_x_eig_vec1[1]);
	printf("%.3f \n", M_x_eig_vec1[2] - eig_val1_x_eig_vec1[2]);
	printf("(M - eig2*I)*V2 = \n");
	printf("%.3f \n", M_x_eig_vec2[0] - eig_val2_x_eig_vec2[0]);
	printf("%.3f \n", M_x_eig_vec2[1] - eig_val2_x_eig_vec2[1]);
	printf("%.3f \n", M_x_eig_vec2[2] - eig_val2_x_eig_vec2[2]);
	printf("(M - eig3*I)*V3 = \n");
	printf("%.3f \n", M_x_eig_vec3[0] - eig_val3_x_eig_vec3[0]);
	printf("%.3f \n", M_x_eig_vec3[1] - eig_val3_x_eig_vec3[1]);
	printf("%.3f \n", M_x_eig_vec3[2] - eig_val3_x_eig_vec3[2]);
	printf("\n");

}

static void test_bisection()
{
    STG_float eps = 0.0001;
    STG_float max_iter = 100;

    bisection_wrap(0, 2, test_func_linear, 0, eps, max_iter, "f(x) = x - 1", 1);

    bisection_wrap(-8, -1, test_func_pow2, 0, eps, max_iter, "f(x) = (x - 2) * (x + 4)", -4);
    bisection_wrap(-1, 5, test_func_pow2, 0, eps, max_iter, "f(x) = (x - 2) * (x + 4)", -2);

    bisection_wrap(-9, -3, test_func_pow3, 0, eps, max_iter, "f(x) = (x - 1) * (x - 8) * (x + 6)", -6);
    bisection_wrap(-3, 5, test_func_pow3, 0, eps, max_iter, "f(x) = (x - 1) * (x - 8) * (x + 6)", 1);
    bisection_wrap(5, 12, test_func_pow3, 0, eps, max_iter, "f(x) = (x - 1) * (x - 8) * (x + 6)", 8);
}

static void test_uniform(STG_int num, STG_float min, STG_float max)
{
	printf("Test uniform distribution. Samples number = %d \n", num);
	printf("Interval:   min = %.2f    max = %.2f \n", min, max);
	STG_float * samp = get_uniform(min, max, num);
	for (STG_int i = 0; i < num; i++)
	{
		printf("%.4f \n", samp[i]);
	}
	free(samp);
	printf("\n");
}

static void test_normal(STG_int num, STG_float mu, STG_float sigma)
{
	printf("Test normal distribution. Samples number = %d \n", num);
	printf("mu = %.2f    sigma = %.2f \n", mu, sigma);
	STG_float * samp = get_normal(mu, sigma, num);
	for (STG_int i = 0; i < num; i++)
	{
		printf("%.4f \n", samp[i]);
	}
	free(samp);
	printf("\n");
}

static void test_trigon(STG_int num)
{
	printf("Test trigon distribution. Samples number = %d \n", num);
	STG_float * samp = get_trigon(num);
	for (STG_int i = 0; i < num; i++)
	{
		printf("%.4f \n", samp[i]);
	}
	free(samp);
	printf("\n");
}


static void test_Smirnov_pulsation(
	STG_float re_uu, STG_float re_vv, STG_float re_ww, STG_float re_uv, STG_float re_uw, STG_float re_vw, STG_int num_modes
)
{
	printf("Test Smirnov pulsations.\n");
	printf("Reynolds stresses:\n");
	printf("%.2f  %.2f  %.2f \n", re_uu, re_uv, re_uw);
	printf("%.2f  %.2f  %.2f \n", re_uv, re_vv, re_vw);
	printf("%.2f  %.2f  %.2f \n", re_uw, re_vw, re_ww);
	printf("Modes number: %d \n", num_modes);
	InitData init_data;
	init_data.i_cnt = 1;
	init_data.j_cnt = 1;
	init_data.k_cnt = 1;
	STG_int num = init_data.i_cnt * init_data.j_cnt * init_data.k_cnt;

	init_data.mesh.x = (STG_float*)malloc(sizeof(STG_float) * num);
	init_data.mesh.y = (STG_float*)malloc(sizeof(STG_float) * num);
	init_data.mesh.z = (STG_float*)malloc(sizeof(STG_float) * num);
	init_data.mesh.x[0] = 1;
	init_data.mesh.y[0] = 1;
	init_data.mesh.z[0] = 1;

	init_data.re.re_uu = (STG_float*)malloc(sizeof(STG_float) * num);
	init_data.re.re_vv = (STG_float*)malloc(sizeof(STG_float) * num);
	init_data.re.re_ww = (STG_float*)malloc(sizeof(STG_float) * num);
	init_data.re.re_uv = (STG_float*)malloc(sizeof(STG_float) * num);
	init_data.re.re_uw = (STG_float*)malloc(sizeof(STG_float) * num);
	init_data.re.re_vw = (STG_float*)malloc(sizeof(STG_float) * num);
	init_data.re.re_uu[0] = re_uu;
	init_data.re.re_vv[0] = re_vv;
	init_data.re.re_ww[0] = re_ww;
	init_data.re.re_uv[0] = re_uv;
	init_data.re.re_uw[0] = re_uw;
	init_data.re.re_vw[0] = re_vw;

	init_data.scales.length_scale = (STG_float*)malloc(sizeof(STG_float) * num);
	init_data.scales.time_scale = (STG_float*)malloc(sizeof(STG_float) * num);
	init_data.scales.length_scale[0] = 1;
	init_data.scales.time_scale[0] = 1;

	SmirnovData data;
	OutData out_data;

	compute_Smirnov_data(init_data, num_modes, 0.1, 0, &data);

	compute_Smirnov_field(init_data, data, &out_data);

	printf("u = %.3f   v = %.3f  w = %.3f \n", out_data.u_p[0][0], out_data.v_p[0][0], out_data.w_p[0][0]);
	printf("\n");

	free_Smirnov_data(&data);
	free_InitData(&init_data);
	free_OutData(&out_data);
}


static void test_Smirnov_pulsation_node(
	STG_float re_uu, STG_float re_vv, STG_float re_ww, STG_float re_uv, STG_float re_uw, STG_float re_vw, STG_int num_modes
)
{
	printf("Test Smirnov pulsations at node.\n");
	printf("Reynolds stresses:\n");
	printf("%.2f  %.2f  %.2f \n", re_uu, re_uv, re_uw);
	printf("%.2f  %.2f  %.2f \n", re_uv, re_vv, re_vw);
	printf("%.2f  %.2f  %.2f \n", re_uw, re_vw, re_ww);
	printf("Modes number: %d \n", num_modes);
	InitData init_data;
	init_data.i_cnt = 1;
	init_data.j_cnt = 1;
	init_data.k_cnt = 1;
	STG_int num = init_data.i_cnt * init_data.j_cnt * init_data.k_cnt;

	init_data.mesh.x = (STG_float*)malloc(sizeof(STG_float) * num);
	init_data.mesh.y = (STG_float*)malloc(sizeof(STG_float) * num);
	init_data.mesh.z = (STG_float*)malloc(sizeof(STG_float) * num);
	init_data.mesh.x[0] = 1;
	init_data.mesh.y[0] = 1;
	init_data.mesh.z[0] = 1;

	init_data.re.re_uu = (STG_float*)malloc(sizeof(STG_float) * num);
	init_data.re.re_vv = (STG_float*)malloc(sizeof(STG_float) * num);
	init_data.re.re_ww = (STG_float*)malloc(sizeof(STG_float) * num);
	init_data.re.re_uv = (STG_float*)malloc(sizeof(STG_float) * num);
	init_data.re.re_uw = (STG_float*)malloc(sizeof(STG_float) * num);
	init_data.re.re_vw = (STG_float*)malloc(sizeof(STG_float) * num);
	init_data.re.re_uu[0] = re_uu;
	init_data.re.re_vv[0] = re_vv;
	init_data.re.re_ww[0] = re_ww;
	init_data.re.re_uv[0] = re_uv;
	init_data.re.re_uw[0] = re_uw;
	init_data.re.re_vw[0] = re_vw;

	init_data.scales.length_scale = (STG_float*)malloc(sizeof(STG_float) * num);
	init_data.scales.time_scale = (STG_float*)malloc(sizeof(STG_float) * num);
	init_data.scales.length_scale[0] = 1;
	init_data.scales.time_scale[0] = 1;

	SmirnovData data;
	OutDataNode out_data;

	compute_Smirnov_data(init_data, num_modes, 0.1, 0, &data);

	compute_Smirnov_field_node(init_data, data, &out_data, 0, 0, 0);

	printf("u = %.3f   v = %.3f  w = %.3f \n", out_data.u_p[0], out_data.v_p[0], out_data.w_p[0]);
	printf("\n");

	free_Smirnov_data(&data);
	free_InitData(&init_data);
	free_OutData(&out_data);
}

int main(int argc, char * argv[])
{
	init_rand();
    test_bisection();
	test_uniform(20, -5, 5);
	test_normal(20, 0, 1);
	test_trigon(20);

	test_eig_values_and_vectors_computing(1, 1, 1, 1, 1, 1);
	test_eig_values_and_vectors_computing(1, 1, 1, 2, 1, 1);
	test_eig_values_and_vectors_computing(1, 1, 3, 3, 1, 4);
	test_eig_values_and_vectors_computing(2, 3, 9, 0, 0, 4);
	test_eig_values_and_vectors_computing(5, -3, 0, 0, -4, 4);
	test_eig_values_and_vectors_computing(1, 1, 3, 0, 0, 0);

	test_Smirnov_pulsation(1, 3, 2, 0, 0, 0, 150);
	test_Smirnov_pulsation(1, 1, 1, 0, 0, 0, 150);
	test_Smirnov_pulsation(1, 1, 1, 2, 1, -3, 150);
	test_Smirnov_pulsation_node(1, 1, 1, 2, 1, -3, 150);

	return 0;
}
