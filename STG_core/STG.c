#include <stdio.h>
#include <stdlib.h>
#ifdef __linux__
    #include "lib/Smirnov.h"
    #include "lib/common.h"
	#include "lib/Davidson.h"
#else
    #include "lib\common.h"
    #include "lib\Smirnov.h"
	#include "lib\Davidson.h"
#endif


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


static void compute_mesh(
	STG_float x_size, STG_float y_size, STG_float z_size, 
	STG_float i_cnt, STG_float j_cnt, STG_float k_cnt, STG_Mesh * mesh
)
{
	STG_int size = i_cnt * j_cnt * k_cnt;
	mesh->x = (STG_float *)malloc(sizeof(STG_float) * size);
	mesh->y = (STG_float *)malloc(sizeof(STG_float) * size);
	mesh->z = (STG_float *)malloc(sizeof(STG_float) * size);
	STG_float x_step = x_size / (i_cnt - 1);
	STG_float y_step = y_size / (j_cnt - 1);
	STG_float z_step = z_size / (k_cnt - 1);
	for (STG_int i = 0; i < i_cnt; i++)
	{
		for (STG_int j = 0; j < j_cnt; j++)
		{
			for (STG_int k = 0; k < k_cnt; k++)
			{
				STG_int num = GET_INDEX(i, j, k, i_cnt, j_cnt, k_cnt);
				mesh->x[num] = i * x_step;
				mesh->y[num] = j * y_step;
				mesh->z[num] = k * z_step;
			}
		}
	}
}

static void fill_re(
	STG_float i_cnt, STG_float j_cnt, STG_float k_cnt,
	STG_float re_uu, STG_float re_vv, STG_float re_ww, 
	STG_float re_uv, STG_float re_uw, STG_float re_vw,
	STG_ReStress * re
)
{
	STG_int size = i_cnt * j_cnt * k_cnt;
	re->re_uu = (STG_float *)malloc(sizeof(STG_float) * size);
	re->re_vv = (STG_float *)malloc(sizeof(STG_float) * size);
	re->re_ww = (STG_float *)malloc(sizeof(STG_float) * size);
	re->re_uv = (STG_float *)malloc(sizeof(STG_float) * size);
	re->re_uw = (STG_float *)malloc(sizeof(STG_float) * size);
	re->re_vw = (STG_float *)malloc(sizeof(STG_float) * size);

	for (STG_int i = 0; i < size; i++)
	{
		re->re_uu[i] = re_uu;
		re->re_vv[i] = re_vv;
		re->re_ww[i] = re_ww;
		re->re_uv[i] = re_uv;
		re->re_uw[i] = re_uw;
		re->re_vw[i] = re_vw;
	}
}

static void fill_scales(
	STG_float i_cnt, STG_float j_cnt, STG_float k_cnt, STG_float ls_i,
	STG_float ls_ux, STG_float ls_uy, STG_float ls_uz,
	STG_float ls_vx, STG_float ls_vy, STG_float ls_vz,
	STG_float ls_wx, STG_float ls_wy, STG_float ls_wz,
	STG_float ts_i, STG_float ts_u, STG_float ts_v, STG_float ts_w, STG_Scales * scales
)
{
	STG_int size = i_cnt * j_cnt * k_cnt;
	scales->ls_i = (STG_float *)malloc(sizeof(STG_float) * size);
	scales->ls_ux = (STG_float *)malloc(sizeof(STG_float) * size);
	scales->ls_uy = (STG_float *)malloc(sizeof(STG_float) * size);
	scales->ls_uz = (STG_float *)malloc(sizeof(STG_float) * size);
	scales->ls_vx = (STG_float *)malloc(sizeof(STG_float) * size);
	scales->ls_vy = (STG_float *)malloc(sizeof(STG_float) * size);
	scales->ls_vz = (STG_float *)malloc(sizeof(STG_float) * size);
	scales->ls_wx = (STG_float *)malloc(sizeof(STG_float) * size);
	scales->ls_wy = (STG_float *)malloc(sizeof(STG_float) * size);
	scales->ls_wz = (STG_float *)malloc(sizeof(STG_float) * size);
	scales->ts_i = (STG_float *)malloc(sizeof(STG_float) * size);
	scales->ts_u = (STG_float *)malloc(sizeof(STG_float) * size);
	scales->ts_v = (STG_float *)malloc(sizeof(STG_float) * size);
	scales->ts_w = (STG_float *)malloc(sizeof(STG_float) * size);
	for (STG_int i = 0; i < size; i++)
	{
		scales->ls_i[i] = ls_i;
		scales->ls_ux[i] = ls_ux;
		scales->ls_uy[i] = ls_uy;
		scales->ls_uz[i] = ls_uz;
		scales->ls_vx[i] = ls_vx;
		scales->ls_vy[i] = ls_vy;
		scales->ls_vz[i] = ls_vz;
		scales->ls_wx[i] = ls_wx;
		scales->ls_wy[i] = ls_wy;
		scales->ls_wz[i] = ls_wz;
		scales->ts_i[i] = ts_i;
		scales->ts_u[i] = ts_u;
		scales->ts_v[i] = ts_v;
		scales->ts_w[i] = ts_w;
	}
	
}

static void test_Davidson_mom_filed(
	STG_int node_cnt, STG_float length, STG_float re_uu, STG_float re_vv, STG_float re_ww,
	STG_float re_uv, STG_float re_uw, STG_float re_vw, 
	STG_float ls_i, STG_int num_modes
)
{
	printf("Test Davidson moment field.\n");
	printf("Reynolds stresses:\n");
	printf("%.2f  %.2f  %.2f \n", re_uu, re_uv, re_uw);
	printf("%.2f  %.2f  %.2f \n", re_uv, re_vv, re_vw);
	printf("%.2f  %.2f  %.2f \n", re_uw, re_vw, re_ww);
	printf("Modes number: %d \n", num_modes);
	STG_InitData init_data;
	init_data.i_cnt = node_cnt;
	init_data.j_cnt = node_cnt;
	init_data.k_cnt = node_cnt;
	STG_int num = init_data.i_cnt * init_data.j_cnt * init_data.k_cnt;
	STG_float ts_i = 0.01;
	compute_mesh(length, length, length, node_cnt, node_cnt, node_cnt, &(init_data.mesh));
	fill_re(node_cnt, node_cnt, node_cnt, re_uu, re_vv, re_ww, re_uv, re_uw, re_vw, &(init_data.re));
	fill_scales(node_cnt, node_cnt, node_cnt, ls_i, 0, 0, 0, 0, 0, 0, 0, 0, 0, ts_i, 0, 0, 0, &(init_data.scales));
	

	STG_float time = 0.;
	STG_float visc = 1.8e-5;
	STG_float dissip = pow(0.09, 0.75) * pow((re_uu + re_vv + re_ww) / 2, 1.5) / ls_i;
	STG_DavidsonData_Stationary stat_data;
	STG_DavidsonData_Transient trans_data;
	STG_VelMomField mom_field;

	STG_alloc_Davidson_trans_data(init_data, num_modes, &trans_data);
	STG_compute_Davidson_stat_data(init_data, num_modes, dissip, visc, 0.1, &stat_data);
	STG_compute_Davidson_trans_data(stat_data, num_modes, 0.1, 1, &trans_data);
	STG_compute_Davidson_moment_field(init_data, stat_data, &trans_data, time, &mom_field);

	printf("First 5 u_abs: \n");
	for (STG_int i = 0; i < 5; i++)
	{
		printf("i_mode = %d  u_abs = %f \n", i, stat_data.u_abs[i]);
	}
	printf("\n");
	printf("u = %.3f   v = %.3f  w = %.3f \n", mom_field.u_p[0], mom_field.v_p[0], mom_field.w_p[0]);
	printf("\n");

	STG_free_Davidson_stat_data(&stat_data);
	STG_free_Davidson_trans_data(&trans_data);
	STG_free_InitData(&init_data);
	STG_free_VelMomField(&mom_field);
}


static void test_Davidson_node_hist(
	STG_float ts, STG_int num_ts, STG_float ts_i, STG_float re_uu, STG_float re_vv, STG_float re_ww,
	STG_float re_uv, STG_float re_uw, STG_float re_vw, STG_int num_modes
)
{
	STG_int node_cnt = 10;
	STG_float length = 0.1;
	STG_float ls_i = 0.01;
	printf("Test Davidson node history.\n");
	printf("Reynolds stresses:\n");
	printf("%.2f  %.2f  %.2f \n", re_uu, re_uv, re_uw);
	printf("%.2f  %.2f  %.2f \n", re_uv, re_vv, re_vw);
	printf("%.2f  %.2f  %.2f \n", re_uw, re_vw, re_ww);
	printf("Modes number: %d \n", num_modes);
	STG_InitData init_data;
	init_data.i_cnt = node_cnt;
	init_data.j_cnt = node_cnt;
	init_data.k_cnt = node_cnt;
	STG_int num = init_data.i_cnt * init_data.j_cnt * init_data.k_cnt;
	compute_mesh(length, length, length, node_cnt, node_cnt, node_cnt, &(init_data.mesh));
	fill_re(node_cnt, node_cnt, node_cnt, re_uu, re_vv, re_ww, re_uv, re_uw, re_vw, &(init_data.re));
	fill_scales(node_cnt, node_cnt, node_cnt, ls_i, 0, 0, 0, 0, 0, 0, 0, 0, 0, ts_i, 0, 0, 0, &(init_data.scales));


	STG_float visc = 1.8e-5;
	STG_float dissip = pow(0.09, 0.75) * pow((re_uu + re_vv + re_ww) / 2, 1.5) / ls_i;
	STG_DavidsonData_Stationary stat_data;
	STG_DavidsonData_Transient trans_data;
	STG_VelNodeHist node_hist;

	STG_alloc_Davidson_trans_data(init_data, num_modes, &trans_data);
	STG_compute_Davidson_stat_data(init_data, num_modes, dissip, visc, ts, &stat_data);
	STG_compute_Davidson_node_hist(init_data, stat_data, ts, num_ts, &trans_data, &node_hist, 1, 1, 1);

	printf("Velosity history\n");
	for (STG_int i = 0; i < num_ts + 1; i++)
	{
		printf("time = %.3f  u = %.3f  v = %.3f  w = %.3f \n", node_hist.time[i], node_hist.u_p[i], node_hist.v_p[i], node_hist.w_p[i]);
	}
	printf("\n");

	STG_free_Davidson_stat_data(&stat_data);
	STG_free_Davidson_trans_data(&trans_data);
	STG_free_InitData(&init_data);
	STG_free_VelNodeHist(&node_hist);
}

static void test_Smirnov_mom_field(
	STG_int node_cnt, STG_float length, STG_float re_uu, STG_float re_vv, STG_float re_ww,
	STG_float re_uv, STG_float re_uw, STG_float re_vw,
	STG_float ls_i, STG_int num_modes
)
{
	printf("Test Smirnov pulsations.\n");
	printf("Reynolds stresses:\n");
	printf("%.2f  %.2f  %.2f \n", re_uu, re_uv, re_uw);
	printf("%.2f  %.2f  %.2f \n", re_uv, re_vv, re_vw);
	printf("%.2f  %.2f  %.2f \n", re_uw, re_vw, re_ww);
	printf("Modes number: %d \n", num_modes);
	STG_InitData init_data;
	init_data.i_cnt = node_cnt;
	init_data.j_cnt = node_cnt;
	init_data.k_cnt = node_cnt;
	STG_int num = init_data.i_cnt * init_data.j_cnt * init_data.k_cnt;
	STG_float ts_i = 0.01;
	compute_mesh(length, length, length, node_cnt, node_cnt, node_cnt, &(init_data.mesh));
	fill_re(node_cnt, node_cnt, node_cnt, re_uu, re_vv, re_ww, re_uv, re_uw, re_vw, &(init_data.re));
	fill_scales(node_cnt, node_cnt, node_cnt, ls_i, 0, 0, 0, 0, 0, 0, 0, 0, 0, ts_i, 0, 0, 0, &(init_data.scales));


    STG_SmirnovData data;
    STG_VelMomField mom_field;
    STG_float time = 0.0;

    STG_compute_Smirnov_data(init_data, num_modes, &data);

    STG_compute_Smirnov_moment_field(init_data, data, time, &mom_field);

    printf("u = %.3f   v = %.3f  w = %.3f \n", mom_field.u_p[0], mom_field.v_p[0], mom_field.w_p[0]);
	printf("\n");

    STG_free_Smirnov_data(&data);
	STG_free_InitData(&init_data);
    STG_free_VelMomField(&mom_field);
}


static void test_Smirnov_node_hist(
	STG_float re_uu, STG_float re_vv, STG_float re_ww, STG_float re_uv, STG_float re_uw, STG_float re_vw, STG_int num_modes
)
{
	printf("Test Smirnov pulsations at node.\n");
	printf("Reynolds stresses:\n");
	printf("%.2f  %.2f  %.2f \n", re_uu, re_uv, re_uw);
	printf("%.2f  %.2f  %.2f \n", re_uv, re_vv, re_vw);
	printf("%.2f  %.2f  %.2f \n", re_uw, re_vw, re_ww);
	printf("Modes number: %d \n", num_modes);
	STG_InitData init_data;
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

    init_data.scales.ls_i = (STG_float*)malloc(sizeof(STG_float) * num);
    init_data.scales.ls_ux = (STG_float*)malloc(sizeof(STG_float) * num);
    init_data.scales.ls_uy = (STG_float*)malloc(sizeof(STG_float) * num);
    init_data.scales.ls_uz = (STG_float*)malloc(sizeof(STG_float) * num);
    init_data.scales.ls_vx = (STG_float*)malloc(sizeof(STG_float) * num);
    init_data.scales.ls_vy = (STG_float*)malloc(sizeof(STG_float) * num);
    init_data.scales.ls_vz = (STG_float*)malloc(sizeof(STG_float) * num);
    init_data.scales.ls_wx = (STG_float*)malloc(sizeof(STG_float) * num);
    init_data.scales.ls_wy = (STG_float*)malloc(sizeof(STG_float) * num);
    init_data.scales.ls_wz = (STG_float*)malloc(sizeof(STG_float) * num);

    init_data.scales.ts_i = (STG_float*)malloc(sizeof(STG_float) * num);
    init_data.scales.ts_u = (STG_float*)malloc(sizeof(STG_float) * num);
    init_data.scales.ts_v = (STG_float*)malloc(sizeof(STG_float) * num);
    init_data.scales.ts_w = (STG_float*)malloc(sizeof(STG_float) * num);

    init_data.scales.ls_i[0] = 0.0019;
    init_data.scales.ts_i[0] = 1;

    STG_SmirnovData data;
    STG_float ts = 0.1;
    STG_int num_ts = 0;
    STG_VelNodeHist node_hist;

    STG_compute_Smirnov_data(init_data, num_modes, &data);

    STG_compute_Smirnov_node_hist(init_data, data, ts, num_ts, &node_hist, 0, 0, 0);

    printf("u = %.3f   v = %.3f  w = %.3f \n", node_hist.u_p[0], node_hist.v_p[0], node_hist.w_p[0]);
	printf("\n");

    STG_free_Smirnov_data(&data);
	STG_free_InitData(&init_data);
    STG_free_VelNodeHist(&node_hist);
}

int main(int argc, char * argv[])
{
	STG_init_rand();
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

	test_Smirnov_mom_field(5, 0.1, 1, 3, 2, 0, 0, 0, 0.01, 100);
	test_Smirnov_mom_field(5, 0.1, 1, 1, 1, 0, 0, 0, 0.01, 100);
	test_Smirnov_mom_field(5, 0.1, 32.346, 0.052, 0.747, -0.415, 0.351, -0.02, 0.01, 100);
	test_Smirnov_node_hist(32.346, 0.052, 0.747, -0.415, 0.351, -0.02, 100);
	test_Smirnov_mom_field(5, 0.1, 0.000001, 0.000001, 0.0000, 1.19e-8, -3.08e-7, 0.000000, 0.01, 100);

	test_Davidson_mom_filed(10, 0.1, 32.346, 0.052, 0.747, -0.415, 0.351, -0.02, 0.05, 100);
	test_Davidson_mom_filed(10, 10, 1, 2, 3, 0, 0, 0, 1, 100);
	test_Davidson_node_hist(0.01, 5, 0.03, 32.346, 0.052, 0.747, -0.415, 0.351, -0.02, 300);

	return 0;
}
