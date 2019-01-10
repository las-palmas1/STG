#include "common.h"
#include <stdio.h>
#include <stdlib.h>


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

static void test_eig_values_computing(STG_float m11, STG_float m22, STG_float m33, STG_float m12, STG_float m13, STG_float m23)
{
	printf("Eigvalues test.\n");
	printf("Matrix: \n");
	printf("%.3f  %.3f  %.3f \n", m11, m12, m13);
	printf("%.3f  %.3f  %.3f \n", m12, m22, m23);
	printf("%.3f  %.3f  %.3f \n", m13, m23, m33);
	STG_float * eig = eig_3x3_sym(m11, m22, m33, m12, m13, m23);

	STG_float M1[3][3] = {
		{ m11 - eig[0], m12, m13 },
		{ m12, m22 - eig[1], m23 },
		{ m13, m23, m33 - eig[2] },
	};
	STG_float det = (
		M1[0][0] * M1[1][1] * M1[2][2] - M1[0][0] * M1[1][2] * M1[2][1] - M1[0][1] * M1[1][0] * M1[2][2] +
		M1[0][1] * M1[1][2] * M1[2][0] + M1[0][2] * M1[1][0] * M1[2][1] - M1[0][2] * M1[1][1] * M1[2][0]
	);
	printf("\n");
	printf("eig1 = %.3f\n", eig[0]);
	printf("eig2 = %.3f\n", eig[1]);
	printf("eig3 = %.3f\n", eig[2]);
	printf("det(M - eig*I) = %.4f\n", det);
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

int main(int argc, char * argv[])
{

    test_bisection();
	test_uniform(20, -5, 5);
	test_normal(20, 0, 1);
	test_trigon(20);

	test_eig_values_computing(1, 1, 1, 1, 1, 1);
	test_eig_values_computing(1, 1, 1, 2, 1, 1);
	test_eig_values_computing(1, 1, 3, 3, 1, 4);
	test_eig_values_computing(2, 3, 9, 0, 0, 4);

	return 0;
}
