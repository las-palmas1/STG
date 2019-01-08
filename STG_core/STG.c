#include "common.h"
#include <stdio.h>


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
        STG_float x1, STG_float x2, STG_float (*func)(STG_float),
        STG_float eps, STG_int max_iter, char * func_name, STG_float expect_res
    )
{
    printf("Test bisection. Function:     %s\n", func_name);
    printf("Inteval:   x1 = %.3f,   x2 = %.3f\n", x1, x2);
    printf("eps = %.4f,  max_iter = %d\n", eps, max_iter);

    STG_float res = bisection(x1, x2, func, eps, max_iter);

    printf("Result:   %.4f\n", res);
    printf("Expected:   %.4f\n", expect_res);
    printf("\n");

}

static void test_bisection()
{
    STG_float eps = 0.0001;
    STG_float max_iter = 100;

    bisection_wrap(0, 2, test_func_linear, eps, max_iter, "f(x) = x - 1", 1);

    bisection_wrap(-8, -1, test_func_pow2, eps, max_iter, "f(x) = (x - 2) * (x + 4)", -4);
    bisection_wrap(-1, 5, test_func_pow2, eps, max_iter, "f(x) = (x - 2) * (x + 4)", -2);

    bisection_wrap(-9, -3, test_func_pow3, eps, max_iter, "f(x) = (x - 1) * (x - 8) * (x + 6)", -6);
    bisection_wrap(-3, 5, test_func_pow3, eps, max_iter, "f(x) = (x - 1) * (x - 8) * (x + 6)", 1);
    bisection_wrap(5, 12, test_func_pow3, eps, max_iter, "f(x) = (x - 1) * (x - 8) * (x + 6)", 8);

}

int main(int argc, char* argv[])
{

    test_bisection();

	return 0;
}
