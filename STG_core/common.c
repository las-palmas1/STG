#include "common.h"
#include <stdlib.h>
#include <math.h>


// Функция возвращает массив заданного числа точек равномерно распределенных на заданном интервале
STG_float * get_uniform(STG_float min, STG_float max, STG_int num)
{
    STG_float * res = (STG_float*)malloc(sizeof (STG_float) * num);
    for (int i = 0; i < num; i++)
    {
        srand(time(0));
        res[i] = min + (STG_float)rand() / (STG_float)RAND_MAX * (max - min);
    }
    return res;
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

// Функция для нахождения с заданной точностью нуля заданной монотонной на заданном промежутке функции по алгоритму бисекции.
// Задается также максимальное число итераций.
STG_float bisection(STG_float x1, STG_float x2, STG_float (*mono_func)(STG_float), STG_float eps, STG_int max_iter)
{
    STG_float res;
    if (mono_func(x1) == 0.){
        res = mono_func(x1);
    }
    if (mono_func(x2) == 0.){
        res = mono_func(x2);
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
            if (sign(mono_func(x_inter)) != sign(mono_func(x_start))){
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
