#include "precompiled.h"
#include "common.h"
#include <stdlib.h>
#include <math.h>


// Функция возвращает массив заданного числа точек равномерно распределенных на заданном интервале
STG_float * get_uniform(STG_float min, STG_float max, STG_int num)
{
    STG_float * res = (STG_float*)malloc(sizeof (STG_float) * num);
	srand(time(NULL));
	for (int i = 0; i < num; i++)
    {
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
	srand(time(NULL));
	for (STG_int i = 0; i < num; i++)
	{
		res[i] = mu + sigma * box_muller();
	}
	return res;
}

// Рассчитывает массив случайных чисел с распределением вероятностей p(x) = 1/2 * sin(x) по методу обратного преобразования 
STG_float * get_trigon(STG_int num)
{
	STG_float * res = (STG_float*)malloc(sizeof(STG_float) * num);
	srand(time(NULL));
	for (STG_int i = 0; i < num; i++)
	{
		STG_float u = (STG_float)rand() / (STG_float)RAND_MAX;
		res[i] = acos(1 - 2 * u);
	}
	return res;
}