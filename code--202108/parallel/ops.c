/*******************************************************
    > File   : ops.c
    > Author : Yuntong
    > Mail   : 171840067@smail.nju.edu.cn
    > Date   : 2021/8/20
*******************************************************/

#include <math.h>
#include "ops.h"

double norm_inf(double *x, int length) {
    int i;
    double max = fabs(x[0]);
    for (i = 1; i < length; i++) {
        if (max < fabs(x[i])) {
            max = fabs(x[i]);
        }
    }
    return max;
}

double norm_2(double *x, int length) {
    int i;
    double norm = 0;
    for (i = 0; i < length; i++) {
        norm += pow(x[i], 2);
    }
    norm = sqrt(norm);
    
    return norm;
}

double dot(double *x, double *y, int length) {
    int i;
    double result = 0;
    for (i = 0; i < length; i++) {
        result += x[i] * y[i];
    }
    return result;
}

void vector_addition(double *x, double *y, double *result, int length) {
    int i;
    double temp;
    for (i = 0; i < length; i++) {
        temp = x[i] + y[i];
        result[i] = temp;
    }
    
    return ;
}

void vector_subtraction(double *x, double *y, double *result, int length) {
    int i;
    double temp;
    for (i = 0; i < length; i++) {
        temp = x[i] - y[i];
        result[i] = temp;
    }
    
    return ;
}

void const_vector_multiply(double c, double *x, double *result, int length) {
    int i;
    double temp;
    for (i = 0; i < length; i++) {
        temp = c * x[i];
        result[i] = temp;
    }
    
    return ;
}

//algebra
double average(double x, double y) {
    double a1 = 0.5 * (x + y);
    double a2 = 2 * (x * y) / (x + y);
    return a1 < a2 ? a1 : a2;
}


void MV_stencil(Stencil *matrix, double *vec, double *west, double *east, double *north, double *south, 
                double *result, int nlx, int nly) {
    int i, j, num;
    double temp[4];
    for (j = 0; j < nly; j++) {
        for (i = 0; i < nlx; i++) {
            num = j * nlx + i;
            temp[0] = (i == 0) ? west[j] : vec[num - 1];             //west
            temp[1] = (i == nlx - 1) ? east[j] : vec[num + 1];       //east
            temp[2] = (j == nly - 1) ? north[i] : vec[num + nlx];    //north
            temp[3] = (j == 0) ? south[i] : vec[num - nlx];          //south
            
            result[num] = dot(matrix[num].value, temp, 4);
        }
    }
    return ;
}

void MV_stencil_zero(Stencil *matrix, double *vec, double *result, int nlx, int nly) {
    int i, j, num;
    double temp[4];
    for (j = 0; j < nly; j++) {
        for (i = 0; i < nlx; i++) {
            num = j * nlx + i;
            temp[0] = (i == 0) ? 0 : vec[num - 1];             //west
            temp[1] = (i == nlx - 1) ? 0 : vec[num + 1];       //east
            temp[2] = (j == nly - 1) ? 0 : vec[num + nlx];     //north
            temp[3] = (j == 0) ? 0 : vec[num - nlx];           //south
            
            result[num] = dot(matrix[num].value, temp, 4);
        }
    }
    return ;
}
