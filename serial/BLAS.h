/*******************************************************
    > File   : BLAS.h
    > Author : Yuntong
    > Mail   : 171840067@smail.nju.edu.cn
    > Date   : 2021/3/23
*******************************************************/

#ifndef BLAS_h
#define BLAS_h

//Compute infinity-norm of x
//  Input  : x, length
//  Output : result = ||x||_inf
double norm_inf(double *x, int length);

//Compute 2-norm of x
//  Input  : x, length
//  Output : result = ||x||_2
double norm_2(double *x, int length);

//Compute inner product of two vectors
//  Input : x, y, length
//  Output : result = < x, y >
double dot(double *x, double *y, int length);

void vector_addition(double *x, double *y, double *result, int length);

void vector_subtraction(double *x, double *y, double *result, int length);

void const_vector_multiply(double c, double *x, double *result, int length);


double average(double x, double y);

#endif /* BLAS_h */
