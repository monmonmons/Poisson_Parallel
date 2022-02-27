/*******************************************************
    > File   : ops.h
    > Author : Yuntong
    > Mail   : 171840067@smail.nju.edu.cn
    > Date   : 2021/8/20
*******************************************************/

#ifndef OPS_h
#define OPS_h

#include "stencil.h"

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

//Matrix-Vector Multiply
//  Input: stencil, 4 boundary vector, vec, result, nlx, nly
//  Modified : result
void MV_stencil(Stencil *matrix, double *vec, double *west, double *east, double *north, double *south, 
                double *result, int nlx, int nly);

void MV_stencil_zero(Stencil *matrix, double *vec, double *result, int nlx, int nly);

#endif /* OPS_h */
