/*******************************************************
    > File   : stencil.c
    > Author : Yuntong
    > Mail   : 171840067@smail.nju.edu.cn
    > Date   : 2021/3/30
*******************************************************/

#ifndef stencil_h
#define stencil_h

#include "BLAS.h"

//Data structure
typedef struct stencil {
    int number;
    double value[4]; //0-west 1-east 2-north 3-south
} Stencil;

//Initialize
//  Input : stencil, number, 4 values
//  Modified : stencil
void initial_stencil(Stencil *stencil, int num, double v1, double v2, double v3, double v4);

// void print_stencil(Stencil *stencil);

void print_line(Stencil *stencil);

//Matrix-Vector Multiply
//  Input: stencil, 4 boundary vector, vec, result, nlx, nly
//  Modified : result
void MV_stencil(Stencil *matrix, double *vec, double *west, double *east, double *north, double *south, 
                double *result, int nlx, int nly);

void MV_stencil_zero(Stencil *matrix, double *vec, double *result, int nlx, int nly);
#endif /* stencil_h */
