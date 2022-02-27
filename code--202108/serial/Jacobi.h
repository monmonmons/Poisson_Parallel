/*******************************************************
    > File   : Jacobi.h
    > Author : Yuntong
    > Mail   : 171840067@smail.nju.edu.cn
    > Date   : 2021/8/20
*******************************************************/

#ifndef Jacobi_h
#define Jacobi_h

#include "ops.h"

//Jacobi method : compute x_new = D^(-1) * (D - A) * x_old + D^(-1) * b
//  Input 	 : A, diag, b, x, EPSILON, size, maxit, numit
//  Modified : A, x, numit
void Jacobi(Stencil *A, double *diag, double *rhs, double *x,
                    double EPSILON, int size, int maxit, int *numit);

//Gauss-Seidel method : compute x_new = (D - L)^(-1) * U * x_old + (D - L)^(-1) * b
//  Input    : A, diag, b, x, EPSILON, size, maxit, numit
//  Modified : A, x, numit
void GS(Stencil *A, double *diag, double *rhs, double *x,
                    double EPSILON, int size, int maxit, int *numit);

#endif /* Jacobi_h */
