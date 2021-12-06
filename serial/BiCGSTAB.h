/*******************************************************
    > File   : BiCGSTAB.h
    > Author : Yuntong
    > Mail   : 171840067@smail.nju.edu.cn
    > Date   : 2021/4/15
*******************************************************/

#ifndef BiCGSTAB_h
#define BiCGSTAB_h

#include "Precond.h"

//BI-Conjugate gradient method
//  Input    : A, diag, b, x, EPSILON, size, numit
//  Modified : A, x, numit
void BICGSTAB(Stencil *A, double *diag, double *rhs, double *x,
        double EPSILON, int size, int maxit, int *numit);

//Preconditioned conjugate gradient method
//  Preconditioned -- right
void P_BICGSTAB(Stencil *A, double *diag, double *rhs, double *x, int precond, int precond_iter,
        double EPSILON, int size, int maxit, int *numit);

#endif /* BiCGSTAB_h */

