/*******************************************************
    > File   : Precond.h
    > Author : Yuntong
    > Mail   : 171840067@smail.nju.edu.cn
    > Date   : 2021/8/15
*******************************************************/

#ifndef Precond_h
#define Precond_h

#include "ops.h"

void Precond_JAC(Stencil *A, double *diag, double *rhs, double *x,
        double EPSILON, int size, int maxit);

void Precond_GS(Stencil *A, double *diag, double *rhs, double *x,
        double EPSILON, int size, int maxit);

void Precond_CG(Stencil *A, double *diag, double *rhs, double *x,
        double EPSILON, int size, int maxit);

#endif /* Precond_h */

