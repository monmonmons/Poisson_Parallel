/*******************************************************
    > File   : Precond.h
    > Author : Yuntong
    > Mail   : 171840067@smail.nju.edu.cn
    > Date   : 
*******************************************************/

#ifndef Precond_h
#define Precond_h

#include "stencil.h"

void Precond_JAC(Stencil *A, double *diag, double *rhs, double *x,
        double EPSILON, int nlx, int nly, int maxit);

void Precond_CG(Stencil *A, double *diag, double *rhs, double *x,
        double EPSILON, int nlx, int nly, int maxit);

#endif /* Precond_h */

