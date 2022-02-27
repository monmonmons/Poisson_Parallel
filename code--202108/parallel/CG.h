/*******************************************************
    > File   : CG.h
    > Author : Yuntong
    > Mail   : 171840067@smail.nju.edu.cn
    > Date   : 2021/4/13
*******************************************************/

#ifndef CG_h
#define CG_h

#include "Precond.h"

//Conjugate gradient method
//  Input    : A, diag, b, x, EPSILON, size, numit
//  Modified : A, x, numit
void CG(Stencil *A, double *diag, double *rhs, double *x,
        double EPSILON, int maxit, int *numit,
        int nlx, int nly, int myidx, int myidy, int npx, int npy);

//Preconditioned conjugate gradient method
//  Preconditoner : D^(-1)

void PCG(Stencil *A, double *diag, double *rhs, double *x, 
        int precond, int precond_iter, double EPSILON, int maxit, int *numit,
        int nlx, int nly, int myidx, int myidy, int npx, int npy);

#endif /* CG_h */
