/*******************************************************
    > File   : CG.h
    > Author : Yuntong
    > Mail   : 171840067@smail.nju.edu.cn
    > Date   : 2021/3/23
*******************************************************/

#ifndef CG_h
#define CG_h

#include "Precond.h"

//Conjugate gradient method
//  Input    : A, diag, b, x, EPSILON, size, numit
//  Modified : A, x, numit
void CG(Stencil *A, double *diag, double *rhs, double *x,
        double EPSILON, int size, int maxit, int *numit);

//Preconditioned conjugate gradient method
//  Preconditioned -- left
void PCG(Stencil *A, double *diag, double *rhs, double *x, int precond, int precond_iter,
        double EPSILON, int size, int maxit, int *numit);

// //  Preconditioned -- split
// void SPCG(Stencil *A, double *diag, double *b, double *x,
//         double EPSILON, int size, int maxit, int *numit);

#endif /* CG_h */
