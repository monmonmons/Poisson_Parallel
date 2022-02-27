/*******************************************************
    > File   : Jacobi.c
    > Author : Yuntong
    > Mail   : 171840067@smail.nju.edu.cn
    > Date   : 2021/8/20
*******************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Jacobi.h"


void Jacobi(Stencil *A, double *diag, double *rhs, double *x,
                   double EPSILON, int size, int maxit, int *numit) {
    int i, nits;
    int mesh_size = size * size;
    double residual, residual0, descent;
    double *tmp, *c, *res;
    Stencil *D_A;
    tmp = (double *)  calloc(mesh_size, sizeof(double));
    c   = (double *)  calloc(mesh_size, sizeof(double));
    res = (double *)  calloc(mesh_size, sizeof(double));
    D_A = (Stencil *) calloc(mesh_size, sizeof(Stencil));
    
    for (i = 0; i < mesh_size; i++) {
        D_A[i].value[0] = A[i].value[0] / diag[i];    // D^(-1) * (A - D)
        D_A[i].value[1] = A[i].value[1] / diag[i];
        D_A[i].value[2] = A[i].value[2] / diag[i];
        D_A[i].value[3] = A[i].value[3] / diag[i];
          c[i]          = rhs[i]        / diag[i];    // c = D^(-1) * b
          x[i]          = 0;
    }
    residual0 = dot(rhs, rhs, mesh_size);
    printf("0 : \t%.8e\n", residual0);
    printf("------------------------------------------------------\n");
    
    for (nits = 1; nits <= maxit; nits++) {
        // tmp = D^(-1)*(A-D) * x_old
        MV_stencil(D_A, x, tmp, size);
        // x_new = c - tmp
        vector_subtraction(c, tmp, x, mesh_size);
        
        //residual
        //  r = b - D * [D^(-1)*(A-D) * x_new + x_new]
        MV_stencil(D_A, x, res, size);
        vector_addition(res, x, res, mesh_size);
        for (i = 0; i < mesh_size; i++) {
            res[i] = diag[i] * res[i];
        }
        vector_subtraction(rhs, res, res, mesh_size);
        residual = dot(res, res, mesh_size);
        
        // Termination criterion
        descent = sqrt(residual / residual0);
        if (nits % 30 == 0 | nits == 1) printf("%d : \t%.8e\n", nits, descent);
        // printf("%d : \t%.8e\n", nits, descent);
        if (descent < EPSILON) {
            printf("%d : \t%.8e\n", nits, descent);
            break;
        }
    }
    
    printf("\n");
    *numit = nits;
    
    free(tmp);  free(c);  free(res);  free(D_A);
    return ;
}


void GS(Stencil *A, double *diag, double *rhs, double *x,
                    double EPSILON, int size, int maxit, int *numit) {
    int i, nits;
    int mesh_size = size * size;
    double residual, residual0, descent;
    double *tmp, *c, *res, *vec;
    Stencil *D_A;
    c   = (double *)  calloc(mesh_size, sizeof(double));
    tmp = (double *)  calloc(mesh_size, sizeof(double));
    res = (double *)  calloc(mesh_size, sizeof(double));
    vec = (double *)  calloc(4,         sizeof(double));
    D_A = (Stencil *) calloc(mesh_size, sizeof(Stencil));
    
    for (i = 0; i < mesh_size; i++) {
        D_A[i].value[0] = A[i].value[0] / diag[i];    // D^(-1) * (A - D)
        D_A[i].value[1] = A[i].value[1] / diag[i];
        D_A[i].value[2] = A[i].value[2] / diag[i];
        D_A[i].value[3] = A[i].value[3] / diag[i];
          c[i]          = rhs[i]        / diag[i];    // c = D^(-1) * b
          x[i]          = 0;
          tmp[i]        = 0;
    }
    residual0 = dot(rhs, rhs, mesh_size);
    printf("0 : \t%.8e\n", residual0);
    printf("------------------------------------------------------\n");
    
    for (nits = 1; nits <= maxit; nits++) {
        //x_new = (D - L)^(-1) * U * x_old + (D - L)^(-1) * b
        for (i = 0; i < mesh_size; i++) {
            vec[0] = ((i % size) == 0) ? 0 : x[i - 1];              //west
            vec[1] = ((i % size) == size - 1) ? 0 : tmp[i + 1];     //east
            vec[2] = ((i / size) == size - 1) ? 0 : tmp[i + size];  //north
            vec[3] = ((i / size) == 0) ? 0 : x[i - size];           //south
            x[i]   = c[i] - dot(D_A[i].value, vec, 4);
        }
        
        //residual
        //  r = b - D * [D^(-1)*(A-D) * x_new + x_new]
        MV_stencil(D_A, x, res, size);
        vector_addition(res, x, res, mesh_size);
        for (i = 0; i < mesh_size; i++) {
            res[i] = diag[i] * res[i];
        }
        vector_subtraction(rhs, res, res, mesh_size);
        residual = dot(res, res, mesh_size);
        
        // Termination criterion
        descent = sqrt(residual / residual0);
        if (nits % 30 == 0 | nits == 1) printf("%d : \t%.8e\n", nits, descent);
        // printf("%d : \t%.8e\n", nits, descent);
        if (descent < EPSILON) {
            printf("%d : \t%.8e\n", nits, descent);
            break;
        }

        //x_old = x_new
        for (i = 0; i < mesh_size; i++) {
            tmp[i] = x[i];
        }
    }
    
    printf("\n");
    *numit = nits;
    
    free(tmp);   free(c);  free(res);  free(vec);  free(D_A);
    return ;
}
