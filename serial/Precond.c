/*******************************************************
    > File   : Precond.c
    > Author : Yuntong
    > Mail   : 171840067@smail.nju.edu.cn
    > Date   : 
*******************************************************/

#include <stdlib.h>
#include <math.h>
#include "Precond.h"


void Precond_JAC(Stencil *A, double *diag, double *rhs, double *x,
        double EPSILON, int size, int maxit) {
    int i, nits;
    int mesh_size = size * size;
    double residual0;
    double *tmp, *c;
    Stencil *D_A;
    c   = (double *)  calloc(mesh_size, sizeof(double));
    tmp = (double *)  calloc(mesh_size, sizeof(double));
    D_A = (Stencil *) calloc(mesh_size, sizeof(Stencil));

    
    for (i = 0; i < mesh_size; i++) {
        D_A[i].value[0] = A[i].value[0] / diag[i];    // D^(-1) * (A - D)
        D_A[i].value[1] = A[i].value[1] / diag[i];
        D_A[i].value[2] = A[i].value[2] / diag[i];
        D_A[i].value[3] = A[i].value[3] / diag[i];
          c[i] 		 	= rhs[i] 		/ diag[i];    // c = D^(-1) * b
          x[i]			= c[i];
    }
    residual0 = dot(rhs, rhs, mesh_size);
    
    for (nits = 2; nits <= maxit; nits++) {
        // tmp = D^(-1)*(A-D) * x_old
        MV_stencil(D_A, x, tmp, size);
        // x_new = c - tmp
        vector_subtraction(c, tmp, x, mesh_size);
    }

    free(tmp);  free(c);  free(D_A);
    return ;
}

void Precond_CG(Stencil *A, double *diag, double *rhs, double *x,
        double EPSILON, int size, int maxit) {
    int i, nits;
    int mesh_size = size * size;
    double alpha, beta;
    double dot_res1, dot_res2, dot_d;
    double residual0;
    double *tmp, *d, *Ad, *res;
    tmp  = (double *) calloc(mesh_size, sizeof(double));
    d    = (double *) calloc(mesh_size, sizeof(double));
    Ad   = (double *) calloc(mesh_size, sizeof(double));
    res  = (double *) calloc(mesh_size, sizeof(double));
    
    //d_0 = r_0 = b
    for (i = 0; i < mesh_size; i++) {
        d[i]   = rhs[i];
        res[i] = rhs[i];
        x[i]   = 0;
    }
    residual0 = dot(res, res, mesh_size);
    dot_res2 = residual0;

    for (nits = 1; nits <= maxit; nits++) {
        // Ad = (A - D) * d + D * d
        MV_stencil(A, d, Ad, size);
        for (i = 0; i < mesh_size; i++) {
            Ad[i] += diag[i] * d[i];
        }

        //inner product
        dot_res1 = dot_res2;
        dot_d    = dot(d, Ad, mesh_size);

        //compute x_new
        alpha = dot_res1 / dot_d;
        const_vector_multiply(alpha, d, tmp, mesh_size);
        vector_addition(x, tmp, x, mesh_size);
        
        //update residul
        const_vector_multiply(alpha, Ad, Ad, mesh_size);
        vector_subtraction(res, Ad, res, mesh_size);
        dot_res2 = dot(res, res, mesh_size);

        //update direction
        beta = dot_res2 / dot_res1;
        const_vector_multiply(beta, d, d, mesh_size);
        vector_addition(res, d, d, mesh_size); 
    }
    
    free(tmp);  free(d);  free(Ad);  free(res);
    
    return ;
}

