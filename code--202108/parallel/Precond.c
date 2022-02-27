/*******************************************************
    > File   : Precond.c
    > Author : Yuntong
    > Mail   : 171840067@smail.nju.edu.cn
    > Date   : 2021/8/25
*******************************************************/

#include <stdlib.h>
#include <math.h>
#include "Precond.h"

void Precond_JAC(Stencil *A, double *diag, double *rhs, double *x,
        double EPSILON, int nlx, int nly, int maxit) {
    int i, nits;
    int mesh_size = nly * nlx;
    double *tmp, *c;
    Stencil *D_A;
    tmp = (double *)  calloc(mesh_size, sizeof(double));
    c   = (double *)  calloc(mesh_size, sizeof(double));
    D_A = (Stencil *) calloc(mesh_size, sizeof(Stencil));

    
    for (i = 0; i < mesh_size; i++) {
        D_A[i].value[0] = A[i].value[0] / diag[i];    // D^(-1) * (A - D)
        D_A[i].value[1] = A[i].value[1] / diag[i];
        D_A[i].value[2] = A[i].value[2] / diag[i];
        D_A[i].value[3] = A[i].value[3] / diag[i];
          c[i] 		 	= rhs[i] 		/ diag[i];    // c = D^(-1) * b
          x[i]			= c[i];
    }
    
    for (nits = 2; nits <= maxit; nits++) {
        // tmp = D^(-1)*(A-D) * x_old
        MV_stencil_zero(D_A, x, tmp, nlx, nly);
        // x_new = c - tmp
        vector_subtraction(c, tmp, x, mesh_size);
    }

    
    free(tmp);  free(c);   free(D_A);
    return ;
}

void Precond_GS(Stencil *A, double *diag, double *rhs, double *x,
        double EPSILON, int nlx, int nly, int maxit) {
    int i, nits;
    int mesh_size = nly * nlx;
    double *tmp, *c, *vec;
    Stencil *D_A;
    tmp = (double *)  calloc(mesh_size, sizeof(double));
    c   = (double *)  calloc(mesh_size, sizeof(double));
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
    
    for (nits = 2; nits <= maxit; nits++) {
        //x_new = (D - L)^(-1) * U * x_old + (D - L)^(-1) * b
        for (i = 0; i < mesh_size; i++) {
            vec[0] = ((i % nlx) == 0) ? 0 : x[i - 1];              //west
            vec[1] = ((i % nlx) == nlx - 1) ? 0 : tmp[i + 1];     //east
            vec[2] = ((i / nlx) == nlx - 1) ? 0 : tmp[i + nlx];  //north
            vec[3] = ((i / nlx) == 0) ? 0 : x[i - nlx];           //south
            x[i] = c[i] - dot(D_A[i].value, vec, 4);
        }
        //x_old = x_new
        for (i = 0; i < mesh_size; i++) {
            tmp[i] = x[i];
        }
    }

    
    free(tmp);  free(c);   free(vec);    free(D_A);
    return ;
}

void Precond_CG(Stencil *A, double *diag, double *rhs, double *x,
        double EPSILON, int nlx, int nly, int maxit) {
    int i, nits;
    int mesh_size = nly * nlx;
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
        MV_stencil_zero(A, d, Ad, nlx, nly);
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
