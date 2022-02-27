/*******************************************************
    > File   : CG.c
    > Author : Yuntong
    > Mail   : 171840067@smail.nju.edu.cn
    > Date   : 2021/8/10
*******************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "CG.h"

void CG(Stencil *A, double *diag, double *rhs, double *x,
        double EPSILON, int size, int maxit, int *numit) {
    int i, nits;
    int mesh_size = size * size;
    double alpha, beta;
    double dot_res1, dot_res2, dot_d;
    double residual, residual0, descent;
    double *tmp, *d, *Ad, *res;
    tmp  = (double *) calloc(mesh_size, sizeof(double));
    d    = (double *) calloc(mesh_size, sizeof(double));
    Ad   = (double *) calloc(mesh_size, sizeof(double));
    res  = (double *) calloc(mesh_size, sizeof(double));
        
    //x_0
    for (i = 0; i < mesh_size; i++)  x[i] = 0;

    //d_0 = r_0 = b - A * x_0
    MV_stencil(A, x, tmp, size);
    for (i = 0; i < mesh_size; i++) {
        tmp[i] += diag[i] * x[i];
    }
    vector_subtraction(rhs, tmp, res, mesh_size);
    for (i = 0; i < mesh_size; i++)  d[i] = res[i];

    residual0 = dot(res, res, mesh_size);
    printf("0 : \t%.8e\n", residual0);
    printf("------------------------------------------------------\n");
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

        
        // Termination criterion
        residual = dot_res2;
        descent = sqrt(residual / residual0);
        if (nits % 20 == 0 | nits == 1) printf("%d : \t%.8e\n", nits, descent);
        if (descent < EPSILON) {
            printf("%d : \t%.8e\n", nits, descent);
            break;
        }
        
    }

    printf("\n");
    *numit = nits;
    
    free(tmp);  free(d);  free(Ad);  free(res);
    
    return ;
}


void PCG(Stencil *A, double *diag, double *rhs, double *x, int precond, int precond_iter,
                double EPSILON, int size, int maxit, int *numit) {
    int i, nits;
    int mesh_size = size * size;
    double alpha, beta;
    double dot_res1, dot_res2, dot_d;
    double residual, residual0, descent;
    double *tmp, *d, *Ad, *res, *res_p;
    tmp = (double *) calloc(mesh_size, sizeof(double));
    d   = (double *) calloc(mesh_size, sizeof(double));
    Ad  = (double *) calloc(mesh_size, sizeof(double));
    res = (double *) calloc(mesh_size, sizeof(double));
    res_p = (double *) calloc(mesh_size, sizeof(double));
    
    //Preconditioner : M^{-1}
    // z_0 = M^{-1} * rhs;
    switch (precond) {
        case 1: Precond_JAC(A, diag, rhs, res_p, EPSILON, size, precond_iter); break;
        case 2: Precond_GS(A, diag, rhs, res_p, EPSILON, size, precond_iter); break;
        case 3: Precond_CG(A, diag, rhs, res_p, EPSILON, size, precond_iter);  break;
    }

    //x_0
    for (i = 0; i < mesh_size; i++)  x[i] = 0;

    //r_0 = b - A * x_0; d_0 = z_0
    MV_stencil(A, x, tmp, size);
    for (i = 0; i < mesh_size; i++) {
        tmp[i] += diag[i] * x[i];
    }
    vector_subtraction(rhs, tmp, res, mesh_size);
    for (i = 0; i < mesh_size; i++)  d[i] = res_p[i];

    residual0 = dot(res, res, mesh_size);
    printf("0 : \t%.8e\n", residual0);
    printf("------------------------------------------------------\n");
    dot_res2 = dot(res, res_p, mesh_size);
    
    for (nits = 1; nits <= maxit; nits++) {
        // Ad = (A - D) * d + D * d
        MV_stencil(A, d, Ad, size);
        for (i = 0; i < mesh_size; i++) {
            Ad[i] += diag[i] * d[i];
        }
        
        //inner product
        dot_res1 = dot_res2;
        dot_d    = dot(d, Ad, mesh_size);

        //update x
        // x_j = x_{j-1} + alpha * d
        alpha = dot_res1 / dot_d;
        const_vector_multiply(alpha, d, tmp, mesh_size);
        vector_addition(x, tmp, x, mesh_size);
        
        //update residul
        // r_j = r_{j-1} - alpha * A * d
        const_vector_multiply(alpha, Ad, Ad, mesh_size);
        vector_subtraction(res, Ad, res, mesh_size);
        // z_j = M^(-1) * r_j
        switch (precond) {
            case 1: Precond_JAC(A, diag, res, res_p, EPSILON, size, precond_iter); break;
            case 2: Precond_CG(A, diag, res, res_p, EPSILON, size, precond_iter);  break;
        }
        dot_res2 = dot(res, res_p, mesh_size);
        
        //update direction
        beta = dot_res2 / dot_res1;
        const_vector_multiply(beta, d, d, mesh_size);
        vector_addition(res_p, d, d, mesh_size);
        
        // Termination criterion
        residual = dot(res, res, mesh_size);
        descent = sqrt(residual / residual0);
        if (nits % 20 == 0 | nits == 1) printf("%d : \t%.8e\n", nits, descent);
        if (descent < EPSILON) {
            printf("%d : \t%.8e\n", nits, descent);
            break;
        }
        
    }

    printf("\n");
    *numit = nits;
    
    
    free(tmp);  free(d);  free(Ad);  free(res); free(res_p);
    
    return ;
}



