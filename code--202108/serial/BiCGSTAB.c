/*******************************************************
    > File   : BiCGSTAB.c
    > Author : Yuntong
    > Mail   : 171840067@smail.nju.edu.cn
    > Date   : 2021/8/10
*******************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "BiCGSTAB.h"

void BICGSTAB(Stencil *A, double *diag, double *rhs, double *x,
        double EPSILON, int size, int maxit, int *numit) {
	int i, nits;
    int mesh_size = size * size;
    double alpha, beta, omega;
    double dot_res1, dot_res2, dot_d;
    double dot_s1, dot_s2;
    double residual, residual0, descent;
    double *tmp, *d, *s, *Ad, *As, *res;
    tmp  = (double *) calloc(mesh_size, sizeof(double));
    d    = (double *) calloc(mesh_size, sizeof(double));
    s    = (double *) calloc(mesh_size, sizeof(double));
    Ad   = (double *) calloc(mesh_size, sizeof(double));
    As   = (double *) calloc(mesh_size, sizeof(double));
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
        
        //update s
        dot_res1 = dot_res2;
        dot_d = dot(Ad, rhs, mesh_size);
        alpha = dot_res1 / dot_d;
        const_vector_multiply(alpha, Ad, tmp, mesh_size);
        vector_subtraction(res, tmp, s, mesh_size);

        // Termination criterion
        residual = dot(s, s, mesh_size);
        descent  = sqrt(residual / residual0);
        if ((nits % 50 == 0) | (nits == 1)) printf("%d : \t%.8e\n", nits, descent);
        if (descent < EPSILON) {
            printf("%d : \t%.8e\n", nits, descent);
            break;
        }

        // As = (A - D) * s + D * s
        MV_stencil(A, s, As, size);
        for (i = 0; i < mesh_size; i++) {
            As[i] += diag[i] * s[i];
        }

        //update x
        dot_s1 = dot(s,  As,  mesh_size);
        dot_s2 = dot(As, As,  mesh_size);
        omega  = dot_s1  / dot_s2;
        const_vector_multiply(alpha, d, tmp, mesh_size);
        vector_addition(x, tmp, x, mesh_size);
        const_vector_multiply(omega, s, tmp, mesh_size);
        vector_addition(x, tmp, x, mesh_size);
        
        //update residul
        const_vector_multiply(omega, As, tmp, mesh_size);
        vector_subtraction(s, tmp, res, mesh_size);
        dot_res2 = dot(res, rhs, mesh_size);

        //update direction
        beta = (dot_res2 / dot_res1) * (alpha / omega);
        const_vector_multiply(omega, Ad, tmp, mesh_size);
        vector_subtraction(d, tmp, tmp, mesh_size);
        const_vector_multiply(beta, tmp, tmp, mesh_size);
        vector_addition(res, tmp, d, mesh_size);
        
    }

    printf("\n");
    *numit = nits;
    
    free(tmp);  free(d);   free(s);   
    free(res);  free(Ad);  free(As);
    
    return ;
}

void P_BICGSTAB(Stencil *A, double *diag, double *rhs, double *x, int precond, int precond_iter,
        double EPSILON, int size, int maxit, int *numit) {
	int i, nits;
    int mesh_size = size * size;
    double alpha, beta, omega;
    double dot_res1, dot_res2, dot_d;
    double dot_s1, dot_s2;
    double residual, residual0, descent;
    double *tmp, *d, *s, *y, *z, *Ay, *Az, *res;
    tmp  = (double *) calloc(mesh_size, sizeof(double));
    d    = (double *) calloc(mesh_size, sizeof(double));
    s    = (double *) calloc(mesh_size, sizeof(double));
    y    = (double *) calloc(mesh_size, sizeof(double));
    z    = (double *) calloc(mesh_size, sizeof(double));
    Ay   = (double *) calloc(mesh_size, sizeof(double));
    Az   = (double *) calloc(mesh_size, sizeof(double));
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
    	// y = M^{-1}*d
	    switch (precond) {
	        case 1: Precond_JAC(A, diag, d, y, EPSILON, size, precond_iter); break;
            case 2: Precond_GS(A, diag, d, y, EPSILON, size, precond_iter); break;
	        case 3: Precond_CG(A, diag, d, y, EPSILON, size, precond_iter);  break;
	    }
        MV_stencil(A, y, Ay, size);
        for (i = 0; i < mesh_size; i++) {
            Ay[i] += diag[i] * y[i];
        }
        
        //update s
        dot_res1 = dot_res2;
        dot_d = dot(Ay, rhs, mesh_size);
        alpha = dot_res1 / dot_d;
        const_vector_multiply(alpha, Ay, tmp, mesh_size);
        vector_subtraction(res, tmp, s, mesh_size);

        // Termination criterion
        residual = dot(s, s, mesh_size);
        descent  = sqrt(residual / residual0);
        if ((nits % 20 == 0) | (nits == 1)) printf("%d : \t%.8e\n", nits, descent);
        if (descent < EPSILON) {
            printf("%d : \t%.8e\n", nits, descent);
            break;
        }

        // z = M^{-1}*s
	    switch (precond) {
	        case 1: Precond_JAC(A, diag, s, z, EPSILON, size, precond_iter); break;
	        case 2: Precond_GS(A, diag, s, z, EPSILON, size, precond_iter);  break;
            case 3: Precond_CG(A, diag, s, z, EPSILON, size, precond_iter);  break;
	    }
        MV_stencil(A, z, Az, size);
        for (i = 0; i < mesh_size; i++) {
            Az[i] += diag[i] * z[i];
        }

        //update x
        dot_s1 = dot(s,  Az,  mesh_size);
        dot_s2 = dot(Az, Az,  mesh_size);
        omega = dot_s1 	 / dot_s2;
        const_vector_multiply(alpha, y, tmp, mesh_size);
        vector_addition(x, tmp, x, mesh_size);
        const_vector_multiply(omega, z, tmp, mesh_size);
        vector_addition(x, tmp, x, mesh_size);
        
        //update residul
        const_vector_multiply(omega, Az, tmp, mesh_size);
        vector_subtraction(s, tmp, res, mesh_size);
        dot_res2 = dot(res, rhs, mesh_size);

        //update direction
        beta = (dot_res2 / dot_res1) * (alpha / omega);
        const_vector_multiply(omega, Ay, tmp, mesh_size);
        vector_subtraction(d, tmp, tmp, mesh_size);
        const_vector_multiply(beta, tmp, tmp, mesh_size);
        vector_addition(res, tmp, d, mesh_size);
        
    }

    printf("\n");
    *numit = nits;
    
    free(d);    free(s);    free(y);   free(z);   
    free(tmp);  free(res);  free(Ay);  free(Az);
    
    return ;
}




