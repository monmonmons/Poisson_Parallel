/*******************************************************
    > File   : CG.c
    > Author : Yuntong
    > Mail   : 171840067@smail.nju.edu.cn
    > Date   : 2021/4/13
*******************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "CG.h"

void CG(Stencil *A, double *diag, double *rhs, double *x,
        double EPSILON, int maxit, int *numit,
        int nlx, int nly, int myidx, int myidy, int npx, int npy) {
    // compute dest and src (in four directions)
    int left, right, above, below;
    left = myidx - 1;
    if (left >= 0) left = myidy * npx + left;
    else left = MPI_PROC_NULL;

    right = myidx + 1;
    if (right < npx) right = myidy * npx + right;
    else right = MPI_PROC_NULL;

    above = myidy + 1;
    if (above < npy) above = above * npx + myidx;
    else above = MPI_PROC_NULL;

    below = myidy - 1;
    if (below >= 0) below = below * npx + myidx;
    else below = MPI_PROC_NULL;

    MPI_Status status;
    // define three mpi datatype
    MPI_Datatype col_d, col, row;
    MPI_Type_vector(nly, 1, nlx, MPI_DOUBLE, &col_d);
    MPI_Type_contiguous(nly, MPI_DOUBLE, &col);
    MPI_Type_contiguous(nlx, MPI_DOUBLE, &row);

    MPI_Type_commit(&col_d);
    MPI_Type_commit(&col);
    MPI_Type_commit(&row);

    int i, nits;
    int mesh_size = nly * nlx;
    double alpha, beta;
    double dot_res1, dot_res2, dot_d;     //global
    double residual, residual0, descent;  //global
    double dot_res2_tmp, dot_d_tmp;       //local
    double *tmp, *d, *Ad, *res;
    double *west, *east, *north, *south;
    tmp     = (double *) calloc(mesh_size, sizeof(double));
    d     = (double *) calloc(mesh_size, sizeof(double));
    Ad    = (double *) calloc(mesh_size, sizeof(double));
    res   = (double *) calloc(mesh_size, sizeof(double));
    west  = (double *) calloc(nly,       sizeof(double));
    east  = (double *) calloc(nly,       sizeof(double));
    north = (double *) calloc(nlx,       sizeof(double));
    south = (double *) calloc(nlx,       sizeof(double));
    
    //x_0
    for (i = 0; i < mesh_size; i++)  x[i] = 0;
    // west d(0, 0:nly)      -->   left   -->  west
    MPI_Sendrecv(x,     1, col_d, left, 111,
                 west,  1, col,   left, 111,
                 MPI_COMM_WORLD, &status);

    // east d(nlx, 0:nly-1)  -->   right  -->  east
    MPI_Sendrecv(&x[nlx-1], 1, col_d, right, 111,
                 east,      1, col,   right, 111,
                 MPI_COMM_WORLD, &status);

    // north d(0:nlx, nly-1) -->   above  -->  north
    MPI_Sendrecv(&x[mesh_size - nlx], 1, row, above, 111,
                 north,               1, row, above, 111,
                 MPI_COMM_WORLD, &status);

    // south d(0:nlx, 0)     -->   below  --> south
    MPI_Sendrecv(x,      1, row, below, 111,
                 south,  1, row, below, 111,
                 MPI_COMM_WORLD, &status);

    //d_0 = r_0 = b - A * x_0
    MV_stencil(A, x, west, east, north, south, tmp, nlx, nly);
    for (i = 0; i < mesh_size; i++) {
        tmp[i] += diag[i] * x[i];
    }
    vector_subtraction(rhs, tmp, res, mesh_size);
    for (i = 0; i < mesh_size; i++)  d[i] = res[i];

    dot_res2_tmp = dot(res, res, mesh_size);
    MPI_Allreduce(&dot_res2_tmp, &residual0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if (myidy * npx + myidx == 0) {
        printf("0 : \t%.8e\n", residual0);
        printf("------------------------------------------------------\n");
    }
    dot_res2 = residual0;
    
    for (nits = 1; nits <= maxit; nits++) {
        // update inner boundary values
        // west d(0, 0:nly)      -->   left   -->  west
        MPI_Sendrecv(d,     1, col_d, left, 111,
                     west,  1, col,   left, 111,
                     MPI_COMM_WORLD, &status);

        // east d(nlx, 0:nly-1)  -->   right  -->  east
        MPI_Sendrecv(&d[nlx-1], 1, col_d, right, 111,
                     east,      1, col,   right, 111,
                     MPI_COMM_WORLD, &status);

        // north d(0:nlx, nly-1) -->   above  -->  north
        MPI_Sendrecv(&d[mesh_size - nlx], 1, row, above, 111,
                     north,               1, row, above, 111,
                     MPI_COMM_WORLD, &status);

        // south d(0:nlx, 0)     -->   below  --> south
        MPI_Sendrecv(d,      1, row, below, 111,
                     south,  1, row, below, 111,
                     MPI_COMM_WORLD, &status);

        // Ad = (A - D) * d + D * d
        MV_stencil(A, d, west, east, north, south, Ad, nlx, nly);
        for (i = 0; i < mesh_size; i++) {
            Ad[i] += diag[i] * d[i];
        }

        //inner product
        dot_res1 = dot_res2;
        dot_d_tmp = dot(d, Ad, mesh_size);
        MPI_Allreduce(&dot_d_tmp, &dot_d, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        //compute x_new
        alpha = dot_res1 / dot_d;
        const_vector_multiply(alpha, d, tmp, mesh_size);
        vector_addition(x, tmp, x, mesh_size);
        
        //update residul
        const_vector_multiply(alpha, Ad, Ad, mesh_size);
        vector_subtraction(res, Ad, res, mesh_size);
        dot_res2_tmp = dot(res, res, mesh_size);
        MPI_Allreduce(&dot_res2_tmp, &dot_res2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        //update direction
        beta = dot_res2 / dot_res1;
        const_vector_multiply(beta, d, d, mesh_size);
        vector_addition(res, d, d, mesh_size);


        // Termination criterion
        residual = dot_res2;
        descent = sqrt(residual / residual0);
        if ((myidy * npx + myidx == 0) & (nits % 20 == 0 | nits == 1)) {
            printf("%d : \t%.8e\n", nits, descent);
        }
        if (descent < EPSILON) {
            if (myidy * npx + myidx == 0) printf("%d : \t%.8e\n", nits, descent);
            break;
        }
    }
    
    if (myidy * npx + myidx == 0) printf("\n");
    *numit = nits;

    MPI_Type_free(&col_d);
    MPI_Type_free(&col);
    MPI_Type_free(&row);
    
    free(tmp);  free(d);    free(Ad);    free(res);
    free(west); free(east); free(north); free(south);
    return ;
}


void PCG(Stencil *A, double *diag, double *rhs, double *x, 
        int precond, int precond_iter, double EPSILON, int maxit, int *numit,
        int nlx, int nly, int myidx, int myidy, int npx, int npy) {
    // compute dest and src (in four directions)
    int left, right, above, below;
    left = myidx - 1;
    if (left >= 0) left = myidy * npx + left;
    else left = MPI_PROC_NULL;

    right = myidx + 1;
    if (right < npx) right = myidy * npx + right;
    else right = MPI_PROC_NULL;

    above = myidy + 1;
    if (above < npy) above = above * npx + myidx;
    else above = MPI_PROC_NULL;

    below = myidy - 1;
    if (below >= 0) below = below * npx + myidx;
    else below = MPI_PROC_NULL;

    MPI_Status status;
    // define three mpi datatype
    MPI_Datatype col_d, col, row;
    MPI_Type_vector(nly, 1, nlx, MPI_DOUBLE, &col_d);
    MPI_Type_contiguous(nly, MPI_DOUBLE, &col);
    MPI_Type_contiguous(nlx, MPI_DOUBLE, &row);

    MPI_Type_commit(&col_d);
    MPI_Type_commit(&col);
    MPI_Type_commit(&row);

    int i, nits;
    int mesh_size = nly * nlx;
    double alpha, beta;
    double dot_res1, dot_res2, dot_d;     //global
    double residual, residual0, descent;  //global
    double dot_res2_tmp, dot_d_tmp;       //local
    double *tmp, *d, *Ad, *res, *res_p;
    double *west, *east, *north, *south;
    tmp   = (double *) calloc(mesh_size, sizeof(double));
    d     = (double *) calloc(mesh_size, sizeof(double));
    Ad    = (double *) calloc(mesh_size, sizeof(double));
    res   = (double *) calloc(mesh_size, sizeof(double));
    res_p = (double *) calloc(mesh_size, sizeof(double));
    west  = (double *) calloc(nly,       sizeof(double));
    east  = (double *) calloc(nly,       sizeof(double));
    north = (double *) calloc(nlx,       sizeof(double));
    south = (double *) calloc(nlx,       sizeof(double));
    
    //Preconditioner : M^{-1}
    // z_0 = M^{-1} * rhs;
    switch (precond) {
        case 1: Precond_JAC(A, diag, rhs, res_p, EPSILON, nlx, nly, precond_iter); break; 
        case 2: Precond_GS(A, diag, rhs, res_p, EPSILON, nlx, nly, precond_iter);  break; 
        case 3: Precond_CG(A, diag, rhs, res_p, EPSILON, nlx, nly, precond_iter);  break; 
    }

    //x_0
    for (i = 0; i < mesh_size; i++)  x[i] = 0;
    // west d(0, 0:nly)      -->   left   -->  west
    MPI_Sendrecv(x,     1, col_d, left, 111,
                 west,  1, col,   left, 111,
                 MPI_COMM_WORLD, &status);

    // east d(nlx, 0:nly-1)  -->   right  -->  east
    MPI_Sendrecv(&x[nlx-1], 1, col_d, right, 111,
                 east,      1, col,   right, 111,
                 MPI_COMM_WORLD, &status);

    // north d(0:nlx, nly-1) -->   above  -->  north
    MPI_Sendrecv(&x[mesh_size - nlx], 1, row, above, 111,
                 north,               1, row, above, 111,
                 MPI_COMM_WORLD, &status);

    // south d(0:nlx, 0)     -->   below  --> south
    MPI_Sendrecv(x,      1, row, below, 111,
                 south,  1, row, below, 111,
                 MPI_COMM_WORLD, &status);

    //r_0 = b - A * x_0; d_0 = z_0
    MV_stencil(A, x, west, east, north, south, tmp, nlx, nly);
    for (i = 0; i < mesh_size; i++) {
        tmp[i] += diag[i] * x[i];
    }
    vector_subtraction(rhs, tmp, res, mesh_size);
    for (i = 0; i < mesh_size; i++)  d[i] = res_p[i];

    dot_res2_tmp = dot(res, res, mesh_size);
    MPI_Allreduce(&dot_res2_tmp, &residual0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if (myidy * npx + myidx == 0) {
        printf("0 : \t%.8e\n", residual0);
        printf("------------------------------------------------------\n");
    }
    dot_res2_tmp = dot(res, res_p, mesh_size);
    MPI_Allreduce(&dot_res2_tmp, &dot_res2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    for (nits = 1; nits <= maxit; nits++) {
        // update inner boundary values
        // west d(0, 0:nly)      -->   left   -->  west
        MPI_Sendrecv(d,     1, col_d, left, 111,
                     west,  1, col,   left, 111,
                     MPI_COMM_WORLD, &status);

        // east d(nlx, 0:nly-1)  -->   right  -->  east
        MPI_Sendrecv(&d[nlx-1], 1, col_d, right, 111,
                     east,      1, col,   right, 111,
                     MPI_COMM_WORLD, &status);

        // north d(0:nlx, nly-1) -->   above  -->  north
        MPI_Sendrecv(&d[mesh_size - nlx], 1, row, above, 111,
                     north,               1, row, above, 111,
                     MPI_COMM_WORLD, &status);

        // south d(0:nlx, 0)     -->   below  --> south
        MPI_Sendrecv(d,      1, row, below, 111,
                     south,  1, row, below, 111,
                     MPI_COMM_WORLD, &status);

        // Ad = (A - D) * d + D * d
        MV_stencil(A, d, west, east, north, south, Ad, nlx, nly);
        for (i = 0; i < mesh_size; i++) {
            Ad[i] += diag[i] * d[i];
        }

        //inner product
        dot_res1  = dot_res2;
        dot_d_tmp = dot(d, Ad, mesh_size);
        MPI_Allreduce(&dot_d_tmp, &dot_d, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

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
            case 1: Precond_JAC(A, diag, res, res_p, EPSILON, nlx, nly, precond_iter); break; 
            case 2: Precond_CG(A, diag, res, res_p, EPSILON, nlx, nly, precond_iter);  break; 
        }
        dot_res2_tmp = dot(res, res_p, mesh_size);
        MPI_Allreduce(&dot_res2_tmp, &dot_res2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        //update direction
        // d_j = z_j + beta * d_{j-1}
        beta = dot_res2 / dot_res1;
        const_vector_multiply(beta, d, d, mesh_size);
        vector_addition(res_p, d, d, mesh_size);


        // Termination criterion
        dot_res2_tmp = dot(res, res, mesh_size);
        MPI_Allreduce(&dot_res2_tmp, &residual, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        descent = sqrt(residual / residual0);
        if ((myidy * npx + myidx == 0) & (nits % 50 == 0 | nits == 1)) {
            printf("%d : \t%.8e\n", nits, descent);
        }
        if (descent < EPSILON) {
            if (myidy * npx + myidx == 0) printf("%d : \t%.8e\n", nits, descent);
            break;
        }
    }
    
    if (myidy * npx + myidx == 0) printf("\n");
    *numit = nits;

    MPI_Type_free(&col_d);
    MPI_Type_free(&col);
    MPI_Type_free(&row);
    
    free(tmp);  free(d);    free(Ad);    free(res);    free(res_p);
    free(west); free(east); free(north); free(south);
    return ;    
}



