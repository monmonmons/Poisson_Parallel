/*******************************************************
    > File   : Jacobi.c
    > Author : Yuntong
    > Mail   : 171840067@smail.nju.edu.cn
    > Date   : 2021/4/11
*******************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include "Jacobi.h"


void Jacobi(Stencil *A, double *diag, double *rhs, double *x,
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
    int mesh_size = nly * nlx;           //local
    double residual_tmp;                 //local
    double residual, residual0, descent; //global
    double *y, *c, *res;
    Stencil *D_A;
    double *west, *east, *north, *south;
    y     = (double *) calloc(mesh_size, sizeof(double));
    c     = (double *) calloc(mesh_size, sizeof(double));
    res   = (double *) calloc(mesh_size, sizeof(double));
    D_A   =(Stencil *) calloc(mesh_size, sizeof(Stencil));
    west  = (double *) calloc(nly,       sizeof(double));
    east  = (double *) calloc(nly,       sizeof(double));
    north = (double *) calloc(nlx,       sizeof(double));
    south = (double *) calloc(nlx,       sizeof(double));

    
    for (i = 0; i < mesh_size; i++) {
        D_A[i].value[0] = A[i].value[0] / diag[i];    // D^(-1) * (A - D)
        D_A[i].value[1] = A[i].value[1] / diag[i];
        D_A[i].value[2] = A[i].value[2] / diag[i];
        D_A[i].value[3] = A[i].value[3] / diag[i];
          c[i]          = rhs[i]        / diag[i];    // c = D^(-1) * b
          x[i]          = 0;
    }
    residual_tmp = dot(rhs, rhs, mesh_size);
    MPI_Allreduce(&residual_tmp, &residual0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if (myidy * npx + myidx == 0) {
        printf("0 : \t%.8e\n", residual0);
        printf("------------------------------------------------------\n");
    }
    
    for (nits = 1; nits <= maxit; nits++) {
        // y = D^(-1)*(A-D) * x_old
        MV_stencil(D_A, x, west, east, north, south, y, nlx, nly);
        // x = c - y
        vector_subtraction(c, y, x, mesh_size);

        // update inner boundary values
        // west x(0, 0:nly)      -->   left   -->  west
        MPI_Sendrecv(x,     1, col_d, left, 111,
                     west,  1, col,   left, 111,
                     MPI_COMM_WORLD, &status);

        // east x(nlx, 0:nly-1)  -->   right  -->  east
        MPI_Sendrecv(&x[nlx-1], 1, col_d, right, 111,
                     east,      1, col,   right, 111,
                     MPI_COMM_WORLD, &status);

        // north x(0:nlx, nly-1) -->   above  -->  north
        MPI_Sendrecv(&x[mesh_size - nlx], 1, row, above, 111,
                     north,               1, row, above, 111,
                     MPI_COMM_WORLD, &status);

        // south x(0:nlx, 0)     -->   below  --> south
        MPI_Sendrecv(x,      1, row, below, 111,
                     south,  1, row, below, 111,
                     MPI_COMM_WORLD, &status);

        // update residual
        MV_stencil(D_A, x, west, east, north, south, res, nlx, nly);
        vector_addition(res, x, res, mesh_size);
        for (i = 0; i < mesh_size; i++) {
            res[i] = diag[i] * res[i];
        }
        vector_subtraction(rhs, res, res, mesh_size);
        residual_tmp = dot(res, res, mesh_size);
        MPI_Allreduce(&residual_tmp, &residual, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        // Termination criterion
        descent = sqrt(residual / residual0);
        if ((myidy * npx + myidx == 0) & (nits % 500 == 0 | nits == 1)) {
            printf("%d : \t%.8e\n", nits, descent);
        }
        if (descent < EPSILON) break;
    }

    if (myidy * npx + myidx == 0) printf("\n");
    *numit = nits;

    MPI_Type_free(&col_d);  MPI_Type_free(&col);  MPI_Type_free(&row);

    
    free(y);    free(c);    free(res);   free(D_A);
    free(west); free(east); free(north); free(south);
    return ;
}
