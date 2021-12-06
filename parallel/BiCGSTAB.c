/*******************************************************
    > File   : BiCGSTAB.c
    > Author : Yuntong
    > Mail   : 171840067@smail.nju.edu.cn
    > Date   : 2021/4/21
*******************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "BiCGSTAB.h"

void BICGSTAB(Stencil *A, double *diag, double *rhs, double *x,
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
    double alpha, beta, omega;
    double dot_res1, dot_res2, dot_d;	//global
    double dot_s1, dot_s2;				//global
    double dot_tmp;						//local
    double residual, residual0, descent;
    double *tmp, *d, *s, *Ad, *As, *res;
    double *west, *east, *north, *south;
    tmp  = (double *) calloc(mesh_size, sizeof(double));
    d    = (double *) calloc(mesh_size, sizeof(double));
    s    = (double *) calloc(mesh_size, sizeof(double));
    Ad   = (double *) calloc(mesh_size, sizeof(double));
    As   = (double *) calloc(mesh_size, sizeof(double));
    res  = (double *) calloc(mesh_size, sizeof(double));
    west = (double *) calloc(nly,       sizeof(double));
    east = (double *) calloc(nly,       sizeof(double));
    north= (double *) calloc(nlx,       sizeof(double));
    south= (double *) calloc(nlx,       sizeof(double));


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

    dot_tmp = dot(res, res, mesh_size);
    MPI_Allreduce(&dot_tmp, &residual0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
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
        
        //update s
        dot_res1 = dot_res2;
        dot_tmp = dot(Ad, rhs, mesh_size);
        MPI_Allreduce(&dot_tmp, &dot_d, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        alpha = dot_res1 / dot_d;
        const_vector_multiply(alpha, Ad, tmp, mesh_size);
        vector_subtraction(res, tmp, s, mesh_size);

        // Termination criterion
        dot_tmp = dot(s, s, mesh_size);
        MPI_Allreduce(&dot_tmp, &residual, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        descent  = sqrt(residual / residual0);
        if ((myidy * npx + myidx == 0) & (nits % 50 == 0 | nits == 1)) {
        	printf("%d : \t%.8e\n", nits, descent);
        }
        if (descent < EPSILON) {
            if (myidy * npx + myidx == 0) printf("%d : \t%.8e\n", nits, descent);
            break;
        }

        // west s(0, 0:nly)      -->   left   -->  west
        MPI_Sendrecv(s,     1, col_d, left, 111,
                     west,  1, col,   left, 111,
                     MPI_COMM_WORLD, &status);

        // east s(nlx, 0:nly-1)  -->   right  -->  east
        MPI_Sendrecv(&s[nlx-1], 1, col_d, right, 111,
                     east,      1, col,   right, 111,
                     MPI_COMM_WORLD, &status);

        // north s(0:nlx, nly-1) -->   above  -->  north
        MPI_Sendrecv(&s[mesh_size - nlx], 1, row, above, 111,
                     north,               1, row, above, 111,
                     MPI_COMM_WORLD, &status);

        // south s(0:nlx, 0)     -->   below  --> south
        MPI_Sendrecv(s,      1, row, below, 111,
                     south,  1, row, below, 111,
                     MPI_COMM_WORLD, &status);

        // As = (A - D) * s + D * s
        MV_stencil(A, s, west, east, north, south, As, nlx, nly);
        for (i = 0; i < mesh_size; i++) {
            As[i] += diag[i] * s[i];
        }

        //update x
        dot_tmp = dot(s,  As,  mesh_size);
        MPI_Allreduce(&dot_tmp, &dot_s1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        dot_tmp = dot(As, As,  mesh_size);
        MPI_Allreduce(&dot_tmp, &dot_s2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        omega  = dot_s1  / dot_s2;
        const_vector_multiply(alpha, d, tmp, mesh_size);
        vector_addition(x, tmp, x, mesh_size);
        const_vector_multiply(omega, s, tmp, mesh_size);
        vector_addition(x, tmp, x, mesh_size);
        
        //update residul
        const_vector_multiply(omega, As, tmp, mesh_size);
        vector_subtraction(s, tmp, res, mesh_size);
        dot_tmp = dot(res, rhs, mesh_size);
        MPI_Allreduce(&dot_tmp, &dot_res2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        //update direction
        beta = (dot_res2 / dot_res1) * (alpha / omega);
        const_vector_multiply(omega, Ad, tmp, mesh_size);
        vector_subtraction(d, tmp, tmp, mesh_size);
        const_vector_multiply(beta, tmp, tmp, mesh_size);
        vector_addition(res, tmp, d, mesh_size);
        
    }

    if (myidy * npx + myidx == 0) printf("\n");
    *numit = nits;

    MPI_Type_free(&col_d);
    MPI_Type_free(&col);
    MPI_Type_free(&row);
    
    free(tmp);  free(d);   free(s);   
    free(res);  free(Ad);  free(As);
    free(west); free(east); free(north); free(south);

    return ;
}

void P_BICGSTAB(Stencil *A, double *diag, double *rhs, double *x, 
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
    double alpha, beta, omega;
    double dot_res1, dot_res2, dot_d;	//global
    double dot_s1, dot_s2;				//global
    double dot_tmp;						//local
    double residual, residual0, descent;
    double *tmp, *d, *s, *y, *z, *Ay, *Az, *res;
    double *west, *east, *north, *south;
    tmp  = (double *) calloc(mesh_size, sizeof(double));
    d    = (double *) calloc(mesh_size, sizeof(double));
    s    = (double *) calloc(mesh_size, sizeof(double));
    y    = (double *) calloc(mesh_size, sizeof(double));
    z    = (double *) calloc(mesh_size, sizeof(double));
    Ay   = (double *) calloc(mesh_size, sizeof(double));
    Az   = (double *) calloc(mesh_size, sizeof(double));
    res  = (double *) calloc(mesh_size, sizeof(double));
    west = (double *) calloc(nly,       sizeof(double));
    east = (double *) calloc(nly,       sizeof(double));
    north= (double *) calloc(nlx,       sizeof(double));
    south= (double *) calloc(nlx,       sizeof(double));


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
        
    dot_tmp = dot(res, res, mesh_size);
    MPI_Allreduce(&dot_tmp, &residual0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if (myidy * npx + myidx == 0) {
        printf("0 : \t%.8e\n", residual0);
        printf("------------------------------------------------------\n");
    }
    dot_res2 = residual0;

    for (nits = 1; nits <= maxit; nits++) {
    	// y = M^{-1}*d
	    switch (precond) {
	        case 1: Precond_JAC(A, diag, d, y, EPSILON, nlx, nly, precond_iter); break;
	        case 2: Precond_CG(A, diag, d, y, EPSILON, nlx, nly, precond_iter);  break;
	    }
        // update inner boundary values
        // west y(0, 0:nly)      -->   left   -->  west
        MPI_Sendrecv(y,     1, col_d, left, 111,
                     west,  1, col,   left, 111,
                     MPI_COMM_WORLD, &status);

        // east y(nlx, 0:nly-1)  -->   right  -->  east
        MPI_Sendrecv(&y[nlx-1], 1, col_d, right, 111,
                     east,      1, col,   right, 111,
                     MPI_COMM_WORLD, &status);

        // north y(0:nlx, nly-1) -->   above  -->  north
        MPI_Sendrecv(&y[mesh_size - nlx], 1, row, above, 111,
                     north,               1, row, above, 111,
                     MPI_COMM_WORLD, &status);

        // south y(0:nlx, 0)     -->   below  --> south
        MPI_Sendrecv(y,      1, row, below, 111,
                     south,  1, row, below, 111,
                     MPI_COMM_WORLD, &status);

        // Ay = (A - D) * y + D * y
        MV_stencil(A, y, west, east, north, south, Ay, nlx, nly);
        for (i = 0; i < mesh_size; i++) {
            Ay[i] += diag[i] * y[i];
        }
        
        //update s
        dot_res1 = dot_res2;
        dot_tmp = dot(Ay, rhs, mesh_size);
        MPI_Allreduce(&dot_tmp, &dot_d, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        alpha = dot_res1 / dot_d;
        const_vector_multiply(alpha, Ay, tmp, mesh_size);
        vector_subtraction(res, tmp, s, mesh_size);

        // Termination criterion
        dot_tmp = dot(s, s, mesh_size);
        MPI_Allreduce(&dot_tmp, &residual, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        descent  = sqrt(residual / residual0);
        if ((myidy * npx + myidx == 0) & (nits % 50 == 0 | nits == 1)) {
        	printf("%d : \t%.8e\n", nits, descent);
        }
        if (descent < EPSILON) {
            if (myidy * npx + myidx == 0) printf("%d : \t%.8e\n", nits, descent);
            break;
        }

        // z = M^{-1}*s
	    switch (precond) {
	        case 1: Precond_JAC(A, diag, s, z, EPSILON, nlx, nly, precond_iter); break;
	        case 2: Precond_CG(A, diag, s, z, EPSILON, nlx, nly, precond_iter);  break;
	    }
        // west z(0, 0:nly)      -->   left   -->  west
        MPI_Sendrecv(z,     1, col_d, left, 111,
                     west,  1, col,   left, 111,
                     MPI_COMM_WORLD, &status);

        // east z(nlx, 0:nly-1)  -->   right  -->  east
        MPI_Sendrecv(&z[nlx-1], 1, col_d, right, 111,
                     east,      1, col,   right, 111,
                     MPI_COMM_WORLD, &status);

        // north z(0:nlx, nly-1) -->   above  -->  north
        MPI_Sendrecv(&z[mesh_size - nlx], 1, row, above, 111,
                     north,               1, row, above, 111,
                     MPI_COMM_WORLD, &status);

        // south z(0:nlx, 0)     -->   below  --> south
        MPI_Sendrecv(z,      1, row, below, 111,
                     south,  1, row, below, 111,
                     MPI_COMM_WORLD, &status);

        // As = (A - D) * s + D * s
        MV_stencil(A, z, west, east, north, south, Az, nlx, nly);
        for (i = 0; i < mesh_size; i++) {
            Az[i] += diag[i] * z[i];
        }

        //update x
        dot_tmp = dot(s,  Az,  mesh_size);
        MPI_Allreduce(&dot_tmp, &dot_s1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        dot_tmp = dot(Az, Az,  mesh_size);
        MPI_Allreduce(&dot_tmp, &dot_s2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        omega  = dot_s1  / dot_s2;
        const_vector_multiply(alpha, y, tmp, mesh_size);
        vector_addition(x, tmp, x, mesh_size);
        const_vector_multiply(omega, z, tmp, mesh_size);
        vector_addition(x, tmp, x, mesh_size);
        
        //update residul
        const_vector_multiply(omega, Az, tmp, mesh_size);
        vector_subtraction(s, tmp, res, mesh_size);
        dot_tmp = dot(res, rhs, mesh_size);
        MPI_Allreduce(&dot_tmp, &dot_res2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        //update direction
        beta = (dot_res2 / dot_res1) * (alpha / omega);
        const_vector_multiply(omega, Ay, tmp, mesh_size);
        vector_subtraction(d, tmp, tmp, mesh_size);
        const_vector_multiply(beta, tmp, tmp, mesh_size);
        vector_addition(res, tmp, d, mesh_size);
        
    }

    if (myidy * npx + myidx == 0) printf("\n");
    *numit = nits;

    MPI_Type_free(&col_d);
    MPI_Type_free(&col);
    MPI_Type_free(&row);
    
    free(d);    free(s);    free(y);     free(z); 
    free(tmp);  free(res);  free(Ay);    free(Az);
    free(west); free(east); free(north); free(south);

    return ;
}

