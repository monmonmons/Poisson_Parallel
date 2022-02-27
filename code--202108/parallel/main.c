/*******************************************************
    > File   : main.c
    > Author : Yuntong
    > Mail   : 171840067@smail.nju.edu.cn
    > Date   : 2021/8/25 
*******************************************************/

// This program runs iteration method parallelly on linear system A * x = b;
// which is discretized from a 2d-poisson-like equation :
//   - div (k (div u)) + tau u = f
//   k = (k1, k2)

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <math.h>
#include "poisson.h"
#include "Jacobi.h"
#include "BiCGSTAB.h"
#include "CG.h"

const double PI = 3.1415926;
double EPSILON = 1e-10;
int    MAXIT;

int npx;  
int npy;

double f(double, double);     //Source term
double g(double, double);     //BC

double u(double, double);     //exact solution

double k1(double, double);    //kappa -- x
double k2(double, double);	  //kappa -- y

int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);

	double tau, tau0; 		// parameters of poisson equation
	double a1, a2;
	int boundary_condition; // 1--Dirichlet, 2--Newmann, 
						    // 3--mixed (west--Dirichlet, others--Newmann)

							// parameters of discretization
	int size;				// n
	double mesh_step;		// h = 1.0 / n

	int method;
	int num_iter;
	int preconditioner, precond_iter;
    int i, j, num;

	/***************************************************
    	Initialize parameters from file
	 ***************************************************/
	char s[20], ch;
	FILE *in_fp = fopen("input.txt", "r");
	if (in_fp == NULL) {
		printf("无法打开文件。\n");
		return -1;
	}
	fscanf(in_fp, "%s %c %d",  s, &ch, &MAXIT);
	fscanf(in_fp, "%s %c %lf", s, &ch, &EPSILON);
	fscanf(in_fp, "%s %c %lf", s, &ch, &tau0);
	fscanf(in_fp, "%s %c %lf", s, &ch, &a1);
	fscanf(in_fp, "%s %c %lf", s, &ch, &a2);
	fscanf(in_fp, "%s %c %d",  s, &ch, &boundary_condition);
	fscanf(in_fp, "%s %c %d",  s, &ch, &size);
	fscanf(in_fp, "%s %c %d",  s, &ch, &method);
	fscanf(in_fp, "%s %c %d",  s, &ch, &preconditioner);
	fscanf(in_fp, "%s %c %d",  s, &ch, &precond_iter);
	fscanf(in_fp, "%s %c %d",  s, &ch, &npx);
	fscanf(in_fp, "%s %c %d",  s, &ch, &npy);
	fclose(in_fp);

	/***************************************************
    	Initialize processors
	 ***************************************************/
	int np, myid;	   //number of processors
	int myidx, myidy;  //myidx = myid % npx; myidy = myid % npy
	int nlx, nly;      //mesh size of subdomain
	int mesh_size;

	MPI_Comm_size(MPI_COMM_WORLD, &np);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if (npx*npy != np) {
        if ( myid == 0)
           printf("Error: npx*npy not equal to np!\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    myidx = myid % npx;
    myidy = myid / npx;

    nlx = size / npx ;
    nly = size / npy ;
    mesh_size = nlx * nly;
    mesh_step = 1.0 / size;
    tau = tau0 * size * size;

	/***************************************************
    	Generate stencil : (A - D) & D
		Generate vector  : b
	 ***************************************************/
	Stencil *matrix = (Stencil *) calloc(mesh_size, sizeof(Stencil));
	double  *diag   = (double *)  calloc(mesh_size, sizeof(double));
	double  *b 		= (double *)  calloc(mesh_size, sizeof(double));
	double  *x 		= (double *)  calloc(mesh_size, sizeof(double));
	double  *exact 	= (double *)  calloc(mesh_size, sizeof(double));
    if (!matrix | !diag | !b | !x | !exact) {	
        if ( myid == 0)
            printf("Out of memory!\n");
        MPI_Abort(MPI_COMM_WORLD, 2);
    }

	if (boundary_condition == 1) {
		matrix_stencil_mode(matrix, diag, k1, k2, tau, size, mesh_step, nlx, nly, myidx, myidy, npx, npy, 1);

	    RHS_vector_Dirichlet(b, k1, k2, f, g, size, mesh_step, nlx, nly, myidx, myidy, npx, npy);

	} else if (boundary_condition == 2) {
		matrix_stencil_mode(matrix, diag, k1, k2, tau, size, mesh_step, nlx, nly, myidx, myidy, npx, npy, 2);

		RHS_vector_Neumann(b, f, size, mesh_step, nlx, nly, myidx, myidy);

	} else if (boundary_condition == 3) {
		matrix_stencil_mixed(matrix, diag, k1, k2, a1, a2, tau, size, mesh_step, nlx, nly, myidx, myidy, npx, npy);

	    RHS_vector_mixed(b, k1, k2, f, g, a1, a2, size, mesh_step, nlx, nly, myidx, myidy, npx, npy);

	} else {
		printf("Undefined boundary condition!\n");
     	return -1;
	}

	/***************************************************
    	Iteration methods
	 ***************************************************/
    double time;
    
    #if 1
	    if (method == 1) {
	    	//CG method
	    	time = MPI_Wtime();
		    CG(matrix, diag, b, x, EPSILON, MAXIT, &num_iter, nlx, nly, myidx, myidy, npx, npy);
		    time = MPI_Wtime() - time; 

		} else if (method == 2) {
			// Preconditioned CG method
		    time = MPI_Wtime();
		    PCG(matrix, diag, b, x, preconditioner, precond_iter, EPSILON, MAXIT, &num_iter, nlx, nly, myidx, myidy, npx, npy);
		    time = MPI_Wtime() - time;

		} else if (method == 3) {
			//BiCGSTAB method
			time = MPI_Wtime();
			BICGSTAB(matrix, diag, b, x, EPSILON, MAXIT, &num_iter, nlx, nly, myidx, myidy, npx, npy);
			time = MPI_Wtime() - time; 

		} else if (method == 4) {
			// Preconditioned BiCGSTAB method
			time = MPI_Wtime();
			P_BICGSTAB(matrix, diag, b, x, preconditioner, precond_iter, EPSILON, MAXIT, &num_iter, nlx, nly, myidx, myidy, npx, npy);
			time = MPI_Wtime() - time; 

		}  else if (method == 5) {
			//Jacobi method
			time = MPI_Wtime();
			Jacobi(matrix, diag, b, x, EPSILON, MAXIT, &num_iter, nlx, nly, myidx, myidy, npx, npy);
			time = MPI_Wtime() - time; 

		} else {
			printf("Undefined iteration method!\n");
			return -1;
		}
	#endif

	//test 10 times
	#if 0
		int count;
		if (method == 1) {
	    	//CG method
	    	for (count = 0; count < 10; count++) {
	    		time = MPI_Wtime();
			    CG(matrix, diag, b, x, EPSILON, MAXIT, &num_iter, nlx, nly, myidx, myidy, npx, npy);
			    time = MPI_Wtime() - time; 
	    	}

		} else if (method == 2) {
			// Preconditioned CG method
			for (count = 0; count < 10; count++) {
	    		time = MPI_Wtime();
			    PCG(matrix, diag, b, x, preconditioner, precond_iter, EPSILON, MAXIT, &num_iter, nlx, nly, myidx, myidy, npx, npy);
			    time = MPI_Wtime() - time;
	    	}

		} else if (method == 3) {
			//BiCGSTAB method
			for (count = 0; count < 10; count++) {
	    		time = MPI_Wtime();
				BICGSTAB(matrix, diag, b, x, EPSILON, MAXIT, &num_iter, nlx, nly, myidx, myidy, npx, npy);
				time = MPI_Wtime() - time; 
	    	}

		} else if (method == 4) {
			// Preconditioned BiCGSTAB method
			for (count = 0; count < 10; count++) {
	    		time = MPI_Wtime();
				P_BICGSTAB(matrix, diag, b, x, preconditioner, precond_iter, EPSILON, MAXIT, &num_iter, nlx, nly, myidx, myidy, npx, npy);
				time = MPI_Wtime() - time; 
	    	}

		}  else if (method == 5) {
			//Jacobi method
			for (count = 0; count < 10; count++) {
	    		time = MPI_Wtime();
				Jacobi(matrix, diag, b, x, EPSILON, MAXIT, &num_iter, nlx, nly, myidx, myidy, npx, npy);
				time = MPI_Wtime() - time; 
	    	}

		} else {
			printf("Undefined iteration method!\n");
			return -1;
		}
		time /= 10;
	#endif

	/***************************************************
    	Correctness verification
	 ***************************************************/

	double error_tmp = 0, error = 0;
	//exact solution
    int x0 = myidx * nlx;
    int y0 = myidy * nly;
    int xi, yj;
    for (j = 0; j < nly; j++) {
    	yj = y0 + j;
        for (i = 0; i < nlx; i++) {
        	xi = x0 + i;
            num = j * nlx + i;
            exact[num] = u((xi + 0.5) * mesh_step, (yj + 0.5) * mesh_step);
        }
    }
    vector_subtraction(x, exact, exact, nlx * nly);

    error_tmp = norm_inf(exact, nlx * nly);
    MPI_Allreduce(&error_tmp, &error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	/***************************************************
    	Print result
	 ***************************************************/
	
	if(myid == 0) {
		printf("------------------------------------------------------\n");
		printf("Mesh size :\t%d = %d * %d\n", size*size, size, size);
		printf("------------------------------------------------------\n");
		printf("Iteration method: ");
		char *method_s = NULL, *precond_s = NULL;
		if (method == 1) {
			method_s = "CG";
		} else if (method == 2) {
			if (preconditioner == 1) {
				precond_s = "Jacobi";
			} else if (preconditioner == 2) {
				precond_s = "GS";
			} else if (preconditioner == 3) {
				precond_s = "CG";
			}
			method_s = "CG";
		} else if (method == 3) {
			method_s = "BiCGSTAB";
		} else if (method == 4) {
			if (preconditioner == 1) {
				precond_s = "Jacobi";
			} else if (preconditioner == 2) {
				precond_s = "GS";
			} else if (preconditioner == 3) {
				precond_s = "CG";
			}
			method_s = "BiCGSTAB";
		} else if (method == 5) {
			method_s = "Jacobi";
		}
		if (precond_s) printf("%s_%s\t", precond_s, method_s);
		else printf("%s\t", method_s);
		printf("tau = %.2e\n", tau0);
		printf("------------------------------------------------------\n");
		printf("Number of iterations  : \t");
		if (precond_s) printf("%d_%d\n", precond_iter, num_iter);
		else printf("%d\n", num_iter);
		if (np != 1) {
			printf("Number of processors  : \t%d = %d * %d\n", np, npx, npy);
		}
	    printf("Time of computation   : \t%.8e s\n", time);
	    printf("Discretization error  : \t%.8e\n", error);
	    printf("\n");
	
	    #if 1
	    	in_fp = fopen("result.txt", "a");
		    if (in_fp == NULL) {
				printf("无法打开文件。\n");
				return -1;
			}
			if (precond_s) fprintf(in_fp, "%s_%s\t", precond_s, method_s);
			else fprintf(in_fp, "%s\t", method_s);
			fprintf(in_fp, "%d\t", size);
			fprintf(in_fp, "%.2e\t", tau0);
			if (precond_s) fprintf(in_fp, "%d_%d\t", precond_iter,num_iter);
			else fprintf(in_fp, "%d\t", num_iter);
			fprintf(in_fp, "%d\t", np);
			fprintf(in_fp, "%d * %d\t", npx, npy);
			fprintf(in_fp, "%.8e\t", time);
			fprintf(in_fp, "%.8e\t", error);
			fprintf(in_fp, "\n");
			fclose(in_fp);
	    #endif
	}
	/***************************************************
    	free pointer
	 ***************************************************/
	free(matrix);  free(diag);  free(b);  free(x);  free(exact);

	MPI_Finalize();
	return 0;
}

//Source term
double f(double x, double y) {
    // return (2 * PI * PI + 0) * sin(PI * x) * sin(PI * y);
    return 1.0;
}

//Dirichlet BC
double g(double x, double y) {
    // return sin(PI * x) * sin(PI * y);
    return - ((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5)) / 4.0;
    // return 0;
}

double u(double x, double y) {
    // return sin(PI * x) * sin(PI * y);
    return - ((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5)) / 4.0;
}

double k1(double x, double y) {
    return 1.0;
}

double k2(double x, double y) {
    return 1.0;
}

