/*******************************************************
    > File   : main.c
    > Author : Yuntong
    > Mail   : 171840067@smail.nju.edu.cn
    > Date   : 2021/3/18
*******************************************************/

// This program runs iteration method on linear system A * x = b;
// which is discretized from a 2d-poisson-like equation :
//   - div (k (div u)) + tau u = f
//   k = (k1, k2)

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "poisson.h"
#include "Jacobi.h"
#include "CG.h"
#include "BiCGSTAB.h"

const double PI = 3.1415926;
double EPSILON;
int    MAXIT;

double f(double, double);     //Source term
double g(double, double);     //BC

double u(double, double);     //exact solution

double k1(double, double);    //kappa -- x
double k2(double, double);	  //kappa -- y

int main(){
	double tau, tau0; 	     // parameters of poisson equation
	double a1, a2;
	int boundary_condition;  // 1--Dirichlet, 2--Newmann, 
						     // 3--mixed (a1 * u + a2 * du/dn = g)
	
							 // parameters of discretization
	int size;				 // n
	int mesh_size;			 // N = n * n
	double mesh_step;		 // h = 1.0 / n

	int method;
	int num_iter = 0;
	int preconditioner, precond_iter;
    int i, j, num;
    
    clock_t start, end;

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
	fscanf(in_fp, "%s %c %lf",  s, &ch, &EPSILON);
	fscanf(in_fp, "%s %c %lf", s, &ch, &tau0);
	fscanf(in_fp, "%s %c %lf", s, &ch, &a1);
	fscanf(in_fp, "%s %c %lf", s, &ch, &a2);
	fscanf(in_fp, "%s %c %d",  s, &ch, &boundary_condition);
	fscanf(in_fp, "%s %c %d",  s, &ch, &size);
	fscanf(in_fp, "%s %c %d",  s, &ch, &method);
	fscanf(in_fp, "%s %c %d",  s, &ch, &preconditioner);
	fscanf(in_fp, "%s %c %d",  s, &ch, &precond_iter);
	fclose(in_fp);


	/***************************************************
    	Generate stencil : (A - D) & D
		Generate vector  : b
	 ***************************************************/
	mesh_size = size * size;   tau = tau0 * mesh_size;
	mesh_step = 1.0 / size;

	Stencil *matrix = (Stencil *) calloc(mesh_size, sizeof(Stencil));
	double  *diag   = (double *)  calloc(mesh_size, sizeof(double));
	double  *b 		= (double *)  calloc(mesh_size, sizeof(double));
	double  *x 		= (double *)  calloc(mesh_size, sizeof(double));
	double  *exact 	= (double *)  calloc(mesh_size, sizeof(double));
	if (!matrix | !diag | !b | !x | !exact) {	
        printf("Out of memory!\n");
    }

	if (boundary_condition == 1) {
		matrix_stencil_mode(matrix, diag, k1, k2, tau, size, mesh_step, 1);

	    RHS_vector_Dirichlet(b, k1, k2, f, g, size, mesh_step);

	} else if (boundary_condition == 2) {
		matrix_stencil_mode(matrix, diag, k1, k2, tau, size, mesh_step, 2);

		RHS_vector_Neumann(b, f, size, mesh_step);

	} else if (boundary_condition == 3) {
		matrix_stencil_mixed(matrix, diag, k1, k2, a1, a2, tau, size, mesh_step);

	    RHS_vector_mixed(b, k1, k2, f, g, a1, a2, size, mesh_step);

	} else {
		printf("Undefined boundary condition!\n");
        return -1;
	}


	// for (i = 0; i < mesh_size; i++) {
	// 	printf("num : %d, b = %.4lf \n", matrix[i].number, b[i]);
	// 	// printf("b = %.4lf ", b[i]);
	// 	// printf("diag = %.4lf  ", diag[i]);
	// 	// print_line(&matrix[i]);
	// }

	/***************************************************
    	Iteration methods
	 ***************************************************/
	double time = 0;

	for (j = 0; j < size; j++) {
        for (i = 0; i < size; i++) {
            num = j * size + i;
            exact[num] = u((i + 0.5) * mesh_step, (j + 0.5) * mesh_step);
        }
    }

	#if 1
		if (method == 1) {
			//CG method
			start = clock();
			CG(matrix, diag, b, x, EPSILON, size, MAXIT, &num_iter);
			end = clock();

		} else if (method == 2) {
			//Preconditioned CG method
		    start = clock();
		    PCG(matrix, diag, b, x, preconditioner, precond_iter, EPSILON, size, MAXIT, &num_iter);
		    end = clock();

		} else if (method == 3) {
			//BICGSTAB method
		    start = clock();
		    BICGSTAB(matrix, diag, b, x, EPSILON, size, MAXIT, &num_iter);
		    end = clock();

		} else if (method == 4) {
			//Preconditioned BICGSTAB method
		    start = clock();
		    P_BICGSTAB(matrix, diag, b, x, preconditioner, precond_iter, EPSILON, size, MAXIT, &num_iter);
		    end = clock();

		} else if (method == 5) {
			//Jacobi method
		    start = clock();
		    Jacobi(matrix, diag, b, x, EPSILON, size, MAXIT, &num_iter);
		    end = clock();

		}else {
			printf("Undefined iteration method!\n");
			return -1;
		}
		time = (double) (end - start);
	#endif

	//test 10 times
	#if 0
		int count;
		if (method == 1) {
			for (count = 0; count < 10; count++) {
				start = clock();
				CG(matrix, diag, b, x, EPSILON, size, MAXIT, &num_iter);
				end = clock();
				time += (double) (end - start);
			}

		} else if (method == 2) {
			for (count = 0; count < 10; count++) {
				start = clock();
				PCG(matrix, diag, b, x, preconditioner, precond_iter, EPSILON, size, MAXIT, &num_iter);
				end = clock();
				time += (double) (end - start);
			}

		} else if (method == 3) {
			for (count = 0; count < 10; count++) {
				start = clock();
				BICGSTAB(matrix, diag, b, x, EPSILON, size, MAXIT, &num_iter);
				end = clock();
				time += (double) (end - start);
			}

		} else if (method == 4) {
		    for (count = 0; count < 10; count++) {
				start = clock();
				P_BICGSTAB(matrix, diag, b, x, preconditioner, precond_iter, EPSILON, size, MAXIT, &num_iter);
				end = clock();
				time += (double) (end - start);
			}

		} else if (method == 5) {
			for (count = 0; count < 10; count++) {
				start = clock();
				Jacobi(matrix, diag, b, x, EPSILON, size, MAXIT, &num_iter);
				end = clock();
				time += (double) (end - start);
			}
		}else {
			printf("Undefined iteration method!\n");
			return -1;
		}

		time /= 10;
	#endif

	/***************************************************
    	Correctness verification
	 ***************************************************/

    double error = 0;
    //exact solution
    // for (j = 0; j < size; j++) {
    //     for (i = 0; i < size; i++) {
    //         num = j * size + i;
    //         exact[num] = u((i + 0.5) * mesh_step, (j + 0.5) * mesh_step);
    //     }
    // }

    vector_subtraction(x, exact, exact, mesh_size);
    error = norm_inf(exact, mesh_size);


	/***************************************************
    	Print result
	 ***************************************************/

    char *method_s = NULL, *precond_s = NULL;
    // double time = (double)(end - start) / CLOCKS_PER_SEC;
    time /= CLOCKS_PER_SEC;
    double flops = mesh_size * num_iter / time;
    printf("------------------------------------------------------\n");
	printf("Mesh size :\t%d = %d * %d\n", mesh_size, size, size);
	printf("------------------------------------------------------\n");
	printf("Iteration method: ");
	if (method == 1) {
		method_s = "CG";
		flops *= 19;
	} else if (method == 2) {
		if (preconditioner == 1) {
			precond_s = "Jacobi";
		} else if (preconditioner == 2) {
			precond_s = "CG";
		}
		method_s = "CG";
	} else if (method == 3) {
		method_s = "BiCGSTAB";
		flops *= 34;
	} else if (method == 4) {
		if (preconditioner == 1) {
			precond_s = "Jacobi";
		} else if (preconditioner == 2) {
			precond_s = "CG";
		}
		method_s = "BiCGSTAB";
	} else if (method == 5) {
		method_s = "Jacobi";
		flops *= 19;
	}
	if (precond_s) printf("%s_%s\t", precond_s, method_s);
	else printf("%s\t", method_s);
	printf("tau = %.2e\n", tau0);
	printf("------------------------------------------------------\n");
	printf("Number of iterations  : \t");
	if (precond_s) printf("%d_%d\n", precond_iter, num_iter);
	else printf("%d\n", num_iter);
    printf("Time of computation   :\t %.8e s\n", time);
    // printf("Number of flop per sec:\t %.8e \n", flops);
    printf("Discretization error  :\t %.8e\n", error);
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
		fprintf(in_fp, "%.6e\t",  time);
		fprintf(in_fp, "%.8e\t", error);
		fprintf(in_fp, "\n");
		fclose(in_fp);
    #endif

	/***************************************************
    	free pointer
	 ***************************************************/
	free(matrix);  free(diag);  free(b);  free(x);  free(exact);
	return 0;
}

//Source term
double f(double x, double y) {
    // return (2 * PI * PI) * sin(PI * x) * sin(PI * y);
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
