/*******************************************************
    > File   : poisson.h
    > Author : Yuntong
    > Mail   : 171840067@smail.nju.edu.cn
    > Date   : 2021/8/15
*******************************************************/

#ifndef poisson_h
#define poisson_h

#include "ops.h"

//Generate stiffness matrix
//	Input    : k1, k2, tau, size, h, bc
//	Modified : matrix, diag

// 1--Dirichlet BC (u = g on gamma)
// 2--Newmann BC   (du / dn = 0)
void matrix_stencil_mode_const(Stencil *matrix, double *diag,
                         double k1, double k2, double tau, int size, double h, int bc);


void matrix_stencil_mode(Stencil *matrix, double *diag, double (*k1)(double, double), 
                         double (*k2)(double, double), double tau, int size, double h, int bc);


// mixed BC  (a1*u + a2* du/dn = g)
void matrix_stencil_mixed(Stencil *matrix, double *diag, 
                          double (*k1)(double, double), double (*k2)(double, double), 
                          double a1, double a2, double tau, int size, double h);




//Generate right-hand-side vector
//	Input    : diag, f, g, size, h, bc
//	Modified : b

// Dirichlet BC (u = g on gamma)
void RHS_vector_Dirichlet_const(double *b, double k1, double k2, 
                                double (*f)(double, double), double (*g)(double, double), 
                                int size, double h);

void RHS_vector_Dirichlet(double *b, double (*k1)(double, double), double (*k2)(double, double),
                                     double (*f)(double, double),  double (*g)(double, double), 
                                     int size, double h);

// Newmann BC   (du / dn = 0)
void RHS_vector_Neumann(double *b, double (*f)(double, double), int size, double h);



// mixed BC  (a1*u + a2* du/dn = g)
void RHS_vector_mixed(double *b, double (*k1)(double, double), double (*k2)(double, double),
                      double (*f)(double, double), double (*g)(double, double), 
                      double a1, double a2, int size, double h);

#endif /* poisson_h */
