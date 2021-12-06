/*******************************************************
    > File   : poisson.c
    > Author : Yuntong
    > Mail   : 171840067@smail.nju.edu.cn
    > Date   : 2021/4/8
*******************************************************/

#include <stdio.h>
#include "poisson.h"

void matrix_stencil_mode(Stencil *matrix, double *diag, double (*k1)(double, double), 
                      double (*k2)(double, double), double tau, int size, double h,
                      int nlx, int nly, int myidx, int myidy, int npx, int npy, int bc) {
    double center = h * h * tau;
    double k1_w, k1_e, k2_n, k2_s;
    double d = 0;			//diag
    int i, j, num_l;		//local
    int xi, yj, num;		//global

    int x0 = myidx * nlx;
    int y0 = myidy * nly;

    //inner mesh
    for (j = 0; j < nly; j++) {
    	yj = y0 + j;
        for (i = 0; i < nlx; i++) {
        	xi = x0 + i;
            k1_w = average(k1((xi + 0.5) * h, (yj + 0.5) * h), k1((xi - 0.5) * h, (yj + 0.5) * h));
            k1_e = average(k1((xi + 0.5) * h, (yj + 0.5) * h), k1((xi + 1.5) * h, (yj + 0.5) * h));
            k2_n = average(k2((xi + 0.5) * h, (yj + 0.5) * h), k2((xi + 0.5) * h, (yj + 1.5) * h));
            k2_s = average(k2((xi + 0.5) * h, (yj + 0.5) * h), k2((xi + 0.5) * h, (yj - 0.5) * h));
            d = center + k1_w + k1_e + k2_n + k2_s;
            num = yj * size + xi;
            num_l = j * nlx + i;
            diag[num_l] = d;
            initial_stencil(&matrix[num_l], num, -k1_w, -k1_e, -k2_n, -k2_s);
        }
    }

    //side
    if (myidx == 0) {
	    // west
    	for (j = 0; j < nly; j++) {
    		yj = y0 + j;
	        k1_w = 0.5 * (3 * k1(0.5 * h, (yj + 0.5) * h) - k1(1.5 * h, (yj + 0.5) * h));
	        k1_e = average(k1(0.5 * h, (yj + 0.5) * h), k1(1.5 * h, (yj + 0.5) * h));
	        k2_n = average(k2(0.5 * h, (yj + 0.5) * h), k2(0.5 * h, (yj + 1.5) * h));
	        k2_s = average(k2(0.5 * h, (yj + 0.5) * h), k2(0.5 * h, (yj - 0.5) * h));
	        switch (bc) {
	            case 1: d = center + 2 * k1_w + k1_e + k2_n + k2_s; break;
	            case 2: d = center + k1_e + k2_n + k2_s; break;
	            default: break;
	        }
	        num = size * yj;
	        num_l = nlx * j;
	        diag[num_l] = d;
	        initial_stencil(&matrix[num_l], num, 0, -k1_e, -k2_n, -k2_s);
    	}
    }
    if (myidx == npx - 1) {
    	//east
    	for (j = 0; j < nly; j++) {
    		yj = y0 + j;
	        k1_w = average(k1(1 - 0.5 * h, (yj + 0.5) * h), k1(1 - 1.5 * h, (yj + 0.5) * h));
	        k1_e = 0.5 * (3 * k1(1 - 0.5 * h, (yj + 0.5) * h) - k1(1 - 1.5 * h, (yj + 0.5) * h));
	        k2_n = average(k2(1 - 0.5 * h, (yj + 0.5) * h), k2(1 - 0.5 * h, (yj + 1.5) * h));
	        k2_s = average(k2(1 - 0.5 * h, (yj + 0.5) * h), k2(1 - 0.5 * h, (yj - 0.5) * h));
	        switch (bc) {
	            case 1: d = center + 2 * k1_w + k1_e + k2_n + k2_s; break;
	            case 2: d = center + k1_e + k2_n + k2_s; break;
	            default: break;
	        }
	        num = size * yj + size - 1;
	        num_l = nlx * j + nlx - 1;
	        diag[num_l] = d;
	        initial_stencil(&matrix[num_l], num, -k1_w, 0, -k2_n, -k2_s);
    	}
    }
    if (myidy == 0) {
    	//south
    	for (i = 0; i < nlx; i++) {
    		xi = x0 + i;
	        k1_w = average(k1((xi + 0.5) * h, 0.5 * h), k1((xi - 0.5) * h, 0.5 * h));
	        k1_e = average(k1((xi + 0.5) * h, 0.5 * h), k1((xi + 1.5) * h, 0.5 * h));
	        k2_n = average(k2((xi + 0.5) * h, 0.5 * h), k2((xi + 0.5) * h, 1.5 * h));
	        k2_s = 0.5 * (3 * k2((xi + 0.5) * h, 0.5 * h) - k2((xi + 0.5) * h, 1.5 * h));
	        switch (bc) {
	            case 1: d = center + k1_w + k1_e + k2_n + 2 * k2_s; break;
	            case 2: d = center + k1_w + k1_e + k2_n; break;
	            default: break;
	        }
	        num = xi;
	        num_l = i;
	        diag[num_l] = d;
	        initial_stencil(&matrix[num_l], num, -k1_w, -k1_e, -k2_n, 0);
    	}
    }
    if (myidy == npy - 1) {
    	//north
    	for (i = 0; i < nlx; i++) {
    		xi = x0 + i;
	        k1_w = average(k1((xi + 0.5) * h, 1 - 0.5 * h), k1((xi - 0.5) * h, 1 - 0.5 * h));
	        k1_e = average(k1((xi + 0.5) * h, 1 - 0.5 * h), k1((xi + 1.5) * h, 1 - 0.5 * h));
	        k2_n = 0.5 * (3 * k2((xi + 0.5) * h, 1 - 0.5 * h) - k2((xi + 0.5) * h, 1 - 1.5 * h));
	        k2_s = average(k2((xi + 0.5) * h, 1 - 0.5 * h), k2((xi + 0.5) * h, 1 - 1.5 * h));
	        switch (bc) {
	            case 1: d = center + k1_w + k1_e + k2_n + 2 * k2_s; break;
	            case 2: d = center + k1_w + k1_e + k2_n; break;
	            default: break;
	        }
	        num = size * size - size + xi;
	        num_l = nlx * nly - nlx + i;
	        diag[num_l] = d;
	        initial_stencil(&matrix[num_l], num, -k1_w, -k1_e, 0, -k2_s);
    	}
    }

    //vertex
    if (myidx == 0 && myidy == 0) {
    	//W-S
	    k1_w = 0.5 * (3 * k1(0.5 * h, 0.5 * h) - k1(1.5 * h, 0.5 * h));
	    k1_e = average(k1(0.5 * h, 0.5 * h), k1(1.5 * h, 0.5 * h));
	    k2_n = average(k2(0.5 * h, 0.5 * h), k2(0.5 * h, 1.5 * h));
	    k2_s = 0.5 * (3 * k2(0.5 * h, 0.5 * h) - k2(0.5 * h, 1.5 * h));
	    switch (bc) {
	        case 1: d = center + 2 * k1_w + k1_e + k2_n + 2 * k2_s; break;
	        case 2: d = center + k1_e + k2_n; break;
	        default: break;
	    }
	    num = 0; 
	    num_l = 0;
	    diag[num_l] = d;
	    initial_stencil(&matrix[num_l], num, 0, -k1_e, -k2_n, 0);    	
    }
    if (myidx == 0 && myidy == npy - 1) {
    	//W-N
	    k1_w = 0.5 * (3 * k1(0.5 * h, 1 - 0.5 * h) - k1(1.5 * h, 1 - 0.5 * h));
		k1_e = average(k1(0.5 * h, 1 - 0.5 * h), k1(1.5 * h, 1 - 0.5 * h));
		k2_n = 0.5 * (3 * k2(0.5 * h, 1 - 0.5 * h) - k2(0.5 * h, 1 - 1.5 * h));
		k2_s = average(k2(0.5 * h, 1 - 0.5 * h), k2(0.5 * h, 1 - 1.5 * h));
		switch (bc) {
		    case 1: d = center + 2 * k1_w + k1_e + 2 * k2_n + k2_s; break;
		    case 2: d = center + k1_e + k2_s; break;
		    default: break;
		}
		num = size * size - size; 
		num_l = nlx * nly - nlx;
		diag[num_l] = d;
		initial_stencil(&matrix[num_l], num, 0, -k1_e, 0, -k2_s);
    }
    if (myidx == npx - 1 && myidy == 0) {
    	//E-S
	    k1_w = average(k1(1 - 0.5 * h, 0.5 * h), k1(1 - 1.5 * h, 0.5 * h));
	    k1_e = 0.5 * (3 * k1(1 - 0.5 * h, 0.5 * h) - k1(1 - 1.5 * h, 0.5 * h));
	    k2_n = average(k2(1 - 0.5 * h, 0.5 * h), k2(1 - 0.5 * h, 1.5 * h));
	    k2_s = 0.5 * (3 * k2(1 - 0.5 * h, 0.5 * h) - k2(1 - 0.5 * h, 1.5 * h));
	    switch (bc) {
            case 1: d = center + k1_w + 2 * k1_e + k2_n + 2 * k2_s; break;
	        case 2: d = center + k1_w + k2_n; break;
	        default: break;
	    }
	    num = size -  1; 
	    num_l = nlx - 1;
	    diag[num_l] = d;
	    initial_stencil(&matrix[num_l], num, -k1_w, 0, -k2_n, 0);       	
    }
    if (myidx == npx - 1 && myidy == npy - 1) {
    	//E-N
 	    k1_w = average(k1(1 - 0.5 * h, 1 - 0.5 * h), k1(1 - 1.5 * h, 1 - 0.5 * h));
	    k1_e = 0.5 * (3 * k1(1 - 0.5 * h, 1 - 0.5 * h) - k1(1 - 1.5 * h, 1 - 0.5 * h));
	    k2_n = 0.5 * (3 * k2(1 - 0.5 * h, 1 - 0.5 * h) - k2(1 - 0.5 * h, 1 - 1.5 * h));
	    k2_s = average(k2(1 - 0.5 * h, 1 - 0.5 * h), k2(1 - 0.5 * h, 1 - 1.5 * h));
	    switch (bc) {
	        case 1: d = center + k1_w + 2 * k1_e + 2 * k2_n + k2_s; break;
	        case 2: d = center + k1_w + k2_s; break;
	        default: break;
	    }
	    num = size * size -  1; 
	    num_l = nlx * nly - 1;
	    diag[num_l] = d;
	    initial_stencil(&matrix[num_l], num, -k1_w, 0, 0, -k2_s);       	
    }

    return ;
}

void matrix_stencil_mixed(Stencil *matrix, double *diag, 
                      double (*k1)(double, double), double (*k2)(double, double), 
                      double a1, double a2, double tau, int size, double h,
                      int nlx, int nly, int myidx, int myidy, int npx, int npy) {
	double center = h * h * tau;
    double k1_w, k1_e, k2_n, k2_s;
    double d = 0;			//diag
    double coef = (a1 - 2 * a2 / h) / (a1 + 2 * a2 / h);
    int i, j, num_l;		//local
    int xi, yj, num;		//global

    int x0 = myidx * nlx;
    int y0 = myidy * nly;

    //inner mesh
    for (j = 0; j < nly; j++) {
    	yj = y0 + j;
        for (i = 0; i < nlx; i++) {
        	xi = x0 + i;
            k1_w = average(k1((xi + 0.5) * h, (yj + 0.5) * h), k1((xi - 0.5) * h, (yj + 0.5) * h));
            k1_e = average(k1((xi + 0.5) * h, (yj + 0.5) * h), k1((xi + 1.5) * h, (yj + 0.5) * h));
            k2_n = average(k2((xi + 0.5) * h, (yj + 0.5) * h), k2((xi + 0.5) * h, (yj + 1.5) * h));
            k2_s = average(k2((xi + 0.5) * h, (yj + 0.5) * h), k2((xi + 0.5) * h, (yj - 0.5) * h));
            d = center + k1_w + k1_e + k2_n + k2_s;
            num = yj * size + xi;
            num_l = j * nlx + i;
            diag[num_l] = d;
            initial_stencil(&matrix[num_l], num, -k1_w, -k1_e, -k2_n, -k2_s);
        }
    }

    //side
    if (myidx == 0) {
	    // west
    	for (j = 0; j < nly; j++) {
    		yj = y0 + j;
	        k1_w = 0.5 * (3 * k1(0.5 * h, (yj + 0.5) * h) - k1(1.5 * h, (yj + 0.5) * h));
	        k1_e = average(k1(0.5 * h, (yj + 0.5) * h), k1(1.5 * h, (yj + 0.5) * h));
	        k2_n = average(k2(0.5 * h, (yj + 0.5) * h), k2(0.5 * h, (yj + 1.5) * h));
	        k2_s = average(k2(0.5 * h, (yj + 0.5) * h), k2(0.5 * h, (yj - 0.5) * h));
	        d = center + (coef + 1) * k1_w + k1_e + k2_n + k2_s;
	        num = size * yj;
	        num_l = nlx * j;
	        diag[num_l] = d;
	        initial_stencil(&matrix[num_l], num, 0, -k1_e, -k2_n, -k2_s);
    	}
    }
    if (myidx == npx - 1) {
    	//east
    	for (j = 0; j < nly; j++) {
    		yj = y0 + j;
	        k1_w = average(k1(1 - 0.5 * h, (yj + 0.5) * h), k1(1 - 1.5 * h, (yj + 0.5) * h));
	        k1_e = 0.5 * (3 * k1(1 - 0.5 * h, (yj + 0.5) * h) - k1(1 - 1.5 * h, (yj + 0.5) * h));
	        k2_n = average(k2(1 - 0.5 * h, (yj + 0.5) * h), k2(1 - 0.5 * h, (yj + 1.5) * h));
	        k2_s = average(k2(1 - 0.5 * h, (yj + 0.5) * h), k2(1 - 0.5 * h, (yj - 0.5) * h));
	        d = center + k1_w + (coef + 1) * k1_e + k2_n + k2_s;
	        num = size * yj + size - 1;
	        num_l = nlx * j + nlx - 1;
	        diag[num_l] = d;
	        initial_stencil(&matrix[num_l], num, -k1_w, 0, -k2_n, -k2_s);
    	}
    }
    if (myidy == 0) {
    	//south
    	for (i = 0; i < nlx; i++) {
    		xi = x0 + i;
	        k1_w = average(k1((xi + 0.5) * h, 0.5 * h), k1((xi - 0.5) * h, 0.5 * h));
	        k1_e = average(k1((xi + 0.5) * h, 0.5 * h), k1((xi + 1.5) * h, 0.5 * h));
	        k2_n = average(k2((xi + 0.5) * h, 0.5 * h), k2((xi + 0.5) * h, 1.5 * h));
	        k2_s = 0.5 * (3 * k2((xi + 0.5) * h, 0.5 * h) - k2((xi + 0.5) * h, 1.5 * h));
	        d = center + k1_w + k1_e + k2_n + (coef + 1) * k2_s;
	        num = xi;
	        num_l = i;
	        diag[num_l] = d;
	        initial_stencil(&matrix[num_l], num, -k1_w, -k1_e, -k2_n, 0);
    	}
    }
    if (myidy == npy - 1) {
    	//north
    	for (i = 0; i < nlx; i++) {
    		xi = x0 + i;
	        k1_w = average(k1((xi + 0.5) * h, 1 - 0.5 * h), k1((xi - 0.5) * h, 1 - 0.5 * h));
	        k1_e = average(k1((xi + 0.5) * h, 1 - 0.5 * h), k1((xi + 1.5) * h, 1 - 0.5 * h));
	        k2_n = 0.5 * (3 * k2((xi + 0.5) * h, 1 - 0.5 * h) - k2((xi + 0.5) * h, 1 - 1.5 * h));
	        k2_s = average(k2((xi + 0.5) * h, 1 - 0.5 * h), k2((xi + 0.5) * h, 1 - 1.5 * h));
	        d = center + k1_w + k1_e + (coef + 1) * k2_n + k2_s;
	        num = size * size - size + xi;
	        num_l = nlx * nly - nlx + i;
	        diag[num_l] = d;
	        initial_stencil(&matrix[num_l], num, -k1_w, -k1_e, 0, -k2_s);
    	}
    }

    //vertex
    if (myidx == 0 && myidy == 0) {
    	//W-S
	    k1_w = 0.5 * (3 * k1(0.5 * h, 0.5 * h) - k1(1.5 * h, 0.5 * h));
	    k1_e = average(k1(0.5 * h, 0.5 * h), k1(1.5 * h, 0.5 * h));
	    k2_n = average(k2(0.5 * h, 0.5 * h), k2(0.5 * h, 1.5 * h));
	    k2_s = 0.5 * (3 * k2(0.5 * h, 0.5 * h) - k2(0.5 * h, 1.5 * h));
	    d = center + (coef + 1) * k1_w + k1_e + k2_n + (coef + 1) * k2_s;
	    num = 0; 
	    num_l = 0;
	    diag[num_l] = d;
	    initial_stencil(&matrix[num_l], num, 0, -k1_e, -k2_n, 0);    	
    }
    if (myidx == 0 && myidy == npy - 1) {
    	//W-N
	    k1_w = 0.5 * (3 * k1(0.5 * h, 1 - 0.5 * h) - k1(1.5 * h, 1 - 0.5 * h));
		k1_e = average(k1(0.5 * h, 1 - 0.5 * h), k1(1.5 * h, 1 - 0.5 * h));
		k2_n = 0.5 * (3 * k2(0.5 * h, 1 - 0.5 * h) - k2(0.5 * h, 1 - 1.5 * h));
		k2_s = average(k2(0.5 * h, 1 - 0.5 * h), k2(0.5 * h, 1 - 1.5 * h));
		d = center + (coef + 1) * k1_w + k1_e + (coef + 1) * k2_n + k2_s;
		num = size * size - size; 
		num_l = nlx * nly - nlx;
		diag[num_l] = d;
		initial_stencil(&matrix[num_l], num, 0, -k1_e, 0, -k2_s);
    }
    if (myidx == npx - 1 && myidy == 0) {
    	//E-S
	    k1_w = average(k1(1 - 0.5 * h, 0.5 * h), k1(1 - 1.5 * h, 0.5 * h));
	    k1_e = 0.5 * (3 * k1(1 - 0.5 * h, 0.5 * h) - k1(1 - 1.5 * h, 0.5 * h));
	    k2_n = average(k2(1 - 0.5 * h, 0.5 * h), k2(1 - 0.5 * h, 1.5 * h));
	    k2_s = 0.5 * (3 * k2(1 - 0.5 * h, 0.5 * h) - k2(1 - 0.5 * h, 1.5 * h));
	    d = center + k1_w + (coef + 1) * k1_e + k2_n + (coef + 1) * k2_s;
	    num = size -  1; 
	    num_l = nlx - 1;
	    diag[num_l] = d;
	    initial_stencil(&matrix[num_l], num, -k1_w, 0, -k2_n, 0);       	
    }
    if (myidx == npx - 1 && myidy == npy - 1) {
    	//E-N
 	    k1_w = average(k1(1 - 0.5 * h, 1 - 0.5 * h), k1(1 - 1.5 * h, 1 - 0.5 * h));
	    k1_e = 0.5 * (3 * k1(1 - 0.5 * h, 1 - 0.5 * h) - k1(1 - 1.5 * h, 1 - 0.5 * h));
	    k2_n = 0.5 * (3 * k2(1 - 0.5 * h, 1 - 0.5 * h) - k2(1 - 0.5 * h, 1 - 1.5 * h));
	    k2_s = average(k2(1 - 0.5 * h, 1 - 0.5 * h), k2(1 - 0.5 * h, 1 - 1.5 * h));
	    d = center + k1_w + (coef + 1) * k1_e + (coef + 1) * k2_n + k2_s;
	    num = size * size -  1; 
	    num_l = nlx * nly - 1;
	    diag[num_l] = d;
	    initial_stencil(&matrix[num_l], num, -k1_w, 0, 0, -k2_s);       	
    }

    return ;
}

// void matrix_stencil_mixed(Stencil *matrix, double *diag, double (*k1)(double, double), 
//                       double (*k2)(double, double), double tau, int size, double h,
//                       int nlx, int nly, int myidx, int myidy, int npx, int npy) {
//     double center = h * h * tau;
//     double k1_w, k1_e, k2_n, k2_s;
//     double d = 0;			//diag
//     int i, j, num_l;		//local
//     int xi, yj, num;		//global

//     int x0 = myidx * nlx;
//     int y0 = myidy * nly;

//     //inner mesh
//     for (j = 0; j < nly; j++) {
//     	yj = y0 + j;
//         for (i = 0; i < nlx; i++) {
//         	xi = x0 + i;
//             k1_w = average(k1((xi + 0.5) * h, (yj + 0.5) * h), k1((xi - 0.5) * h, (yj + 0.5) * h));
//             k1_e = average(k1((xi + 0.5) * h, (yj + 0.5) * h), k1((xi + 1.5) * h, (yj + 0.5) * h));
//             k2_n = average(k2((xi + 0.5) * h, (yj + 0.5) * h), k2((xi + 0.5) * h, (yj + 1.5) * h));
//             k2_s = average(k2((xi + 0.5) * h, (yj + 0.5) * h), k2((xi + 0.5) * h, (yj - 0.5) * h));
//             d = center + k1_w + k1_e + k2_n + k2_s;
//             num = yj * size + xi;
//             num_l = j * nlx + i;
//             diag[num_l] = d;
//             initial_stencil(&matrix[num_l], num, -k1_w, -k1_e, -k2_n, -k2_s);
//         }
//     }

//     //side
//     if (myidx == 0) {
// 	    // west
//     	for (j = 0; j < nly; j++) {
//     		yj = y0 + j;
// 	        k1_w = 0.5 * (3 * k1(0.5 * h, (yj + 0.5) * h) - k1(1.5 * h, (yj + 0.5) * h));
// 	        k1_e = average(k1(0.5 * h, (yj + 0.5) * h), k1(1.5 * h, (yj + 0.5) * h));
// 	        k2_n = average(k2(0.5 * h, (yj + 0.5) * h), k2(0.5 * h, (yj + 1.5) * h));
// 	        k2_s = average(k2(0.5 * h, (yj + 0.5) * h), k2(0.5 * h, (yj - 0.5) * h));
// 	        d = center + 2 * k1_w + k1_e + k2_n + k2_s;
// 	        num = size * yj;
// 	        num_l = nlx * j;
// 	        diag[num_l] = d;
// 	        initial_stencil(&matrix[num_l], num, 0, -k1_e, -k2_n, -k2_s);
//     	}
//     }
//     if (myidx == npx - 1) {
//     	//east
//     	for (j = 0; j < nly; j++) {
//     		yj = y0 + j;
// 	        k1_w = average(k1(1 - 0.5 * h, (yj + 0.5) * h), k1(1 - 1.5 * h, (yj + 0.5) * h));
// 	        k2_n = average(k2(1 - 0.5 * h, (yj + 0.5) * h), k2(1 - 0.5 * h, (yj + 1.5) * h));
// 	        k2_s = average(k2(1 - 0.5 * h, (yj + 0.5) * h), k2(1 - 0.5 * h, (yj - 0.5) * h));
// 	        d = center + k1_e + k2_n + k2_s;
// 	        num = size * yj + size - 1;
// 	        num_l = nlx * j + nlx - 1;
// 	        diag[num_l] = d;
// 	        initial_stencil(&matrix[num_l], num, -k1_w, 0, -k2_n, -k2_s);
//     	}
//     }
//     if (myidy == 0) {
//     	//south
//     	for (i = 0; i < nlx; i++) {
//     		xi = x0 + i;
// 	        k1_w = average(k1((xi + 0.5) * h, 0.5 * h), k1((xi - 0.5) * h, 0.5 * h));
// 	        k1_e = average(k1((xi + 0.5) * h, 0.5 * h), k1((xi + 1.5) * h, 0.5 * h));
// 	        k2_n = average(k2((xi + 0.5) * h, 0.5 * h), k2((xi + 0.5) * h, 1.5 * h));
// 	        d = center + k1_w + k1_e + k2_n;
// 	        num = xi;
// 	        num_l = i;
// 	        diag[num_l] = d;
// 	        initial_stencil(&matrix[num_l], num, -k1_w, -k1_e, -k2_n, 0);
//     	}
//     }
//     if (myidy == npy - 1) {
//     	//north
//     	for (i = 0; i < nlx; i++) {
//     		xi = x0 + i;
// 	        k1_w = average(k1((xi + 0.5) * h, 1 - 0.5 * h), k1((xi - 0.5) * h, 1 - 0.5 * h));
// 	        k1_e = average(k1((xi + 0.5) * h, 1 - 0.5 * h), k1((xi + 1.5) * h, 1 - 0.5 * h));
// 	        k2_s = average(k2((xi + 0.5) * h, 1 - 0.5 * h), k2((xi + 0.5) * h, 1 - 1.5 * h));
// 	        d = center + k1_w + k1_e + k2_n;
// 	        num = size * size - size + xi;
// 	        num_l = nlx * nly - nlx + i;
// 	        diag[num_l] = d;
// 	        initial_stencil(&matrix[num_l], num, -k1_w, -k1_e, 0, -k2_s);
//     	}
//     }

//     //vertex
//     if (myidx == 0 && myidy == 0) {
//     	//W-S
// 	    k1_w = 0.5 * (3 * k1(0.5 * h, 0.5 * h) - k1(1.5 * h, 0.5 * h));
// 	    k1_e = average(k1(0.5 * h, 0.5 * h), k1(1.5 * h, 0.5 * h));
// 	    k2_n = average(k2(0.5 * h, 0.5 * h), k2(0.5 * h, 1.5 * h));
// 	    d = center + 2 * k1_w + k1_e + k2_n;
// 	    num = 0; 
// 	    num_l = 0;
// 	    diag[num_l] = d;
// 	    initial_stencil(&matrix[num_l], num, 0, -k1_e, -k2_n, 0);    	
//     }
//     if (myidx == 0 && myidy == npy - 1) {
//     	//W-N
// 	    k1_w = 0.5 * (3 * k1(0.5 * h, 1 - 0.5 * h) - k1(1.5 * h, 1 - 0.5 * h));
// 		k1_e = average(k1(0.5 * h, 1 - 0.5 * h), k1(1.5 * h, 1 - 0.5 * h));
// 		k2_s = average(k2(0.5 * h, 1 - 0.5 * h), k2(0.5 * h, 1 - 1.5 * h));
// 		d = center + 2 * k1_w + k1_e + k2_s;
// 		num = size * size - size; 
// 		num_l = nlx * nly - nlx;
// 		diag[num_l] = d;
// 		initial_stencil(&matrix[num_l], num, 0, -k1_e, 0, -k2_s);
//     }
//     if (myidx == npx - 1 && myidy == 0) {
//     	//E-S
// 	    k1_w = average(k1(1 - 0.5 * h, 0.5 * h), k1(1 - 1.5 * h, 0.5 * h));
// 	    k2_n = average(k2(1 - 0.5 * h, 0.5 * h), k2(1 - 0.5 * h, 1.5 * h));
// 	    d = center + k1_w + k2_n;
// 	    num = size -  1; 
// 	    num_l = nlx - 1;
// 	    diag[num_l] = d;
// 	    initial_stencil(&matrix[num_l], num, -k1_w, 0, -k2_n, 0);       	
//     }
//     if (myidx == npx - 1 && myidy == npy - 1) {
//     	//E-N
//  	    k1_w = average(k1(1 - 0.5 * h, 1 - 0.5 * h), k1(1 - 1.5 * h, 1 - 0.5 * h));
// 	    k2_s = average(k2(1 - 0.5 * h, 1 - 0.5 * h), k2(1 - 0.5 * h, 1 - 1.5 * h));
// 	    d = center + k1_w + k2_s;
// 	    num = size * size -  1; 
// 	    num_l = nlx * nly - 1;
// 	    diag[num_l] = d;
// 	    initial_stencil(&matrix[num_l], num, -k1_w, 0, 0, -k2_s);       	
//     }

//     return ;
// }


void RHS_vector_Dirichlet(double *b, double (*k1)(double, double), double (*k2)(double, double),
                          double (*f)(double, double), double (*g)(double, double), int size, double h,
                          int nlx, int nly, int myidx, int myidy, int npx, int npy) {
	int i, j, num_l;		//local
	int xi, yj;				//global

    int x0 = myidx * nlx;
    int y0 = myidy * nly;
    double k1_w, k1_e, k2_n, k2_s;
    double h_2 = h * h;

    //inner mesh
    for (j = 0; j < nly; j++) {
    	yj = y0 + j;
    	for (i = 0; i < nlx; i++) {
    		xi = x0 + i;
    		num_l = j * nlx + i;
            b[num_l] = h_2 * f((xi + 0.5) * h ,(yj + 0.5) * h);
        }
    }

    //side
    if (myidx == 0) {
    	//west
    	for (j = 0; j < nly; j++) {
    		yj = y0 + j;
    		k1_w = 0.5 * (3 * k1(0.5 * h, (yj + 0.5) * h) - k1(1.5 * h, (yj + 0.5) * h));
	        num_l = j * nlx;
	        b[num_l] = h_2 * f(0.5 * h , (yj + 0.5) * h) + 2 * k1_w * g(0, (yj + 0.5) * h);
    	}
    }
    if (myidx == npx - 1) {
    	//east
    	for (j = 0; j < nly; j++) {
    		yj = y0 + j;
	        k1_e = 0.5 * (3 * k1(1 - 0.5 * h, 1 - 0.5 * h) - k1(1 - 1.5 * h, 1 - 0.5 * h));
	        num_l = j * nlx + nlx - 1;
	        b[num_l] = h_2 * f(1 - 0.5 * h , (yj + 0.5) * h) + 2 * k1_e * g(1, (yj + 0.5) * h);
    	}
    }
    if (myidy == 0) {
    	//south
    	for (i = 0; i < nlx; i++) {
    		xi = x0 + i;
	        k2_s = 0.5 * (3 * k2((xi + 0.5) * h, 0.5 * h) - k2((xi + 0.5) * h, 1.5 * h));
	        num_l = i;
	        b[num_l] = h_2 * f((xi + 0.5) * h , 0.5 * h) + 2 * k2_s * g((xi + 0.5) * h , 0);
    	}
    }
    if (myidy == npy - 1) {
    	//north
    	for (i = 0; i < nlx; i++) {
    		xi = x0 + i;
    		k2_n = 0.5 * (3 * k2((xi + 0.5) * h, 1 - 0.5 * h) - k2((xi + 0.5) * h, 1 - 1.5 * h));
	        num_l = nlx * nly - nlx + i;
	        b[num_l] = h_2 * f((xi + 0.5) * h , 1 - 0.5 * h) + 2 * k2_n * g((xi + 0.5) * h , 1);
    	}
    }

    //vertex
    if (myidx == 0 && myidy == 0) {
    	//W-S
	    k1_w = 0.5 * (3 * k1(0.5 * h, 0.5 * h) - k1(1.5 * h, 0.5 * h));
	    k2_s = 0.5 * (3 * k2(0.5 * h, 0.5 * h) - k2(0.5 * h, 1.5 * h));
	    num_l = 0;
	    b[num_l] = h_2 * f(0.5 * h, 0.5 * h) + 2 * (k1_w * g(0, 0.5 * h) + k2_s * g(0.5 * h, 0));  	
    }
    if (myidx == 0 && myidy == npy - 1) {
    	//W-N
	    k1_w = 0.5 * (3 * k1(0.5 * h, 1 - 0.5 * h) - k1(1.5 * h, 1 - 0.5 * h));
		k2_n = 0.5 * (3 * k2(0.5 * h, 1 - 0.5 * h) - k2(0.5 * h, 1 - 1.5 * h));
		num_l = nlx * nly - nlx;
		b[num_l] = h_2 * f(0.5 * h, 1 - 0.5 * h) + 2 * (k1_w * g(0, 1 - 0.5 * h) + k2_n * g(0.5 * h, 1));
    }
    if (myidx == npx - 1 && myidy == 0) {
    	//E-S
	    k1_e = 0.5 * (3 * k1(1 - 0.5 * h, 0.5 * h) - k1(1 - 1.5 * h, 0.5 * h));
	    k2_s = 0.5 * (3 * k2(1 - 0.5 * h, 0.5 * h) - k2(1 - 0.5 * h, 1.5 * h));
	    num_l = nlx - 1;
	    b[num_l] = h_2 * f(1 - 0.5 * h, 0.5 * h) + 2 * (k1_e * g(1, 0.5 * h) + k2_s * g(1 - 0.5 * h, 0));
    }
    if (myidx == npx - 1 && myidy == npy - 1) {
    	//E-N
	    k1_e = 0.5 * (3 * k1(1 - 0.5 * h, 1 - 0.5 * h) - k1(1 - 1.5 * h, 1 - 0.5 * h));
	    k2_n = 0.5 * (3 * k2(1 - 0.5 * h, 1 - 0.5 * h) - k2(1 - 0.5 * h, 1 - 1.5 * h));
	    num_l = nlx * nly - 1;
	    b[num_l] = h_2 * f(1 - 0.5 * h, 1 - 0.5 * h) + 2 * (k1_e * g(1, 1 - 0.5 * h) + k2_n * g(1 - 0.5 * h, 1));
    }

    return ;
}


void RHS_vector_Neumann(double *b, double (*f)(double, double), int size, double h,
						int nlx, int nly, int myidx, int myidy) {
	int i, j, num_l;		//local
	int xi, yj;				//global

    int x0 = myidx * nlx;
    int y0 = myidy * nly;
    double h_2 = h * h;

    for (j = 0; j < nly; j++) {
    	yj = y0 + j;
    	for (i = 0; i < nlx; i++) {
    		xi = x0 + i;
    		num_l = j * nlx + i;
    		b[num_l] = h_2 * f((xi + 0.5) * h, (yj + 0.5) * h);
    	}
    }

    return ;
}

void RHS_vector_mixed(double *b, double (*k1)(double, double), double (*k2)(double, double), 
                      double (*f)(double, double), double (*g)(double, double), 
                      double a1, double a2, int size, double h,
                      int nlx, int nly, int myidx, int myidy, int npx, int npy) {
	int i, j, num_l;		//local
	int xi, yj;				//global

    int x0 = myidx * nlx;
    int y0 = myidy * nly;
    double k1_w, k1_e, k2_n, k2_s;
    double h_2 = h * h;
    double coef = 2 / (a1 + 2 * a2 / h);

    //inner mesh
    for (j = 0; j < nly; j++) {
    	yj = y0 + j;
    	for (i = 0; i < nlx; i++) {
    		xi = x0 + i;
    		num_l = j * nlx + i;
            b[num_l] = h_2 * f((xi + 0.5) * h ,(yj + 0.5) * h);
        }
    }

    //side
    if (myidx == 0) {
    	//west
    	for (j = 0; j < nly; j++) {
    		yj = y0 + j;
    		k1_w = 0.5 * (3 * k1(0.5 * h, (yj + 0.5) * h) - k1(1.5 * h, (yj + 0.5) * h));
	        num_l = j * nlx;
	        b[num_l] = h_2 * f(0.5 * h , (yj + 0.5) * h) + coef * k1_w * g(0, (yj + 0.5) * h);
    	}
    }
    if (myidx == npx - 1) {
    	//east
    	for (j = 0; j < nly; j++) {
    		yj = y0 + j;
	        k1_e = 0.5 * (3 * k1(1 - 0.5 * h, 1 - 0.5 * h) - k1(1 - 1.5 * h, 1 - 0.5 * h));
	        num_l = j * nlx + nlx - 1;
	        b[num_l] = h_2 * f(1 - 0.5 * h , (yj + 0.5) * h) + coef * k1_e * g(1, (yj + 0.5) * h);
    	}
    }
    if (myidy == 0) {
    	//south
    	for (i = 0; i < nlx; i++) {
    		xi = x0 + i;
	        k2_s = 0.5 * (3 * k2((xi + 0.5) * h, 0.5 * h) - k2((xi + 0.5) * h, 1.5 * h));
	        num_l = i;
	        b[num_l] = h_2 * f((xi + 0.5) * h , 0.5 * h) + coef * k2_s * g((xi + 0.5) * h , 0);
    	}
    }
    if (myidy == npy - 1) {
    	//north
    	for (i = 0; i < nlx; i++) {
    		xi = x0 + i;
    		k2_n = 0.5 * (3 * k2((xi + 0.5) * h, 1 - 0.5 * h) - k2((xi + 0.5) * h, 1 - 1.5 * h));
	        num_l = nlx * nly - nlx + i;
	        b[num_l] = h_2 * f((xi + 0.5) * h , 1 - 0.5 * h) + coef * k2_n * g((xi + 0.5) * h , 1);
    	}
    }

        //vertex
    if (myidx == 0 && myidy == 0) {
    	//W-S
	    k1_w = 0.5 * (3 * k1(0.5 * h, 0.5 * h) - k1(1.5 * h, 0.5 * h));
	    k2_s = 0.5 * (3 * k2(0.5 * h, 0.5 * h) - k2(0.5 * h, 1.5 * h));
	    num_l = 0;
	    b[num_l] = h_2 * f(0.5 * h, 0.5 * h) + coef * (k1_w * g(0, 0.5 * h) + k2_s * g(0.5 * h, 0));  	
    }
    if (myidx == 0 && myidy == npy - 1) {
    	//W-N
	    k1_w = 0.5 * (3 * k1(0.5 * h, 1 - 0.5 * h) - k1(1.5 * h, 1 - 0.5 * h));
		k2_n = 0.5 * (3 * k2(0.5 * h, 1 - 0.5 * h) - k2(0.5 * h, 1 - 1.5 * h));
		num_l = nlx * nly - nlx;
		b[num_l] = h_2 * f(0.5 * h, 1 - 0.5 * h) + coef * (k1_w * g(0, 1 - 0.5 * h) + k2_n * g(0.5 * h, 1));
    }
    if (myidx == npx - 1 && myidy == 0) {
    	//E-S
	    k1_e = 0.5 * (3 * k1(1 - 0.5 * h, 0.5 * h) - k1(1 - 1.5 * h, 0.5 * h));
	    k2_s = 0.5 * (3 * k2(1 - 0.5 * h, 0.5 * h) - k2(1 - 0.5 * h, 1.5 * h));
	    num_l = nlx - 1;
	    b[num_l] = h_2 * f(1 - 0.5 * h, 0.5 * h) + coef * (k1_e * g(1, 0.5 * h) + k2_s * g(1 - 0.5 * h, 0));
    }
    if (myidx == npx - 1 && myidy == npy - 1) {
    	//E-N
	    k1_e = 0.5 * (3 * k1(1 - 0.5 * h, 1 - 0.5 * h) - k1(1 - 1.5 * h, 1 - 0.5 * h));
	    k2_n = 0.5 * (3 * k2(1 - 0.5 * h, 1 - 0.5 * h) - k2(1 - 0.5 * h, 1 - 1.5 * h));
	    num_l = nlx * nly - 1;
	    b[num_l] = h_2 * f(1 - 0.5 * h, 1 - 0.5 * h) + coef * (k1_e * g(1, 1 - 0.5 * h) + k2_n * g(1 - 0.5 * h, 1));
    }

    return ;
}


// void RHS_vector_mixed(double *b, double (*k1)(double, double), double (*f)(double, double),
//                       double (*g)(double, double), int size, double h,
//                       int nlx, int nly, int myidx, int myidy, int npx, int npy) {
// 	int i, j, num_l;		//local
// 	int xi, yj;				//global

//     int x0 = myidx * nlx;
//     int y0 = myidy * nly;
//     double k1_w;
//     double h_2 = h * h;

//     //inner mesh
//     for (j = 0; j < nly; j++) {
//     	yj = y0 + j;
//     	for (i = 0; i < nlx; i++) {
//     		xi = x0 + i;
//     		num_l = j * nlx + i;
//             b[num_l] = h_2 * f((xi + 0.5) * h ,(yj + 0.5) * h);
//         }
//     }

//     //side
//     if (myidx == 0) {
//     	//west
//     	for (j = 0; j < nly; j++) {
//     		yj = y0 + j;
//     		k1_w = 0.5 * (3 * k1(0.5 * h, (yj + 0.5) * h) - k1(1.5 * h, (yj + 0.5) * h));
// 	        num_l = j * nlx;
// 	        b[num_l] = h_2 * f(0.5 * h , (yj + 0.5) * h) + 2 * k1_w * g(0, (yj + 0.5) * h);
//     	}
//     }

//     return ;
// }







