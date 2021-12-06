/*******************************************************
    > File   : poisson.c
    > Author : Yuntong
    > Mail   : 171840067@smail.nju.edu.cn
    > Date   : 2021/3/18
*******************************************************/

#include <stdio.h>
#include "poisson.h"

void matrix_stencil_mode_const(Stencil *matrix, double *diag,
                        double k1, double k2, double tau, int size, double h, int bc) {
    double center = 2 * (k1 + k2) +  h * h * tau;
    double d = 0.0;
    int i, j, num;
    
    //side
    for (i = 1; i < size - 1; i++) {
        //south
        switch (bc) {
            case 1: d = center + k2; break;
            case 2: d = center - k2; break;
            default: break;
        }
        num = i;
        diag[num] = d;
        initial_stencil(&matrix[num], num, -k1, -k1, -k2, 0);

        //north
        num = size * size - size + i;
        diag[num] = d;
        initial_stencil(&matrix[num], num, -k1, -k1, 0, -k2);

        //west
        switch (bc) {
            case 1: d = center + k1; break;
            case 2: d = center - k1; break;
            default: break;
        }
        num = size * i;
        diag[num] = d;
        initial_stencil(&matrix[num], num, 0, -k1, -k2, -k2);

        //east
        num = size * i + size - 1;
        diag[num] = d;
        initial_stencil(&matrix[num], num, -k1, 0, -k2, -k2);
    }
    
    //vertex
    switch (bc) {
        case 1: d = center + k1 + k2; break;
        case 2: d = center - k1 - k2; break;
        default: break;
    }
    //W-S
    num = 0; diag[num] = d;
    initial_stencil(&matrix[num], num, 0, -k1, -k2, 0);
    //E-S
    num = size - 1; diag[num] = d;
    initial_stencil(&matrix[num], num, -k1, 0, -k2, 0);
    
    //W-N
    num = size * size - size; diag[num] = d;
    initial_stencil(&matrix[num], num, 0, -k1, 0, -k2);
    //E-N
    num = size * size - 1; diag[num] = d;
    initial_stencil(&matrix[num], num, -k1, 0, 0, -k2);
    
    //inner mesh
    for (j = 1; j < size - 1; j++) {
        for (i = 1; i < size - 1; i++) {
            d = center;
            num = j * size + i;
            diag[num] = d;
            initial_stencil(&matrix[num], num, -k1, -k1, -k2, -k2);
        }
    }
    
    return ;
}

void matrix_stencil_mode(Stencil *matrix, double *diag,
                          double (*k1)(double, double), double (*k2)(double, double),
                          double tau, int size, double h, int bc) {
    double center = h * h * tau;
    double k1_w, k1_e, k2_n, k2_s;
    double d = 0;
    int i, j, num;
    
    //side
    for (i = 1; i < size - 1; i++) {
        //south -- i = i; j = 0
        k1_w = average(k1((i + 0.5) * h, 0.5 * h), k1((i - 0.5) * h, 0.5 * h));
        k1_e = average(k1((i + 0.5) * h, 0.5 * h), k1((i + 1.5) * h, 0.5 * h));
        k2_n = average(k2((i + 0.5) * h, 0.5 * h), k2((i + 0.5) * h, 1.5 * h));
        k2_s = 0.5 * (3 * k2((i + 0.5) * h, 0.5 * h) - k2((i + 0.5) * h, 1.5 * h));
        switch (bc) {
            case 1: d = center + k1_w + k1_e + k2_n + 2 * k2_s; break;
            case 2: d = center + k1_w + k1_e + k2_n; break;
            default: break;
        }
        num = i;
        diag[num] = d;
        initial_stencil(&matrix[num], num, -k1_w, -k1_e, -k2_n, 0);

        //north -- i = i; j = size - 1
        k1_w = average(k1((i + 0.5) * h, 1 - 0.5 * h), k1((i - 0.5) * h, 1 - 0.5 * h));
        k1_e = average(k1((i + 0.5) * h, 1 - 0.5 * h), k1((i + 1.5) * h, 1 - 0.5 * h));
        k2_n = 0.5 * (3 * k2((i + 0.5) * h, 1 - 0.5 * h) - k2((i + 0.5) * h, 1 - 1.5 * h));
        k2_s = average(k2((i + 0.5) * h, 1 - 0.5 * h), k2((i + 0.5) * h, 1 - 1.5 * h));
        switch (bc) {
            case 1: d = center + k1_w + k1_e + 2 * k2_n + k2_s; break;
            case 2: d = center + k1_w + k1_e + k2_s; break;
            default: break;
        }
        num = size * size - size + i;
        diag[num] = d;
        initial_stencil(&matrix[num], num, -k1_w, -k1_e, 0, -k2_s);

        //west -- i = 0;  j = i
        k1_w = 0.5 * (3 * k1(0.5 * h, (i + 0.5) * h) - k1(1.5 * h, (i + 0.5) * h));
        k1_e = average(k1(0.5 * h, (i + 0.5) * h), k1(1.5 * h, (i + 0.5) * h));
        k2_n = average(k2(0.5 * h, (i + 0.5) * h), k2(0.5 * h, (i + 1.5) * h));
        k2_s = average(k2(0.5 * h, (i + 0.5) * h), k2(0.5 * h, (i - 0.5) * h));
        switch (bc) {
            case 1: d = center + 2 * k1_w + k1_e + k2_n + k2_s; break;
            case 2: d = center + k1_e + k2_n + k2_s; break;
            default: break;
        }
        num = size * i;
        diag[num] = d;
        initial_stencil(&matrix[num], num, 0, -k1_e, -k2_n, -k2_s);

        //east -- i = size - 1; j = i
        k1_w = average(k1(1 - 0.5 * h, (i + 0.5) * h), k1(1 - 1.5 * h, (i + 0.5) * h));
        k1_e = 0.5 * (3 * k1(1 - 0.5 * h, (i + 0.5) * h) - k1(1 - 1.5 * h, (i + 0.5) * h));
        k2_n = average(k2(1 - 0.5 * h, (i + 0.5) * h), k2(1 - 0.5 * h, (i + 1.5) * h));
        k2_s = average(k2(1 - 0.5 * h, (i + 0.5) * h), k2(1 - 0.5 * h, (i - 0.5) * h));
        switch (bc) {
            case 1: d = center + k1_w + 2 * k1_e + k2_n + k2_s; break;
            case 2: d = center + k1_w + k2_n + k2_s; break;
            default: break;
        }
        num = size * i + size - 1;
        diag[num] = d;
        initial_stencil(&matrix[num], num, -k1_w, 0, -k2_n, -k2_s);
    }
    
    //vertex
    //W-S -- i = 0; j = 0
    k1_w = 0.5 * (3 * k1(0.5 * h, 0.5 * h) - k1(1.5 * h, 0.5 * h));
    k1_e = average(k1(0.5 * h, 0.5 * h), k1(1.5 * h, 0.5 * h));
    k2_n = average(k2(0.5 * h, 0.5 * h), k2(0.5 * h, 1.5 * h));
    k2_s = 0.5 * (3 * k2(0.5 * h, 0.5 * h) - k2(0.5 * h, 1.5 * h));
    switch (bc) {
        case 1: d = center + 2 * k1_w + k1_e + k2_n + 2 * k2_s; break;
        case 2: d = center + k1_e + k2_n; break;
        default: break;
    }
    num = 0; diag[num] = d;
    initial_stencil(&matrix[num], num, 0, -k1_e, -k2_n, 0);
    //E-S -- i = size - 1; j = 0
    k1_w = average(k1(1 - 0.5 * h, 0.5 * h), k1(1 - 1.5 * h, 0.5 * h));
    k1_e = 0.5 * (3 * k1(1 - 0.5 * h, 0.5 * h) - k1(1 - 1.5 * h, 0.5 * h));
    k2_n = average(k2(1 - 0.5 * h, 0.5 * h), k2(1 - 0.5 * h, 1.5 * h));
    k2_s = 0.5 * (3 * k2(1 - 0.5 * h, 0.5 * h) - k2(1 - 0.5 * h, 1.5 * h));
    switch (bc) {
        case 1: d = center + k1_w + 2 * k1_e + k2_n + 2 * k2_s; break;
        case 2: d = center + k1_w + k2_n; break;
        default: break;
    }
    num = size - 1; diag[num] = d;
    initial_stencil(&matrix[num], num, -k1_w, 0, -k2_n, 0);
    
    //W-N -- i = 0; j = size - 1
    k1_w = 0.5 * (3 * k1(0.5 * h, 1 - 0.5 * h) - k1(1.5 * h, 1 - 0.5 * h));
    k1_e = average(k1(0.5 * h, 1 - 0.5 * h), k1(1.5 * h, 1 - 0.5 * h));
    k2_n = 0.5 * (3 * k2(0.5 * h, 1 - 0.5 * h) - k2(0.5 * h, 1 - 1.5 * h));
    k2_s = average(k2(0.5 * h, 1 - 0.5 * h), k2(0.5 * h, 1 - 1.5 * h));
    switch (bc) {
        case 1: d = center + 2 * k1_w + k1_e + 2 * k2_n + k2_s; break;
        case 2: d = center + k1_e + k2_s; break;
        default: break;
    }
    num = size * size - size; diag[num] = d;
    initial_stencil(&matrix[num], num, 0, -k1_e, 0, -k2_s);
    //E-N -- i = size - 1; j = size - 1
    k1_w = average(k1(1 - 0.5 * h, 1 - 0.5 * h), k1(1 - 1.5 * h, 1 - 0.5 * h));
    k1_e = 0.5 * (3 * k1(1 - 0.5 * h, 1 - 0.5 * h) - k1(1 - 1.5 * h, 1 - 0.5 * h));
    k2_n = 0.5 * (3 * k2(1 - 0.5 * h, 1 - 0.5 * h) - k2(1 - 0.5 * h, 1 - 1.5 * h));
    k2_s = average(k2(1 - 0.5 * h, 1 - 0.5 * h), k2(1 - 0.5 * h, 1 - 1.5 * h));
    switch (bc) {
        case 1: d = center + k1_w + 2 * k1_e + 2 * k2_n + k2_s; break;
        case 2: d = center + k1_w + k2_s; break;
        default: break;
    }
    num = size * size - 1; diag[num] = d;
    initial_stencil(&matrix[num], num, -k1_w, 0, 0, -k2_s);
    
    //inner mesh
    for (j = 1; j < size - 1; j++) {
        for (i = 1; i < size - 1; i++) {
            k1_w = average(k1((i + 0.5) * h, (j + 0.5) * h), k1((i - 0.5) * h, (j + 0.5) * h));
            k1_e = average(k1((i + 0.5) * h, (j + 0.5) * h), k1((i + 1.5) * h, (j + 0.5) * h));
            k2_n = average(k2((i + 0.5) * h, (j + 0.5) * h), k2((i + 0.5) * h, (j + 1.5) * h));
            k2_s = average(k2((i + 0.5) * h, (j + 0.5) * h), k2((i + 0.5) * h, (j - 0.5) * h));
            d = center + k1_w + k1_e + k2_n + k2_s;
            num = j * size + i;
            diag[num] = d;
            initial_stencil(&matrix[num], num, -k1_w, -k1_e, -k2_n, -k2_s);
        }
    }
    
    return ;
}

void matrix_stencil_mixed(Stencil *matrix, double *diag, double (*k1)(double, double),
                          double (*k2)(double, double), double a1, double a2,
                          double tau, int size, double h) {
    double center = h * h * tau;
    double k1_w, k1_e, k2_n, k2_s;
    double d = 0;
    int nlx = size, nly = size;
    int i, j, num;

    double coef = (a1 - 2 * a2 / h) / (a1 + 2 * a2 / h);
    // double coef = (2 * a1) / (a1 + 2 * a2 / h);
    
    //side
    for (j = 1; j < nly - 1; j++) {
        //west -- i = 0;  j = 1 : nly - 1
        k1_w = 0.5 * (3 * k1(0.5 * h, (j + 0.5) * h) - k1(1.5 * h, (j + 0.5) * h));
        k1_e = average(k1(0.5 * h, (j + 0.5) * h), k1(1.5 * h, (j + 0.5) * h));
        k2_n = average(k2(0.5 * h, (j + 0.5) * h), k2(0.5 * h, (j + 1.5) * h));
        k2_s = average(k2(0.5 * h, (j + 0.5) * h), k2(0.5 * h, (j - 0.5) * h));
        d = center + (coef + 1) * k1_w + k1_e + k2_n + k2_s;
        num = size * j;
        diag[num] = d;
        initial_stencil(&matrix[num], num, 0, -k1_e, -k2_n, -k2_s);

        //east -- i = nly - 1; j = 1 : nly - 1
        k1_w = average(k1(1 - 0.5 * h, (j + 0.5) * h), k1(1 - 1.5 * h, (j + 0.5) * h));
        k1_e = 0.5 * (3 * k1(1 - 0.5 * h, (j + 0.5) * h) - k1(1 - 1.5 * h, (j + 0.5) * h));
        k2_n = average(k2(1 - 0.5 * h, (j + 0.5) * h), k2(1 - 0.5 * h, (j + 1.5) * h));
        k2_s = average(k2(1 - 0.5 * h, (j + 0.5) * h), k2(1 - 0.5 * h, (j - 0.5) * h));
        d = center + k1_w + (coef + 1) * k1_e + k2_n + k2_s;
        num = size * j + size - 1;
        diag[num] = d;
        initial_stencil(&matrix[num], num, -k1_w, 0, -k2_n, -k2_s);
    }

    for (i = 1; i < nlx - 1; i++) {
        //south -- i = 1 : nlx - 1; j = 0
        k1_w = average(k1((i + 0.5) * h, 0.5 * h), k1((i - 0.5) * h, 0.5 * h));
        k1_e = average(k1((i + 0.5) * h, 0.5 * h), k1((i + 1.5) * h, 0.5 * h));
        k2_n = average(k2((i + 0.5) * h, 0.5 * h), k2((i + 0.5) * h, 1.5 * h));
        k2_s = 0.5 * (3 * k2((i + 0.5) * h, 0.5 * h) - k2((i + 0.5) * h, 1.5 * h));
        d = center + k1_w + k1_e + k2_n + (coef + 1) * k2_s;
        num = i;
        diag[num] = d;
        initial_stencil(&matrix[num], num, -k1_w, -k1_e, -k2_n, 0);

        //north -- i = 1 : nlx - 1; j = size - 1
        k1_w = average(k1((i + 0.5) * h, 1 - 0.5 * h), k1((i - 0.5) * h, 1 - 0.5 * h));
        k1_e = average(k1((i + 0.5) * h, 1 - 0.5 * h), k1((i + 1.5) * h, 1 - 0.5 * h));
        k2_n = 0.5 * (3 * k2((i + 0.5) * h, 1 - 0.5 * h) - k2((i + 0.5) * h, 1 - 1.5 * h));
        k2_s = average(k2((i + 0.5) * h, 1 - 0.5 * h), k2((i + 0.5) * h, 1 - 1.5 * h));
        d = center + k1_w + k1_e + (coef + 1) * k2_n + k2_s;
        num = size * size - size + i;
        diag[num] = d;
        initial_stencil(&matrix[num], num, -k1_w, -k1_e, 0, -k2_s);
    }
    
    //vertex
    //W-S -- i = 0; j = 0
    k1_w = 0.5 * (3 * k1(0.5 * h, 0.5 * h) - k1(1.5 * h, 0.5 * h));
    k1_e = average(k1(0.5 * h, 0.5 * h), k1(1.5 * h, 0.5 * h));
    k2_n = average(k2(0.5 * h, 0.5 * h), k2(0.5 * h, 1.5 * h));
    k2_s = 0.5 * (3 * k2(0.5 * h, 0.5 * h) - k2(0.5 * h, 1.5 * h));
    d = center + (coef + 1) * k1_w + k1_e + k2_n + (coef + 1) * k2_s;
    num = 0; diag[num] = d;
    initial_stencil(&matrix[num], num, 0, -k1_e, -k2_n, 0);
    
    //W-N -- i = 0; j = size - 1
    k1_w = 0.5 * (3 * k1(0.5 * h, 1 - 0.5 * h) - k1(1.5 * h, 1 - 0.5 * h));
    k1_e = average(k1(0.5 * h, 1 - 0.5 * h), k1(1.5 * h, 1 - 0.5 * h));
    k2_n = 0.5 * (3 * k2(0.5 * h, 1 - 0.5 * h) - k2(0.5 * h, 1 - 1.5 * h));
    k2_s = average(k2(0.5 * h, 1 - 0.5 * h), k2(0.5 * h, 1 - 1.5 * h));
    d = center + (coef + 1) * k1_w + k1_e + (coef + 1) * k2_n + k2_s;
    num = size * size - size; diag[num] = d;
    initial_stencil(&matrix[num], num, 0, -k1_e, 0, -k2_s);
    
    //E-S -- i = size - 1; j = 0
    k1_w = average(k1(1 - 0.5 * h, 0.5 * h), k1(1 - 1.5 * h, 0.5 * h));
    k1_e = 0.5 * (3 * k1(1 - 0.5 * h, 0.5 * h) - k1(1 - 1.5 * h, 0.5 * h));
    k2_n = average(k2(1 - 0.5 * h, 0.5 * h), k2(1 - 0.5 * h, 1.5 * h));
    k2_s = 0.5 * (3 * k2(1 - 0.5 * h, 0.5 * h) - k2(1 - 0.5 * h, 1.5 * h));
    d = center + k1_w + (coef + 1) * k1_e + k2_n + (coef + 1) * k2_s;
    num = size - 1; diag[num] = d;
    initial_stencil(&matrix[num], num, -k1_w, 0, -k2_n, 0);

    //E-N -- i = size - 1; j = size - 1
    k1_w = average(k1(1 - 0.5 * h, 1 - 0.5 * h), k1(1 - 1.5 * h, 1 - 0.5 * h));
    k1_e = 0.5 * (3 * k1(1 - 0.5 * h, 1 - 0.5 * h) - k1(1 - 1.5 * h, 1 - 0.5 * h));
    k2_n = 0.5 * (3 * k2(1 - 0.5 * h, 1 - 0.5 * h) - k2(1 - 0.5 * h, 1 - 1.5 * h));
    k2_s = average(k2(1 - 0.5 * h, 1 - 0.5 * h), k2(1 - 0.5 * h, 1 - 1.5 * h));
    d = center + k1_w + (coef + 1) * k1_e + (coef + 1) * k2_n + k2_s;
    num = size * size - 1; diag[num] = d;
    initial_stencil(&matrix[num], num, -k1_w, 0, 0, -k2_s);
    
    //inner mesh
    for (j = 1; j < nly - 1; j++) {
        for (i = 1; i < nlx - 1; i++) {
            k1_w = average(k1((i + 0.5) * h, (j + 0.5) * h), k1((i - 0.5) * h, (j + 0.5) * h));
            k1_e = average(k1((i + 0.5) * h, (j + 0.5) * h), k1((i + 1.5) * h, (j + 0.5) * h));
            k2_n = average(k2((i + 0.5) * h, (j + 0.5) * h), k2((i + 0.5) * h, (j + 1.5) * h));
            k2_s = average(k2((i + 0.5) * h, (j + 0.5) * h), k2((i + 0.5) * h, (j - 0.5) * h));
            d = center + k1_w + k1_e + k2_n + k2_s;
            num = j * size + i;
            diag[num] = d;
            initial_stencil(&matrix[num], num, -k1_w, -k1_e, -k2_n, -k2_s);
        }
    }
    
    return ;
}

// void matrix_stencil_mixed_const(Stencil *matrix, double *diag,
//                         double k1, double k2, double tau, int size, double h) {
//     double center = 2 * (k1 + k2) +  h * h * tau;
//     double d;
//     int i, j, num;
    
//     //side
//     for (i = 1; i < size - 1; i++) {
//         //south
//         d = center - k2;
//         num = i;
//         diag[num] = d;
//         initial_stencil(&matrix[num], num, -k1, -k1, -k2, 0);

//         //north
//         num = size * size - size + i;
//         diag[num] = d;
//         initial_stencil(&matrix[num], num, -k1, -k1, 0, -k2);

//         //west
//         d = center + k1;
//         num = size * i;
//         diag[num] = d;
//         initial_stencil(&matrix[num], num, 0, -k1, -k2, -k2);

//         //east
//         d = center - k1;
//         num = size * i + size - 1;
//         diag[num] = d;
//         initial_stencil(&matrix[num], num, -k1, 0, -k2, -k2);
//     }
    
//     //vertex
//     d = center + k1 - k2;
//     //W-S
//     num = 0; diag[num] = d;
//     initial_stencil(&matrix[num], num, 0, -k1, -k2, 0);
//     //W-N
//     num = size * size - size; diag[num] = d;
//     initial_stencil(&matrix[num], num, 0, -k1, 0, -k2);
    
//     d = center - k1 - k2;
//     //E-S
//     num = size - 1; diag[num] = d;
//     initial_stencil(&matrix[num], num, -k1, 0, -k2, 0);
//     //E-N
//     num = size * size - 1; diag[num] = d;
//     initial_stencil(&matrix[num], num, -k1, 0, 0, -k2);
    
//     //inner mesh
//     for (j = 1; j < size - 1; j++) {
//         for (i = 1; i < size - 1; i++) {
//             d = center;
//             num = j * size + i;
//             diag[num] = d;
//             initial_stencil(&matrix[num], num, -k1, -k1, -k2, -k2);
//         }
//     }

//     return ;
// }

// void matrix_stencil_mixed(Stencil *matrix, double *diag, double (*k1)(double, double),
//                           double (*k2)(double, double), double tau, int size, double h) {
//     double center = h * h * tau;
//     double k1_w, k1_e, k2_n, k2_s;
//     double d = 0;
//     int i, j, num;
    
//     //side
//     for (i = 1; i < size - 1; i++) {
//         //south -- i = i; j = 0
//         k1_w = average(k1((i + 0.5) * h, 0.5 * h), k1((i - 0.5) * h, 0.5 * h));
//         k1_e = average(k1((i + 0.5) * h, 0.5 * h), k1((i + 1.5) * h, 0.5 * h));
//         k2_n = average(k2((i + 0.5) * h, 0.5 * h), k2((i + 0.5) * h, 1.5 * h));
//         d = center + k1_w + k1_e + k2_n;
//         num = i;
//         diag[num] = d;
//         initial_stencil(&matrix[num], num, -k1_w, -k1_e, -k2_n, 0);

//         //north -- i = i; j = size - 1
//         k1_w = average(k1((i + 0.5) * h, 1 - 0.5 * h), k1((i - 0.5) * h, 1 - 0.5 * h));
//         k1_e = average(k1((i + 0.5) * h, 1 - 0.5 * h), k1((i + 1.5) * h, 1 - 0.5 * h));
//         k2_s = average(k2((i + 0.5) * h, 1 - 0.5 * h), k2((i + 0.5) * h, 1 - 1.5 * h));
//         d = center + k1_w + k1_e + k2_s;
//         num = size * size - size + i;
//         diag[num] = d;
//         initial_stencil(&matrix[num], num, -k1_w, -k1_e, 0, -k2_s);

//         //west -- i = 0;  j = i
//         k1_w = 0.5 * (3 * k1(0.5 * h, (i + 0.5) * h) - k1(1.5 * h, (i + 0.5) * h));
//         k1_e = average(k1(0.5 * h, (i + 0.5) * h), k1(1.5 * h, (i + 0.5) * h));
//         k2_n = average(k2(0.5 * h, (i + 0.5) * h), k2(0.5 * h, (i + 1.5) * h));
//         k2_s = average(k2(0.5 * h, (i + 0.5) * h), k2(0.5 * h, (i - 0.5) * h));
//         d = center + 2 * k1_w + k1_e + k2_n + k2_s;
//         num = size * i;
//         diag[num] = d;
//         initial_stencil(&matrix[num], num, 0, -k1_e, -k2_n, -k2_s);

//         //east -- i = size - 1; j = i
//         k1_w = average(k1(1 - 0.5 * h, (i + 0.5) * h), k1(1 - 1.5 * h, (i + 0.5) * h));
//         k2_n = average(k2(1 - 0.5 * h, (i + 0.5) * h), k2(1 - 0.5 * h, (i + 1.5) * h));
//         k2_s = average(k2(1 - 0.5 * h, (i + 0.5) * h), k2(1 - 0.5 * h, (i - 0.5) * h));
//         d = center + k1_w + k2_n + k2_s;
//         num = size * i + size - 1;
//         diag[num] = d;
//         initial_stencil(&matrix[num], num, -k1_w, 0, -k2_n, -k2_s);
//     }
    
//     //vertex
//     //W-S -- i = 0; j = 0
//     k1_w = 0.5 * (3 * k1(0.5 * h, 0.5 * h) - k1(1.5 * h, 0.5 * h));
//     k1_e = average(k1(0.5 * h, 0.5 * h), k1(1.5 * h, 0.5 * h));
//     k2_n = average(k2(0.5 * h, 0.5 * h), k2(0.5 * h, 1.5 * h));
//     d = center + 2 * k1_w + k1_e + k2_n;
//     num = 0; diag[num] = d;
//     initial_stencil(&matrix[num], num, 0, -k1_e, -k2_n, 0);
    
//     //W-N -- i = 0; j = size - 1
//     k1_w = 0.5 * (3 * k1(0.5 * h, 1 - 0.5 * h) - k1(1.5 * h, 1 - 0.5 * h));
//     k1_e = average(k1(0.5 * h, 1 - 0.5 * h), k1(1.5 * h, 1 - 0.5 * h));
//     k2_s = average(k2(0.5 * h, 1 - 0.5 * h), k2(0.5 * h, 1 - 1.5 * h));
//     d = center + 2 * k1_w + k1_e + k2_s;
//     num = size * size - size; diag[num] = d;
//     initial_stencil(&matrix[num], num, 0, -k1_e, 0, -k2_s);
    
//     //E-S -- i = size - 1; j = 0
//     k1_w = average(k1(1 - 0.5 * h, 0.5 * h), k1(1 - 1.5 * h, 0.5 * h));
//     k2_n = average(k2(1 - 0.5 * h, 0.5 * h), k2(1 - 0.5 * h, 1.5 * h));
//     d = center + k1_w + k2_n;
//     num = size - 1; diag[num] = d;
//     initial_stencil(&matrix[num], num, -k1_w, 0, -k2_n, 0);

//     //E-N -- i = size - 1; j = size - 1
//     k1_w = average(k1(1 - 0.5 * h, 1 - 0.5 * h), k1(1 - 1.5 * h, 1 - 0.5 * h));
//     k2_s = average(k2(1 - 0.5 * h, 1 - 0.5 * h), k2(1 - 0.5 * h, 1 - 1.5 * h));
//     d = center + k1_w + k2_s;
//     num = size * size - 1; diag[num] = d;
//     initial_stencil(&matrix[num], num, -k1_w, 0, 0, -k2_s);
    
//     //inner mesh
//     for (j = 1; j < size - 1; j++) {
//         for (i = 1; i < size - 1; i++) {
//             k1_w = average(k1((i + 0.5) * h, (j + 0.5) * h), k1((i - 0.5) * h, (j + 0.5) * h));
//             k1_e = average(k1((i + 0.5) * h, (j + 0.5) * h), k1((i + 1.5) * h, (j + 0.5) * h));
//             k2_n = average(k2((i + 0.5) * h, (j + 0.5) * h), k2((i + 0.5) * h, (j + 1.5) * h));
//             k2_s = average(k2((i + 0.5) * h, (j + 0.5) * h), k2((i + 0.5) * h, (j - 0.5) * h));
//             d = center + k1_w + k1_e + k2_n + k2_s;
//             num = j * size + i;
//             diag[num] = d;
//             initial_stencil(&matrix[num], num, -k1_w, -k1_e, -k2_n, -k2_s);
//         }
//     }
    
//     return ;
// }


void RHS_vector_Dirichlet_const(double *b, double k1, double k2, double (*f)(double, double), double (*g)(double, double),
                      int size, double h) {
    int i, j, num;
    double h_2 = h * h;

    //side
    for (i = 1; i < size - 1; i++) {
        //south -- i = i; j = 0
        num = i;
        b[num] = h_2 * f((i + 0.5) * h , 0.5 * h) + 2 * k2 * g((i + 0.5) * h , 0);

        //north -- i = i; j = size - 1
        num = size * size - size + i;
        b[num] = h_2 * f((i + 0.5) * h , 1 - 0.5 * h) + 2 * k2 * g((i + 0.5) * h , 1);

        //west -- i = 0;  j = i
        num = i * size;
        b[num] = h_2 * f(0.5 * h , (i + 0.5) * h) + 2 * k1 * g(0, (i + 0.5) * h);

        //east -- i = size - 1; j = i
        num = i * size  + size - 1;
        b[num] = h_2 * f(1 - 0.5 * h , (i + 0.5) * h) + 2 * k1 * g(1, (i + 0.5) * h);
    }

    //vertex
    //W-S -- i = 0; j = 0
    num = 0;
    b[num] = h_2 * f(0.5 * h, 0.5 * h) + 2 * (k1 * g(0, 0.5 * h) + k2 * g(0.5 * h, 0));

    //E-S -- i = size - 1; j = 0
    num = size - 1;
    b[num] = h_2 * f(1 - 0.5 * h, 0.5 * h) + 2 * (k1 * g(1, 0.5 * h) + k2 * g(1 - 0.5 * h, 0));

    //W-N -- i = 0; j = size - 1
    num = size * size - size;
    b[num] = h_2 * f(0.5 * h, 1 - 0.5 * h) + 2 * (k1 * g(0, 1 - 0.5 * h) + k2 * g(0.5 * h, 1));

    //E-N -- i = size - 1; j = size - 1
    num = size * size - 1;
    b[num] = h_2 * f(1 - 0.5 * h, 1 - 0.5 * h) + 2 * (k1 * g(1, 1 - 0.5 * h) + k2 * g(1 - 0.5 * h, 1));

    //inner mesh
    for (j = 1; j < size - 1; j++) {
        for (i = 1; i < size - 1; i++) {
            num = j * size + i;
            b[num] = h_2 * f((i + 0.5) * h ,(j + 0.5) * h);
        }
    }

    return ;
}

void RHS_vector_Dirichlet(double *b, double (*k1)(double, double), double (*k2)(double, double),
                          double (*f)(double, double), double (*g)(double, double), int size, double h) {
    int nlx = size, nly = size;
    int i, j, num;
    double k1_w, k1_e, k2_n, k2_s;
    double h_2 = h * h;
    
    //side
    for (i = 1; i < size - 1; i++) {
        //south -- i = i; j = 0
        k2_s = 0.5 * (3 * k2((i + 0.5) * h, 0.5 * h) - k2((i + 0.5) * h, 1.5 * h));
        num = i;
        b[num] = h_2 * f((i + 0.5) * h , 0.5 * h) + 2 * k2_s * g((i + 0.5) * h , 0);

        //north -- i = i; j = size - 1
        k2_n = 0.5 * (3 * k2((i + 0.5) * h, 1 - 0.5 * h) - k2((i + 0.5) * h, 1 - 1.5 * h));
        num = size * size - size + i;
        b[num] = h_2 * f((i + 0.5) * h , 1 - 0.5 * h) + 2 * k2_n * g((i + 0.5) * h , 1);

        //west -- i = 0;  j = i
        k1_w = 0.5 * (3 * k1(0.5 * h, (i + 0.5) * h) - k1(1.5 * h, (i + 0.5) * h));
        num = i * size;
        b[num] = h_2 * f(0.5 * h , (i + 0.5) * h) + 2 * k1_w * g(0, (i + 0.5) * h);

        //east -- i = size - 1; j = i
        k1_e = 0.5 * (3 * k1(1 - 0.5 * h, 1 - 0.5 * h) - k1(1 - 1.5 * h, 1 - 0.5 * h));
        num = i * size  + size - 1;
        b[num] = h_2 * f(1 - 0.5 * h , (i + 0.5) * h) + 2 * k1_e * g(1, (i + 0.5) * h);
    }
    
    //vertex
    //W-S -- i = 0; j = 0
    k1_w = 0.5 * (3 * k1(0.5 * h, 0.5 * h) - k1(1.5 * h, 0.5 * h));
    k2_s = 0.5 * (3 * k2(0.5 * h, 0.5 * h) - k2(0.5 * h, 1.5 * h));
    num = 0;
    b[num] = h_2 * f(0.5 * h, 0.5 * h) + 2 * (k1_w * g(0, 0.5 * h) + k2_s * g(0.5 * h, 0));

    //E-S -- i = size - 1; j = 0
    k1_e = 0.5 * (3 * k1(1 - 0.5 * h, 0.5 * h) - k1(1 - 1.5 * h, 0.5 * h));
    k2_s = 0.5 * (3 * k2(1 - 0.5 * h, 0.5 * h) - k2(1 - 0.5 * h, 1.5 * h));
    num = size - 1;
    b[num] = h_2 * f(1 - 0.5 * h, 0.5 * h) + 2 * (k1_e * g(1, 0.5 * h) + k2_s * g(1 - 0.5 * h, 0));

    //W-N -- i = 0; j = size - 1
    k1_w = 0.5 * (3 * k1(0.5 * h, 1 - 0.5 * h) - k1(1.5 * h, 1 - 0.5 * h));
    k2_n = 0.5 * (3 * k2(0.5 * h, 1 - 0.5 * h) - k2(0.5 * h, 1 - 1.5 * h));
    num = size * size - size;
    b[num] = h_2 * f(0.5 * h, 1 - 0.5 * h) + 2 * (k1_w * g(0, 1 - 0.5 * h) + k2_n * g(0.5 * h, 1));

    //E-N -- i = size - 1; j = size - 1
    k1_e = 0.5 * (3 * k1(1 - 0.5 * h, 1 - 0.5 * h) - k1(1 - 1.5 * h, 1 - 0.5 * h));
    k2_n = 0.5 * (3 * k2(1 - 0.5 * h, 1 - 0.5 * h) - k2(1 - 0.5 * h, 1 - 1.5 * h));
    num = size * size - 1;
    b[num] = h_2 * f(1 - 0.5 * h, 1 - 0.5 * h) + 2 * (k1_e * g(1, 1 - 0.5 * h) + k2_n * g(1 - 0.5 * h, 1));
    
    //inner mesh
    for (j = 1; j < nly - 1; j++) {
        for (i = 1; i < nlx - 1; i++) {
            num = j * size + i;
            b[num] = h_2 * f((i + 0.5) * h ,(j + 0.5) * h);
        }
    }

    return ;
}


void RHS_vector_Neumann(double *b, double (*f)(double, double), int size, double h) {
    int nlx = size, nly = size;
    int i, j, num;
    double h_2 = h * h;
    
    for (j = 0; j < nly; j++) {
        for (i = 0; i < nlx; i++) {
            num = j * size + i;
            b[num] = h_2 * f((i + 0.5) * h ,(j + 0.5) * h);
        }
    }
    
    return ;
}

void RHS_vector_mixed(double *b, double (*k1)(double, double), double (*k2)(double, double),
                      double (*f)(double, double), double (*g)(double, double), 
                      double a1, double a2, int size, double h) {
    int nlx = size, nly = size;
    int i, j, num;
    double k1_w, k1_e, k2_n, k2_s;
    double h_2 = h * h;
    double coef = 2 / (a1 + 2 * a2 / h);
    
    //side
    for (j = 0; j < nly - 1; j++) {
        //west -- i = 0; j = 1 : nly - 1
        k1_w = 0.5 * (3 * k1(0.5 * h, (j + 0.5) * h) - k1(1.5 * h, (j + 0.5) * h));
        num = j * size;
        b[num] = h_2 * f(0.5 * h , (j + 0.5) * h) + coef * k1_w * g(0, (j + 0.5) * h);

        //east -- i = size - 1; j = 1 : nly - 1
        k1_e = 0.5 * (3 * k1(1 - 0.5 * h, 1 - 0.5 * h) - k1(1 - 1.5 * h, 1 - 0.5 * h));
        num = j * size  + size - 1;
        b[num] = h_2 * f(1 - 0.5 * h , (j + 0.5) * h) + coef * k1_e * g(1, (j + 0.5) * h);
    }

    for (i = 1; i < nlx - 1; i++) {
        //south -- i = 1 : nlx - 1; j = 0
        k2_s = 0.5 * (3 * k2((i + 0.5) * h, 0.5 * h) - k2((i + 0.5) * h, 1.5 * h));
        num = i;
        b[num] = h_2 * f((i + 0.5) * h , 0.5 * h) + coef * k2_s * g((i + 0.5) * h , 0);

        //north -- i = 1 : nlx - 1; j = size - 1
        k2_n = 0.5 * (3 * k2((i + 0.5) * h, 1 - 0.5 * h) - k2((i + 0.5) * h, 1 - 1.5 * h));
        num = size * size - size + i;
        b[num] = h_2 * f((i + 0.5) * h , 1 - 0.5 * h) + coef * k2_n * g((i + 0.5) * h , 1);
    }
    
    //vertex
    //W-S -- i = 0; j = 0
    k1_w = 0.5 * (3 * k1(0.5 * h, 0.5 * h) - k1(1.5 * h, 0.5 * h));
    k2_s = 0.5 * (3 * k2(0.5 * h, 0.5 * h) - k2(0.5 * h, 1.5 * h));
    num = 0;
    b[num] = h_2 * f(0.5 * h, 0.5 * h) + coef * (k1_w * g(0, 0.5 * h) + k2_s * g(0.5 * h, 0));

    //E-S -- i = size - 1; j = 0
    k1_e = 0.5 * (3 * k1(1 - 0.5 * h, 0.5 * h) - k1(1 - 1.5 * h, 0.5 * h));
    k2_s = 0.5 * (3 * k2(1 - 0.5 * h, 0.5 * h) - k2(1 - 0.5 * h, 1.5 * h));
    num = size - 1;
    b[num] = h_2 * f(1 - 0.5 * h, 0.5 * h) + coef * (k1_e * g(1, 0.5 * h) + k2_s * g(1 - 0.5 * h, 0));

    //W-N -- i = 0; j = size - 1
    k1_w = 0.5 * (3 * k1(0.5 * h, 1 - 0.5 * h) - k1(1.5 * h, 1 - 0.5 * h));
    k2_n = 0.5 * (3 * k2(0.5 * h, 1 - 0.5 * h) - k2(0.5 * h, 1 - 1.5 * h));
    num = size * size - size;
    b[num] = h_2 * f(0.5 * h, 1 - 0.5 * h) + coef * (k1_w * g(0, 1 - 0.5 * h) + k2_n * g(0.5 * h, 1));

    //E-N -- i = size - 1; j = size - 1
    k1_e = 0.5 * (3 * k1(1 - 0.5 * h, 1 - 0.5 * h) - k1(1 - 1.5 * h, 1 - 0.5 * h));
    k2_n = 0.5 * (3 * k2(1 - 0.5 * h, 1 - 0.5 * h) - k2(1 - 0.5 * h, 1 - 1.5 * h));
    num = size * size - 1;
    b[num] = h_2 * f(1 - 0.5 * h, 1 - 0.5 * h) + coef * (k1_e * g(1, 1 - 0.5 * h) + k2_n * g(1 - 0.5 * h, 1));

    //inner mesh
    for (j = 1; j < nly - 1; j++) {
        for (i = 1; i < nlx - 1; i++) {
            num = j * size + i;
            b[num] = h_2 * f((i + 0.5) * h ,(j + 0.5) * h);
        }
    }
    
    return ;
}

// void RHS_vector_mixed_const(double *b, double k1, double (*f)(double, double), double (*g)(double, double),
//                       int size, double h){
//     int i, j, num;
//     double h_2 = h * h;
    
//     //west
//     for (j = 0; j < size; j++) {
//         //i = 0;  j = j
//         num = j * size;
//         b[num]  = h_2 * f(0.5 * h , (j + 0.5) * h) + 2 * k1 * g(0, (j + 0.5) * h);
        
//     }
    
//     //inner mesh
//     for (j = 0; j < size; j++) {
//         for (i = 1; i < size; i++) {
//             num = j * size + i;
//             b[num] = h_2 * f((i + 0.5) * h ,(j + 0.5) * h);
//         }
//     }
    
//     return ;
// }

// void RHS_vector_mixed(double *b, double (*k1)(double, double), double (*f)(double, double),
//                       double (*g)(double, double), int size, double h) {
//     int i, j, num;
//     double k1_w;
//     double h_2 = h * h;
    
//     //west
//     for (j = 0; j < size; j++) {
//         //i = 0;  j = j
//         k1_w = 0.5 * (3 * k1(0.5 * h, (j + 0.5) * h) - k1(1.5 * h, (j + 0.5) * h));
//         num = j * size;
//         b[num]  = h_2 * f(0.5 * h , (j + 0.5) * h) + 2 * k1_w * g(0, (j + 0.5) * h);
        
//     }
    
//     //inner mesh
//     for (j = 0; j < size; j++) {
//         for (i = 1; i < size; i++) {
//             num = j * size + i;
//             b[num] = h_2 * f((i + 0.5) * h ,(j + 0.5) * h);
//         }
//     }
    
//     return ;
// }

