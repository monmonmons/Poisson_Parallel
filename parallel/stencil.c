/*******************************************************
    > File   : stencil.c
    > Author : Yuntong
    > Mail   : 171840067@smail.nju.edu.cn
    > Date   : 2021/3/30
*******************************************************/

#include <stdio.h>
#include "stencil.h"

void initial_stencil(Stencil *stencil, int num, double v1, double v2, double v3, double v4) {
    stencil->number = num;
    stencil->value[0] = v1;
    stencil->value[1] = v2;
    stencil->value[2] = v3;
    stencil->value[3] = v4;
    return ;
}

// void print_stencil(Stencil *stencil){
//     printf("num   : %d\n", stencil->number);
//     printf("west  : %lf\n", stencil->value[0]);
//     printf("east  : %lf\n", stencil->value[1]);
//     printf("north : %lf\n", stencil->value[2]);
//     printf("south : %lf\n", stencil->value[3]);
//     return;
// }

void print_line(Stencil *stencil) {
    printf("num   : %d  value:", stencil->number);
    int i;
    for(i = 0; i < 4; i++) {
        printf("%lf ", stencil->value[i]);
    }
    printf("\n");
    return;
}

void MV_stencil(Stencil *matrix, double *vec, double *west, double *east, double *north, double *south, 
                double *result, int nlx, int nly) {
    int i, j, num;
    double temp[4];
    for (j = 0; j < nly; j++) {
        for (i = 0; i < nlx; i++) {
            num = j * nlx + i;
            temp[0] = (i == 0) ? west[j] : vec[num - 1];             //west
            temp[1] = (i == nlx - 1) ? east[j] : vec[num + 1];       //east
            temp[2] = (j == nly - 1) ? north[i] : vec[num + nlx];    //north
            temp[3] = (j == 0) ? south[i] : vec[num - nlx];          //south
            
            result[num] = dot(matrix[num].value, temp, 4);
        }
    }
    return ;
}

void MV_stencil_zero(Stencil *matrix, double *vec, double *result, int nlx, int nly) {
    int i, j, num;
    double temp[4];
    for (j = 0; j < nly; j++) {
        for (i = 0; i < nlx; i++) {
            num = j * nlx + i;
            temp[0] = (i == 0) ? 0 : vec[num - 1];             //west
            temp[1] = (i == nlx - 1) ? 0 : vec[num + 1];       //east
            temp[2] = (j == nly - 1) ? 0 : vec[num + nlx];     //north
            temp[3] = (j == 0) ? 0 : vec[num - nlx];           //south
            
            result[num] = dot(matrix[num].value, temp, 4);
        }
    }
    return ;
}

