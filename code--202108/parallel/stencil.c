/*******************************************************
    > File   : stencil.c
    > Author : Yuntong
    > Mail   : 171840067@smail.nju.edu.cn
    > Date   : 2021/8/20
*******************************************************/

#include <stdio.h>
#include "stencil.h"

void initial_stencil(Stencil *stencil, double v1, double v2, double v3, double v4) {
    stencil->value[0] = v1;
    stencil->value[1] = v2;
    stencil->value[2] = v3;
    stencil->value[3] = v4;
    return ;
}

void print_line(Stencil *stencil) {
    int i;
    for(i = 0; i < 4; i++) {
        printf("%lf ", stencil->value[i]);
    }
    printf("\n");
    return;
}


