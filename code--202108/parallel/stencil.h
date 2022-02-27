/*******************************************************
    > File   : stencil.c
    > Author : Yuntong
    > Mail   : 171840067@smail.nju.edu.cn
    > Date   : 2021/8/20
*******************************************************/

#ifndef stencil_h
#define stencil_h

//Data structure
typedef struct stencil {
    double value[4]; //0-west 1-east 2-north 3-south
} Stencil;

//Initialize
//  Input : stencil, 4 values
//  Modified : stencil
void initial_stencil(Stencil *stencil, double v1, double v2, double v3, double v4);


void print_line(Stencil *stencil);


#endif /* stencil_h */
