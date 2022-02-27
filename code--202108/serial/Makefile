# c语言的编译器
CC = gcc

# -g 增加调试信息
# -Wall 打开大部分警告信息
CFLAGS = -g -Wall -O3

# program依赖哪些目标文件
MAINOBJS = main.o CG.o BiCGSTAB.o Jacobi.o poisson.o Precond.o ops.o stencil.o


program   : $(MAINOBJS)
	$(CC) $(CFLAGS) -o    program $(MAINOBJS) -lm

main.o    : main.c poisson.h ops.h
	$(CC) $(CFLAGS) -c -o main.o main.c

Jacobi.o  : Jacobi.c Jacobi.h ops.h
	$(CC) $(CFLAGS) -c -o Jacobi.o Jacobi.c

BiCGSTAB.o: BiCGSTAB.c BiCGSTAB.h Precond.h
	$(CC) $(CFLAGS) -c -o BiCGSTAB.o BiCGSTAB.c

CG.o 	  : CG.c CG.h Precond.h
	$(CC) $(CFLAGS) -c -o CG.o CG.c

Precond.o :Precond.c Precond.h ops.h
	$(CC) $(CFLAGS) -c -o Precond.o Precond.c

poisson.o : poisson.c poisson.h ops.h
	$(CC) $(CFLAGS) -c -o poisson.o poisson.c

ops.o 	  : ops.c ops.h stencil.h
	$(CC) $(CFLAGS) -c -o ops.o ops.c

stencil.o : stencil.c stencil.h 
	$(CC) $(CFLAGS) -c -o stencil.o stencil.c

.PHONY : clean

clean :
	rm -f $(MAINOBJS) program