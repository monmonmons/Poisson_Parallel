# c语言的编译器
CC = mpicc

# -g 增加调试信息
# -Wall 打开大部分警告信息
CFLAGS = -g -Wall -O2

# program依赖哪些目标文件
MAINOBJS = main.o CG.o BiCGSTAB.o Jacobi.o poisson.o Precond.o stencil.o BLAS.o


program : $(MAINOBJS)
	$(CC) $(CFLAGS) -o program $(MAINOBJS) -lm

main.o : main.c poisson.h stencil.h
	$(CC) $(CFLAGS) -c -o main.o main.c

Jacobi.o : Jacobi.c Jacobi.h stencil.h
	$(CC) $(CFLAGS) -c -o Jacobi.o Jacobi.c

BiCGSTAB.o: BiCGSTAB.c BiCGSTAB.h Precond.h
	$(CC) $(CFLAGS) -c -o BiCGSTAB.o BiCGSTAB.c

CG.o : CG.c CG.h stencil.h
	$(CC) $(CFLAGS) -c -o CG.o CG.c

Precond.o :Precond.c Precond.h stencil.h
	$(CC) $(CFLAGS) -c -o Precond.o Precond.c
	
poisson.o : poisson.c poisson.h stencil.h
	$(CC) $(CFLAGS) -c -o poisson.o poisson.c

stencil.o : stencil.c stencil.h BLAS.h
	$(CC) $(CFLAGS) -c -o stencil.o stencil.c

BLAS.o 	  : BLAS.c BLAS.h
	$(CC) $(CFLAGS) -c -o BLAS.o BLAS.c

.PHONY : clean

clean :
	rm -f $(MAINOBJS) program