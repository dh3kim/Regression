# makefile.

# user-defined variables.
F95 = f95
OPT1 = -c
OPT2 = -o
OPT3 = -free -V -o
MODULE = -M/home/dokim/lib/lapack95_modules
LAPACK95 = /home/dokim/lib/lapack95.a
LAPACK = /home/dokim/lib/lapack.a
BLAS = /home/dokim/lib/blas.a
TMG = /home/dokim/lib/tmglib.a
AR = ar cr
#LIB = ./tmp_lib.a

# main part of makefile.
all: Malawi_OLS ReadMalawiShp Malawi_la ReadColumbusShp

Malawi_OLS: LU.o OLS_LU.o Malawi.o
	${F95} ${OPT2} Malawi_OLS Malawi.o OLS_LU.o LU.o

Malawi_la: Malawi_test1.o OLS_LA.o   
	${F95} ${OPT2} Malawi_la Malawi_test1.o OLS_LA.o \
        ${LAPACK95} ${LAPACK} ${BLAS} 

OLS_LA.o: OLS_LA.f95
	${F95} ${OPT1} ${MODULE} OLS_LA.f95

LU.o: LU.f95 
	${F95} ${OPT1} LU.f95 

OLS_LU.o: OLS_LU.f95 
	${F95} ${OPT1} OLS_LU.f95

Malawi_test1.o: Malawi_test1.f95
	${F95} ${OPT1} Malawi_test1.f95

Malawi.o: Malawi.f95
	${F95} ${OPT1} Malawi.f95 

ReadMalawiShp: byte_order.o ReadSHP.o ReadMalawiShp.o
	${F95} ${OPT2} ReadMalawiShp ReadMalawiShp.o ReadSHP.o byte_order.o

ReadColumbusShp.o: ReadColumbusShp.f95
	${F95} ${OPT1} ReadColumbusShp.f95

ReadColumbusShp: byte_order.o ReadSHP.o ReadColumbusShp.o
	${F95} ${OPT2} ReadColumbusShp ReadColumbusShp.o ReadSHP.o byte_order.o

byte_order.o: byte_order.f95
	${F95} ${OPT1} byte_order.f95

ReadSHP.o: ReadSHP.f95
	${F95} ${OPT1} ReadSHP.f95

ReadMalawiShp.o: ReadMalawiShp.f95
	${F95} ${OPT1} ReadMalawiShp.f95

clear:
	rm -f *.o

