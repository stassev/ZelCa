CPP=g++
#CFLAGS= -O0 -ggdb3
CFLAGS=   -m64 -O3 -Wall -Wextra -pipe -march=native -mtune=native -DEIGEN_NO_DEBUG -fPIC

#LIBS= -L../DEPENDENCIES/Cuba-3.2/lib 
#INCLUDES= -I../DEPENDENCIES/eigen-3.2.0/include/eigen3 -I../DEPENDENCIES/Cuba-3.2/include/

LIBS+= -lm -lgsl -lgslcblas -lcuba 

ZStaC :  auxGammaChi.o  auxNij.o auxXi.o auxCosmo.o auxXiRS.o aux3pt.o auxXiRS2d.o
		${CPP} ${CFLAGS} -o ZelCa main.c  auxGammaChi.o auxNij.o auxXi.o auxCosmo.o auxXiRS.o auxXiRS2d.o aux3pt.o ${LIBS} 

auxGammaChi.o :
		${CPP} ${CFLAGS} -c auxGammaChi.c  ${INCLUDES} 

auxNij.o :
		${CPP} ${CFLAGS} -c auxNij.c   ${INCLUDES} 

auxXi.o :
		${CPP} ${CFLAGS} -c auxXi.c  ${INCLUDES} 

auxXiRS.o :
		${CPP} ${CFLAGS} -c auxXiRS.c    ${INCLUDES} 


auxXiRS2d.o :
		${CPP} ${CFLAGS} -c auxXiRS2d.c    ${INCLUDES} 
		
aux3pt.o :
		${CPP} ${CFLAGS} -c aux3pt.c    ${INCLUDES} 

auxCosmo.o :
		${CPP} ${CFLAGS} -c auxCosmo.c 

clean:
		rm -f *.o
