
CC = g++

CFLAGS = -g3 -O3 -DHAVE_INLINE 

all: lai_regress 

lai_regress: lai_regress.o linear.o
	${CC} ${CFLAGS} lai_regress.o linear.o -lz -Llib/ -lopenblas -llapack -o lai_regress

lai_regress.o: lai_regress.cpp lai_regress.h 
	${CC} ${CFLAGS} -c lai_regress.cpp

linear.o: lai_regress.cpp lai_regress.h linear.cpp linear.h t_distribution.hpp
	${CC} ${CFLAGS} -c -llapack linear.cpp

clean:
	rm -f *.o
