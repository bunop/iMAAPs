all: iMAAPs

iMAAPs: data.o err.o fitting.o par.o tools.o help.o
	g++ data.o err.o fitting.o par.o tools.o help.o iMAAPs.cpp -o iMAAPs  -fopenmp  -lm  -lgsl -lgslcblas -lfftw3

data.o: data.cpp par.o fitting.o par.o tools.o
	g++ -c data.cpp -fopenmp -lm -lpthread -lgsl -lgslcblas -lfftw3
	
err.o: err.cpp
	g++ -c err.cpp
	
fitting.o: fitting.cpp
	g++ -c fitting.cpp -lgsl -lgslcblas -lfftw3

par.o: par.cpp
	g++ -c par.cpp

tools.o: tools.cpp
	g++ -c tools.cpp

help.o: help.cpp
	g++ -c help.cpp

clean:
	rm -rf *.o iMAAPs

.PHONY: all clean
