all:test

test:main.cpp libBoltzmann.a libnrutil.a
	@g++ -o test main.cpp -L. -lBoltzmann -lnrutil

libBoltzmann.a:Boltzmann.o
	@ar rcs libBoltzmann.a Boltzmann.o

Boltzmann.o:Boltzmann.cpp
	@g++ -c Boltzmann.cpp

libnrutil.a:nrutil.o Integ.o
	@ar rcs libnrutil.a nrutil.o Integ.o

nrutil.o:nrutil.cpp
	@g++ -c nrutil.cpp

Integ.o:Integ.cpp
	@g++ -c Integ.cpp

clean:
	@rm -rf *.o *.a test