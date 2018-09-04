GCC=c++ -fPIC
all: lib eph

targets: novas eph_manager novascon nutation solsys1 readeph0

novas:
	${GCC} -c $@.c
eph_manager:
	${GCC} -c $@.c 

novascon: 
	${GCC} -c $@.c

nutation:
	${GCC} -c $@.c

solsys1: 
	${GCC} -c $@.c

readeph0:
	${GCC} -c $@.c

example:
	${GCC}  $@.c -o $@ -lnovas -L.

lib: targets
	ar ru libnovas.a novas.o eph_manager.o novascon.o nutation.o solsys1.o readeph0.o
	ranlib libnovas.a

clean:
	rm -f *.o
	
tidy: clean
	-rm -f JPLEPH
	-rm asc2eph example
	-rm de421/JPLEPH
	-rm de424/JPLEPH

asc2eph:
	gfortran $@.f -o $@

eph: asc2eph
	-rm JPLEPH
	cat de421/header.421 de421/ascp1900.421 | ./asc2eph
	mv JPLEPH de421
	ln -s de421/JPLEPH JPLEPH
