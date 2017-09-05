
all: lib eph

targets: novas eph_manager novascon nutation solsys1 readeph0

novas:
	c++ -c $@.c
eph_manager:
	c++ -c $@.c 

novascon: 
	c++ -c $@.c

nutation:
	c++ -c $@.c

solsys1: 
	c++ -c $@.c

readeph0:
	c++ -c $@.c

example:
	c++  $@.c -o $@ -lnovas -L.

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
