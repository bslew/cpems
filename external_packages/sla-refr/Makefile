
all: lib

targets: 
	gfortran -fPIC -c *.f


lib: targets
	ar ru libslaRefr.a *.o
	ranlib libslaRefr.a

clean:
	rm -f *.o
tidy: clean
	rm *.a
