
all: lib

targets: 
	gfortran -fPIC -c *.f


lib: targets
	ar ru libvelKB.a *.o
	ranlib libvelKB.a

clean:
	rm -f *.o
tidy: clean
	rm *.a
