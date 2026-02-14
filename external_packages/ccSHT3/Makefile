#  Define the C compiler for serial code:
CC = c++

#  Define the C compiler for parallel code:
MPCC = mpicc

#  Define the directory which contains the 
#  FFTW header files:
#FFTWINCLUDE = /usr/include
#FFTWINCLUDE = /home/blew/programy/ccSHTlib/fftw-2.1.5/fftw
#  Define the directory which contains the 
#  FFTW library files:
#FFTWLIBRARY = /usr/lib
#FFTWLIBRARY = /home/blew/programy/ccSHTlib/fftw-2.1.5/fftw-inst
#  If FORTRAN adds an underscore to the object names
#  or uses upper case object names one or both of 
#  the following should be added to the CFLAGS and 
#  MPCFLAGS definitions: 
#
#  -DF77_APPEND_UNDERSCORE
#  -DF77_UPPERCASE_NAMES 

#  Define the flags to be passed to the 
#  serial C complier.  
CFLAGS = -O3 #-Wall
#CFLAGS+=-DFFTW3
#  Define the flags to be passed to the 
#  parallel C compiler.  
MPCFLAGS = -O3

INSTALL_PATH=/usr/local/lib

##################################################
#    DO NOT CHANGE MAKEFILE BELOW THIS POINT     #
##################################################

#.PHONY : all serial library executable clean tidy

all : libccSHT3.a install
#ccSHTx libccSHTmpi.a ccSHTmpix tidy


clean : 
	rm -f lib/libccSHT3.a 

tidy :
	rm -f *.o 

install:
#	sudo install -m 755 -d libccSHT3.a $(INSTALL_PATH)
	sudo install -m 644 -D include/ccSHT3.h $(INSTALL_PATH)/../include/ccSHT3/ccSHT3.h
	sudo install -m 644 -D include/fftw_complex_helper.h $(INSTALL_PATH)/../include/ccSHT3/fftw_complex_helper.h
	sudo install -m 644 -D include/generalTools.h $(INSTALL_PATH)/../include/ccSHT3/generalTools.h
	sudo install -m 644 -D include/sizedType.h $(INSTALL_PATH)/../include/ccSHT3/sizedType.h
	sudo cp libccSHT3.a $(INSTALL_PATH)
#	sudo cp include/ccSHT3.h $(INSTALL_PATH)/../include/ccSHT3
#	sudo cp include/fftw_complex_helper.h $(INSTALL_PATH)/../include/ccSHT3/
#	sudo cp include/generalTools.h $(INSTALL_PATH)/../include/ccSHT3/
#	sudo cp include/sizedType.h $(INSTALL_PATH)/../include/ccSHT3/

FFTW = -I$(FFTWINCLUDE) -L$(FFTWLIBRARY)
CFLAGS2 = $(CFLAGS) -fPIC
#MPCFLAGS2 = $(MPCFLAGS) -D USE_MPI
INCLUDE = -I./include $(FFTW)
LINK = -L./lib
VPATH = src:include:lib


generalTools_ser.o : generalTools.cpp generalTools.h sizedType.h
	$(CC) $(CFLAGS2) $(INCLUDE) -c src/generalTools.cpp 
	mv generalTools.o generalTools_ser.o

fftw_complex_helper.o : fftw_complex_helper.cpp fftw_complex_helper.h generalTools.h sizedType.h
	$(CC) $(CFLAGS2) $(INCLUDE) $(FFTW) -c src/fftw_complex_helper.cpp 

ccSHT3.o : ccSHT3.cpp ccSHT3.h generalTools.h fftw_complex_helper.h sizedType.h
	$(CC) $(CFLAGS2) $(INCLUDE) $(FFTW) -c src/ccSHT3.cpp


#ccSHTfortran.o : ccSHTfortran.c ccSHT3.h generalTools.h fftw_complex_helper.h sizedType.h ccSHTfortran.h
#	$(CC) $(CFLAGS2) $(INCLUDE) $(FFTW) -c src/ccSHTfortran.c


libccSHT3.a : generalTools_ser.o fftw_complex_helper.o ccSHT3.o 
	ar rc libccSHT3.a generalTools_ser.o fftw_complex_helper.o ccSHT3.o
	ranlib libccSHT3.a

