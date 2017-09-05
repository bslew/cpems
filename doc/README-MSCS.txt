 * INFORMATION:

This is the MscsCPEDS combined package dedicated to Cosmic Microwave
Background (CMB) map simulations and statistics computations.


 * LICENCE

This package will be availale under GNU general public licence, 
however some of the libraries it uses have similar licences of use,
however not entirely consistent withg GNU. 


 * REQUIREMENTS

For proper compilation you need a number of libraries (mostly under GNU
licence)

 - CPEDS ( included in this distribution )
 - (c)pgplot
 - proj4 
 - fftw3
 - matpack
 - matplotlib (with basemap toolkit)
 - python
 - gsl
 - gslcblas
 - cfitsio
 - X11, Xpm, g2c, z, 
 - standard c/c++ libraries and header files (so devels generally also might
   be needed)


If you have them installed in non-standard places you have to modify 
the Makefile accordingly. See comment strings in the Makefile file.


 * COMPILATION and INSTALATION

If there is no problem with finding the libraries and header files
all you need to do is the type:

make all  

For making object files only:

make parts

For making object files and general programs:

make parts progs


The binaries are installed in the ./bin directory
so you have to include this directory in the your PATH
variable to make your shell look for them there.

 * DOCUMENTATION

For programming with Mscs, see comments in the header and source files 
For programs usage type option -h or --help for some info.

At the moment there is no manual for this package. It's still strongly
under developent

If you have questions feel free to write me at: blew at th.nao.ac.jp


 * BUGS

For any bug reports, please write authors.

 * TROUBLESHOOTING

(1)
if CPEDS package do not compile because of some implicit function
declarations like round() or trunc(), try add 
#define _ISOC99_SOURCE
#include <features.h>
at the top of cpeds-math.c and cpeds-pixelization.c files

(2)
if making some executables fails due to some undefined references to some
procedures from ccSHT or fftw package make sure that you use fftw2.1.5 and
that these were compiled using c++ compiler. Some modification to thiese
packages might be needed.


(3) mappack compilation:
* Sometimes there seem to be a problem in compilation of matpack giving errorr message of missing strtol functions etc.
The easiest remedy for that is to include the <stdlib.h> in the beginning of that file
* Also if you get messages like: mpbutton.cpp:34:21: error: X11/xpm.h: No such file or directory
you most likely need package
libxpm-dev on Ubuntu systems
and 
libXpm-devel-3.5.7-1.fc8 on fedora etc.

(4) lam compilation:
Depending on the version of your c++ compiler you might encounter the error saying that your compiler does not support 
the bool type. In such cases it might be that the version of the compiler checks the lam more strictly than it supposed to.
In such cases passing option CXXFLAGS=-fpermissive to the configure script might help.





 * AUTHORS

Bartosz Lew <blew at th.nao.ac.jp>


