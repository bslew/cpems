#!/bin/sh
# based on Ubuntu-18.04

# CMAKE
apt install cmake  cmake-data:all
# open MPI
apt install libopenmpi2 openmpi-common:all openmpi-bin mpi-default-bin libopenmpi-dev 
# algebra and math
apt install libblas3 liblapack3 libarmadillo libarmadillo-dev libgslcblas0 libblas-dev libgsl23 libgsl-dev libgmp-dev
# data files
apt install libhdf5-100 libhdf5-cpp-100 libhdf5-dev libhdf5-openmpi-100 libhdf5-openmpi-dev  libjsoncpp1 
# graphics
apt install libpng-dev libpng-tools libjpeg-dev 
# plotting (to become depreciated)
apt install libcpgplot0 libpgplot0 pgplot5  
# projections
apt install proj-bin libproj-dev proj-data:all libproj12 
# fftw
apt install fftw2 libfftw3-long3 libfftw3-quad3 libfftw3-bin libfftw3-dev libfftw3-mpi3 libfftw3-mpi-dev
# interpolation and other stuff
apt install libcgal-dev libcgal-qt5-13
# qt5 
apt install qt5-default 
# fits
apt install libcfitsio-dev libccfits-dev libcfitsio-bin
# astro stuff
apt install libnova-0.16-0 libnova-dev 
# system
apt install gfortran libx11-dev
# parser
apt install libtclap-dev
