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

#
# below stuff updated on ubuntu 20.04
#
# logger
apt install libspdlog1/focal libspdlog-dev/focal
# gmp
apt install libgmp-dev/focal libgmp10/focal
# CGAL
apt install libcgal-dev/focal libcgal-qt5-dev/focal libcgal-ipelets/focal
# math
apt install libarmadillo-dev/focal libarmadillo9/focal libproj-dev/focal libproj15/focal libgsl23/focal libgsl-dev/focal  
# gsl
apt install libgsl23/focal libgsl-dev/focal
#proj
apt install libproj-dev/focal libproj15/focal

# HDF5
apt install libhdf5-103 libhdf5-cpp-103 libhdf5-dev libhdf5-openmpi-103
#cmake 
apt install extra-cmake-modules/focal cmake-extras/focal cmake-data/focal

#boost
apt install libboost-test1.71.0/focal libboost-test1.71-dev/focal libboost-filesystem1.71-dev/focal libboost-filesystem1.71.0/focal

#QT5
sudo apt install build-essential qtcreator qt5-default

# openCV
apt install libopencv-dev/focal python3-pyqt5/focal

# yaml
apt install libyaml-cpp-dev/focal libyaml-dev/focal libyaml-cpp0.6/focal  exiv2/focal-security libexiv2-dev/focal-security libgexiv2-2/focal libgexiv2-dev/focal

apt install libboost-filesystem1.71-dev/focal libboost-filesystem1.71.0/focal libboost-test1.71.0/focal libboost-test1.71-dev/focal

#cmake
apt install cmake-fedora/focal  make-data/focal cmake-extras/focal extra-cmake-modules/focal









