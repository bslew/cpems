cmake_minimum_required (VERSION 2.6)

################################################################
################################################################
################################################################
# Set project name
#
project (cpems CXX C)
#enable_language(<lang> [OPTIONAL] )

################################################################
################################################################
################################################################
# Set version number
#
set (CPEMS_VERSION "0.2")

# Set build type

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

################################################################
################################################################
################################################################
# configure a header file to pass some of the CMake settings
# to the source code
#configure_file (
#  "${PROJECT_SOURCE_DIR}/src/PGEEJ_config.h.in"
#  "${PROJECT_BINARY_DIR}/PGEEJ_config.h"
#  )

set (CMAKE_MODULE_PATH 
#	${CMAKE_ROOT}/Modules
	${PROJECT_SOURCE_DIR}/cmake
	${CMAKE_MODULE_PATH}
	)


################################################################
################################################################
################################################################
# Build external packages first
#
#ADD_CUSTOM_COMMAND(
#        TARGET CPEDS
#        PRE_BUILD
#        WORKING_DIRECTORY external_packages/velKB
#        OUTPUT external_packages/velKB
#        COMMAND env CXXFLAGS=${ROOTS_COMPILER_FLAGS} ./configure ${TARGET_FLAG} ${ROOTS_HOST_OPTION}
#        COMMAND make
#        COMMENT "Building velKB"
#        VERBATIM
#    )

#set (VELKB_LIB_DIR "external_packages/velKB") 

add_custom_target(extdir
    COMMAND /bin/bash -c 'if [ ! -d ${PROJECT_BINARY_DIR}/external_packages ]\; then mkdir ${PROJECT_BINARY_DIR}/external_packages\; fi'
    )
    
add_custom_target(Lib_velKB
    COMMAND cp -rp ${PROJECT_SOURCE_DIR}/external_packages/velKB ${PROJECT_BINARY_DIR}/external_packages
    COMMAND cd ${PROJECT_BINARY_DIR}/external_packages/velKB && make
    COMMENT "Building velKB"
    VERBATIM
    )

add_custom_target(Lib_slaRefr 
    COMMAND cp -rp ${PROJECT_SOURCE_DIR}/external_packages/sla-refr ${PROJECT_BINARY_DIR}/external_packages
    COMMAND cd ${PROJECT_BINARY_DIR}/external_packages/sla-refr && make
    COMMENT "Building sla-refr"
    VERBATIM
    )

add_custom_target(Lib_novas
    COMMAND cp -rp ${PROJECT_SOURCE_DIR}/external_packages/novas ${PROJECT_BINARY_DIR}/external_packages
    COMMAND cd ${PROJECT_BINARY_DIR}/external_packages/novas && make
    COMMENT "Building novas"
    VERBATIM
    )

add_custom_target(Lib_ccSHT3
    COMMAND cp -rp ${PROJECT_SOURCE_DIR}/external_packages/ccSHT3 ${PROJECT_BINARY_DIR}/external_packages
    COMMAND cd ${PROJECT_BINARY_DIR}/external_packages/ccSHT3 && make
    COMMENT "Building ccSHT3"
    VERBATIM
    )

add_custom_target(external
    COMMENT "Building external_packages"
    )

add_dependencies(Lib_velKB extdir)
add_dependencies(Lib_slaRefr extdir)
add_dependencies(Lib_novas extdir)
add_dependencies(Lib_ccSHT3 extdir)
add_dependencies(external extdir Lib_velKB Lib_slaRefr Lib_novas Lib_ccSHT3)

#add_custom_target(gsl18
#    COMMAND cd ${PROJECT_BINARY_DIR}/external_packages/gsl-1.8 && ./configure && make -j10
#    COMMAND 
#    COMMENT "Building gsl18"
#    VERBATIM
#    )

#add_custom_command(
#    TARGET velKB PRE_BUILD 
#    WORKING_DIRECTORY external_packages/velKB
#    COMMAND make
#    COMMENT "Building velKB"
#    VERBATIM
#    )

################################################################
################################################################
################################################################
# Find required packages
#

find_package(CFITSIO REQUIRED)
#  CFITSIO_INCLUDE_DIR - the CFITSIO include directory
#  CFITSIO_LIBRARIES - Link these to use CFITSIO
if (CFITSIO_FOUND)
	MESSAGE(STATUS "CFITSIO_INCLUDE_DIR: " ${CFITSIO_INCLUDE_DIR})
	MESSAGE(STATUS "CFITSIO_LIBRARIES: " ${CFITSIO_LIBRARIES})
endif()

FIND_PACKAGE(GSL 2 REQUIRED)
# GSL_INCLUDE_DIR      = where to find headers 
# GSL_LIBRARIES   = full path to the libraries
if (GSL_FOUND)
	MESSAGE(STATUS "GSL_INCLUDE_DIR: " ${GSL_INCLUDE_DIR})
	MESSAGE(STATUS "GSL_INCLUDES: " ${GSL_INCLUDES})
	MESSAGE(STATUS "GSL_LIBRARIES: " ${GSL_LIBRARIES})
endif()


find_package(HDF5 1.8 REQUIRED)
# defines HDF5_INCLUDE_DIRS HDF5_LIBRARIES
if (HDF5_FOUND)
	MESSAGE(STATUS "HDF5_INCLUDE_DIRS: " ${HDF5_INCLUDE_DIRS})
	MESSAGE(STATUS "HDF5_LIBRARIES: " ${HDF5_LIBRARIES})
endif()

#find_package(Qt5 COMPONENTS OpenGL Gui Core Script ScriptTools QUIET)
find_package(Qt5 COMPONENTS Core REQUIRED)
# the recommended use of Qt5 is different than other packages found see:
# https://doc.qt.io/qt-5/cmake-manual.html

find_package(NOVA REQUIRED)
#  NOVA_INCLUDE_DIR - the NOVA include directory
#  NOVA_LIBRARIES - Link these to use NOVA

#find_package(m REQUIRED)

find_package(PROJ REQUIRED)
#    PROJ_INCLUDE_DIR
#    PROJ_LIBRARY

find_package(FFTW 3 REQUIRED)
#   FFTW_LIBRARIES           ... full path to fftw library
#   FFTW_INCLUDES            ... fftw include directory

if (FFTW_FOUND)
	MESSAGE(STATUS "FFTW_INCLUDE_DIR: " ${FFTW_INCLUDE_DIR})
	MESSAGE(STATUS "FFTW_INCLUDES: " ${FFTW_INCLUDES})
	MESSAGE(STATUS "FFTW_LIBRARIES: " ${FFTW_LIBRARIES})
endif()


find_package(GMP REQUIRED)
# GMP_INCLUDES
# GMP_LIBRARIES

find_package(CGAL REQUIRED)
# CGAL_INCLUDE_DIRS
# CGAL_LIBRARIES_DIR
# CGAL_MODULES_DIR
if (CGAL_FOUND)
	MESSAGE(STATUS "CGAL_INCLUDE_DIRS: " ${CGAL_INCLUDE_DIRS})
	MESSAGE(STATUS "CGAL_LIBRARIES_DIR: " ${CGAL_LIBRARIES_DIR})
	MESSAGE(STATUS "CGAL_MODULES_DIR: " ${CGAL_MODULES_DIR})
endif()

find_package(PGPLOT REQUIRED)
#  PGPLOT_INCLUDES    = Path to the PGPLOT header files
#  PGPLOT_LIBRARIES   = Path to all parts of the PGPLOT library
#  PGPLOT_LIBRARY_DIR = Path to the directory containing the PGPLOT libraries


find_package(Armadillo REQUIRED)
# ARMADILLO_INCLUDE_DIRS - list of required include directories
# ARMADILLO_LIBRARIES - list of libraries to be linked

################################################################
################################################################
################################################################
# Find required headers and libraries
#
#find_path(CPEDS_INCLUDE_DIR cpeds-msgs.h 
#	$ENV{CPEDS_DIR}/include 
#	/usr/local/lib
#	)

if (NOT CFITSIO_FOUND)
	find_path(FITSIO_INCLUDE_DIR fitsio.h 
		/usr/include/cfitsio 
		/usr/local/include/cfitsio 
		$ENV{HOME}/local/include
		)
endif()

if (NOT GSL_FOUND)
	find_path(GSL_INCLUDE_DIR gsl_multifit.h
		$ENV{HOME}/local/include/gsl
		/usr/local/include/gsl 
		/usr/include/gsl
		)

	find_library (GSL_LIBRARIES gsl gslcblas m)
		
endif()


# nova library
if (NOT NOVA_FOUND)
	find_library (NOVA_LIBRARIES nova
		PATHS
		$ENV{HOME}/local
		$ENV{HOME}/anaconda3
		PATH_SUFFIXES lib lib64
	)
endif()

# proj library
if (NOT PROJ_FOUND)
	find_library (PROJ_LIB proj
		PATHS
		$ENV{HOME}/local
		$ENV{HOME}/anaconda3
		PATH_SUFFIXES lib lib64
	)
endif()

# QtCore library
if (NOT Qt5Core_FOUND)
	find_library (QTCORE_LIBRARY Qt5Core
		PATHS
		$ENV{HOME}/local
		$ENV{HOME}/anaconda3
		PATH_SUFFIXES lib lib64
	)
endif()


# FFTW library
if (NOT FFTW_FOUND)
	find_library (FFTW_LIBRARIES fftw3 fftw3l
		PATHS
		$ENV{HOME}/local
		$ENV{HOME}/anaconda3
		PATH_SUFFIXES lib lib64
	)
#	find_library (FFTWL3_LIB fftw3l
#		PATHS
#		$ENV{HOME}/local/lib
#		$ENV{HOME}/anaconda3/lib
#	)
endif()


# HDF5 library
if (NOT HDF5_FOUND)
	find_path(HDF_INCLUDE_DIR hdf5.h
		$ENV{HOME}/anaconda3/include
		$ENV{HOME}/local/include
		/usr/local/include/gsl 
		/usr/include/gsl
		)

	find_library (HDF5_LIBRARIES hdf5 z dl m
		PATHS
		$ENV{HOME}/local
		$ENV{HOME}/anaconda3
		PATH_SUFFIXES lib lib64
	)

endif()


# CGAL library
if (NOT CGAL_FOUND)
	find_library (CGAL_LIBRARIES_DIR CGAL
		PATHS
		$ENV{HOME}/local
		PATH_SUFFIXES lib lib64
	)
endif()


# GMP library
if (NOT GMP_FOUND)
	find_library (GMP_LIBRARIES gmp
		PATHS
		/usr
		/local
		$ENV{HOME}/local
		$ENV{HOME}/anaconda3
		PATH_SUFFIXES lib lib64
	)
endif()


# fitsio library
if (NOT FITSIO_FOUND)
	find_library (FITSIO_LIB cfitsio
		PATHS
		$ENV{HOME}/local
		PATH_SUFFIXES lib lib64
	)
endif()


# PGPLOT library
if (NOT PGPLOT_FOUND)
	find_library (PGPLOT_LIBRARIES pgplot cpgplot 
		PATHS
		$ENV{HOME}/local
		PATH_SUFFIXES lib lib64
	)
endif()


# armadillo library
if (NOT Armadillo_FOUND)
	find_library (ARMADILLO_LIB armadillo
		PATHS
		$ENV{HOME}/local
		PATH_SUFFIXES lib lib64
	)
endif()

#find_path(QTCORE_HEADERS QList
#	$ENV{HOME}/anaconda3/include/qt
#	$ENV{HOME}/anaconda3/include/qt/QtCore
#    /usr/include/qt5/QtCore
#)

find_path (STL_INCLUDE_DIR NAMES 
	stl_iterator.h 
#	c++/7/bits/ 
#	c++/8/bits/stl_iterator.h 
	PATHS
	/
	ENV{GCC_COMPILER_DEFAULTS_PATH}/include/c++/8.2.0/bits
	gcc/8.2.0
#	/opt/exp_soft/local/generic/gcc/8.2.0/include/c++/8.2.0
#	/opt/exp_soft/local/generic/gcc/8.2.0/include/c++/8.2.0/bits	
#	*/c++/7/bits/stl_iterator.h 
	/usr/include/c++/7/bits 
	PATH_SUFFIXES bits
	)

#execute_process(
#  COMMAND "which c++"
#  OUTPUT_VARIABLE PYTHON_SITE
#  OUTPUT_STRIP_TRAILING_WHITESPACE)

################################################################
################################################################
# find python installation path
execute_process(
  COMMAND "python" -c "from distutils import sysconfig as sc
print(sc.get_python_lib(prefix='', plat_specific=True))"
  OUTPUT_VARIABLE PYTHON_SITE
  OUTPUT_STRIP_TRAILING_WHITESPACE)

################################################################
################################################################
################################################################
# Output status
#
# MESSAGE( STATUS "CPEDS_INCLUDE_DIR: " ${CPEDS_INCLUDE_DIR} )
#MESSAGE( STATUS "FITSIO_INCLUDE_DIR: " ${FITSIO_INCLUDE_DIR} )
#MESSAGE( STATUS "PROJECT_BINARY_DIR: " ${PROJECT_BINARY_DIR})

MESSAGE( STATUS "PYTHON_SITE: " ${PYTHON_SITE})
#MESSAGE( STATUS "PYTHON: " ${PYTHON_EXECUTABLE})
################################################################
################################################################
################################################################
# Include project directories to search for headers
#
include_directories(
#	"${PROJECT_BINARY_DIR}/src"
    src
	external_packages/novas
	external_packages/ccSHT3/include
	${CFITSIO_INCLUDE_DIR}
	${GSL_INCLUDE_DIR}
#	${QTCORE_HEADERS}
#	${QTCORE_HEADERS}/..
	${HDF_INCLUDE_DIRS}
#	${STL_INCLUDE_DIR}
#	get_filename_component(PARENT_DIR ${QTCORE_HEADERS} DIRECTORY)
#	${CPEDS_INCLUDE_DIR}
#	${MSCS_INCLUDE_DIR}
    /usr/include/python2.7
	)

#add_subdirectory(external_packages/novas)
#add_subdirectory(external_packages/velKB)

################################################################
################################################################
################################################################
# Include project directories to search for libraries
#
link_directories(
    ${PROJECT_BINARY_DIR}/external_packages/velKB
    ${PROJECT_BINARY_DIR}/external_packages/sla-refr
    ${PROJECT_BINARY_DIR}/external_packages/ccSHT3/lib
    ${PROJECT_BINARY_DIR}/external_packages/novas
    )


################################################################
################################################################
################################################################
# Define source files
#
set (CPEDSLIB_HEADERS 
    src/qtbase.h
    src/cpeds-common.h
	src/cpeds-consts.h
	src/cpeds-rng.h
	src/cpeds-math.h 
	src/cpeds-cosmo.h 
	src/cpeds-pdf.h 
	src/cpedsTime.h
	src/cpeds-msgs.h
	src/cpeds-point_set.h
	src/cpeds-point2d.h
	src/cpeds-point3d.h
	src/cpeds-direction_set.h
	src/cpedsFunctionDistance.h
	src/cpeds-project.h
	src/cpeds-smooth.h
	src/cpeds-direction.h
	src/cpedsMoonDirection.h
	src/cpedsPlanetDirection.h
	src/cpeds-function-cosmo.h
	src/cpedsMC.h
	src/cpedsMCMC.h
	src/cpedsMCMCset.h
	src/Nbody_io.h
	src/cpeds-fits.h
	src/matrix.h
	src/mscsVector.h
	src/Mscs-object.h
	src/cpeds-list.h
	src/MClink.h
	src/MCparameterSpace.h
    src/cpeds-templates.h
	src/MscsDateTimeData.h
	src/RT32pointing_models.h
	src/cpeds_angular_correlation_fn.h
	)

set (CPEDSLIB_SOURCES 
    ${CPEDSLIB_HEADERS}
	src/cpeds-consts.cpp
	src/cpeds-rng.cpp
	src/cpeds-math.cpp 
	src/cpeds-pixelization.cpp
	src/cpeds-cosmo.cpp 
	src/cpeds-pdf.cpp 
	src/cpedsTime.cpp
	src/cpeds-msgs.cpp
	src/cpeds-point_set.cpp
	src/cpeds-direction_set.cpp
	src/cpedsFunctionDistance.cpp
	src/cpeds-project.cpp
	src/cpeds-smooth.cpp
	src/cpeds-direction.cpp
	src/cpedsMoonDirection.cpp
	src/cpedsPlanetDirection.cpp
	src/cpeds-function-cosmo.cpp
	src/cpedsMCMC.cpp 
	src/cpedsMC.cpp
	src/MClink.cpp
	src/MCparameterSpace.cpp
	src/cpedsMCMCset.cpp
	src/Nbody_io.cpp 
	src/cpeds-fits.cpp 
#    src/cpeds-templates.cpp
	src/MscsDateTimeData.cpp
	src/RT32pointing_models.cpp
	src/cpeds_angular_correlation_fn.cpp
	)
#	${CPEDS_INCLUDE_DIR}/cpeds-msgs.h

set (CPEDS_PYTHON_BINDINGS 
    pyth/pyCPEDScommonFunctions/cpedsRotation.cpp
    )

set (MSCSCORELIB_HEADERS
	src/Mscs-global-defs.h
	src/Mscs-common.h
	src/Mscs-map.h
	src/Mscs-alms.h
	src/Mscs-alm.h
	src/Mscs-gaussian_simulation.h
	src/mscsMapINC.h
	src/Mscsgrid.h
	src/Mscs-minkowski-f.h
	)		


set (MSCSCORELIB_SOURCES
    ${MSCSCORELIB_HEADERS} 
	src/Mscs-global-defs.cpp
	src/Mscs-common.cpp 
	src/Mscs-map.cpp 
	src/Mscs-map-io_methods.cpp 
	src/Mscs-alms.cpp 
	src/Mscs-alm.cpp 
	src/Mscs-map-sht.cpp 
	src/Mscs-map-interpolate.cpp 
	src/Mscs-gaussian_simulation.cpp 
	src/mscsMapINC.cpp 
	src/Mscsgrid.cpp 
	src/Mscs-minkowski-f.cpp
	src/Mscs-map-minkowski-functionals.cpp
	src/Mscs-map-differencial-calculus.cpp 
	)		


set (MSCSFNLIB_HEADERS
	src/mscsFunction3dCLEAN.h
	src/Mscs-function.h
	src/Mscs-map-window_function.h
	src/Mscs-map-beam_function.h
	src/Mscs-correlation_function.h
	src/Mscs-function3dreg.h
	src/Mscs-function3dregc.h
	src/Mscs-function3dregc.h
	src/subdomain.h
	src/mscsSurface3dregc.h
	src/mscsLine.h
	src/Mscs-function2dreg.h
	src/Mscs-angular_power_spectrum.h
	src/mandelbrotSet.h
	src/pointsDensity.h
	src/mscsFunctionFit.h
	src/MscsPDF2D.h
	src/MscsPDF1D.h
	)

set (MSCSFNLIB_SOURCES
    ${MSCSFNLIB_HEADERS}
	src/mscsFunction3dCLEAN.cpp
	src/Mscs-function.cpp 
	src/Mscs-function-hdf5.cpp 
	src/Mscs-map-window_function.cpp 
	src/Mscs-map-beam_function.cpp 
	src/Mscs-correlation_function.cpp 
	src/Mscs-function3dreg.cpp 
	src/Mscs-function3dregc.cpp 
	src/Mscs-function3dregc-hdf5.cpp 
	src/subdomain.cpp 
	src/mscsSurface3dregc.cpp 
	src/mscsLine.cpp 
	src/Mscs-function2dreg.cpp 
	src/Mscs-angular_power_spectrum.cpp 
	src/mandelbrotSet.cpp 
	src/pointsDensity.cpp 
	src/mscsFunctionFit.cpp 
	src/MscsPDF2D.cpp 
	src/MscsPDF1D.cpp 
	)
#	src/Mscs-minkowski-f.c src/Mscs-minkowski-f.h
#	src/Mscs-topo-dodecahedron.c Mscs-topo-dodecahedron.h
#	src/Mscs-map-minkowski-functionals.c src/Mscs-map-minkowski-functionals.h
#	src/Mscs-map-differencial-calculus.c src/Mscs-map-differencial-calculus.h
#	src/Mscs-transfer_function.c src/Mscs-transfer_function.h
#	src/Mscs-matter_power_spectrum.cpp src/Mscs-matter_power_spectrum.h

set (MSCSPLOTLIB_HEADERS
	src/Mscs-map_plot.h
	src/Mscs-colormap.h
	src/Mscs-map_proj.h
	)

set (MSCSPLOTLIB_SOURCES
    ${MSCSPLOTLIB_HEADERS}
	src/Mscs-map_plot.cpp
	src/Mscs-colormap.cpp
	src/Mscs-map_proj.cpp
	)


set (MSCSWMAPLIB_HEADERS
	src/Mscs-WMAPspecifications.h
	src/Mscs-WMAPsimulation.h
	src/MscsWMAPdata.h
	)
	
set (MSCSWMAPLIB_SOURCES
    ${MSCSWMAPLIB_HEADERS}
	src/Mscs-WMAPspecifications.cpp 
	src/Mscs-WMAPsimulation.cpp 
	src/MscsWMAPdata.cpp 
	)
	
################################################################
################################################################
################################################################
# Add library targets
#
#add_library(CPEDS STATIC ${CPEDSLIB_SOURCES} )
#add_library(Mscsfn STATIC ${MSCSFNLIB_SOURCES} )
#add_library(Mscscore STATIC ${MSCSCORELIB_SOURCES} )
#add_library(Mscsplot STATIC ${MSCSPLOTLIB_SOURCES} )
#add_library(MscsWMAP STATIC ${MSCSWMAPLIB_SOURCES} )

add_library(CPEDS SHARED  ${CPEDSLIB_SOURCES} )
#add_library(cpedsRotation SHARED ${CPEDS_PYTHON_BINDINGS} )
add_library(Mscsfn SHARED  ${MSCSFNLIB_SOURCES} )
add_library(Mscscore SHARED ${MSCSCORELIB_SOURCES} )
add_library(Mscsplot SHARED ${MSCSPLOTLIB_SOURCES} )
add_library(MscsWMAP SHARED ${MSCSWMAPLIB_SOURCES} )

#add_library(${PROJECT_BINARY_DIR}/external_packages/
#add_dependencies(CPEDS slaRefr novas ccSHT3 velKB)

################################################################
################################################################
################################################################
# Add CPEDS executable targets
#
add_executable(cosmocalc src/cosmocalc.cpp )
add_executable(integrate src/integrate.cpp )
add_executable(interpolate src/interpolate.cpp )
add_executable(compute_CL_ranges src/compute_CL_ranges.cpp )
add_executable(bin_function src/bin_function.cpp )
add_executable(coco src/coco.cpp )

# there are many more programs in the package and much more tests that also should be compiled somewhere

# Add Mscs executable targets
#
add_executable(select-map-pixels src/select-map-pixels.cpp )
add_executable(matrix_oper src/matrix_oper.cpp )
add_executable(calcCMBflux src/calcCMBflux.cpp )
add_executable(mkCMBfield src/mkCMBfield.cpp)
add_executable(smoothPoints src/smoothPoints.cpp)
add_executable(draw_maps src/draw_maps.cpp)
add_executable(calculate_minkowski_v src/calculate_minkowski_v.cpp )

#draw_mapsMPI	:
#	MPIFLAG = -DGOMPI
#	$(CCMPI) $(OFLAGS) $(MPIFLAG) $(DISTR_DIRcpp)/draw_maps.c -o $(OUTDIR)/draw_mapsMPI $(IFLAGS) $(LLFLAGS) $(LFLAGS) 

add_executable(LMdataFit src/LMdataFit.cpp)
add_executable(function_fft src/function_fft.cpp)
add_executable(mkFractal src/mkFractal.cpp)
add_executable(mk_gauss_map src/mk_gauss_map.cpp)
add_executable(map_oper src/map_oper.cpp)
add_executable(calculate_Cth src/calculate_Cth.cpp)

################################################################
################################################################
################################################################
# Add Mscs executable test targets
#
add_executable(test-mscsfn tests/test-mscsfn.cpp )
add_executable(test-function-finter tests/test-function-finter.cpp )
add_executable(test-mscsfn-hdf5 tests/test-mscsfn-hdf5.cpp )
add_executable(test-domain-decomposition tests/test-domain-decomposition.cpp )
add_executable(test-function-fft tests/test-function-fft.cpp )
add_executable(test-function-derivative tests/test-function-derivative.cpp )
add_executable(test-fits-read tests/test-fits-read.cpp )
add_executable(test-cpedsMCMC-parabola_2param tests/test-cpedsMCMC-parabola_2param.cpp )


add_executable(test-hdf-openMP2 tests/test-hdf-openMP2.cpp)

#test-hdf-openMP2-omp:
#	$(CC) -DENABLE_OMP $(OFLAGS)   $(DISTR_DIRcpp)/test-hdf-openMP2.cpp -o $(OUTDIR)/$@ -lhdf5 -L/usr/local/lib -I/usr/local/include

add_executable(test-hdf-openMP tests/test-hdf-openMP.cpp)
#test-hdf-openMP-omp	:
#	$(CC) -DENABLE_OMP $(OFLAGS)   $(DISTR_DIRcpp)/test-hdf-openMP.cpp -o $(OUTDIR)/$@ $(IFLAGS) $(LLFLAGS) $(LFLAGS) 



#-------------------------------------------
# test programs
# YET NEED TO BE PORTED FROM MSCS
#------------------------------------------------------------------------------------------------------------------------
# WMAP dedicated programs -- these go through all WMAP files by default
# mk_gauss_map_simWMAP	: 
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/mk_gauss_map_simWMAP.cpp -o $(OUTDIR)/mk_gauss_map_simWMAP $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# 
# mk_gauss_map_simWMAP_MPI	: 
# 	$(CCMPI) $(OFLAGS) $(MPIFLAG)   $(DISTR_DIRcpp)/mk_gauss_map_simWMAP.cpp -o $(OUTDIR)/mk_gauss_map_simWMAP_MPI $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# 
# mk_gauss_map_simWMAPexp	: 
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/mk_gauss_map_simWMAP.c -o $(OUTDIR)/mk_gauss_map_simWMAPexp $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# 
# mk_gauss_map_simWMAP2	: 
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/mk_gauss_map_simWMAP.c -o $(OUTDIR)/mk_gauss_map_simWMAP2 $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# 
# mk_INC_WMAP	: 
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/mk_INC_WMAP.c -o $(OUTDIR)/mk_INC_WMAP $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# 
# mk_INC_sim	: 
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/mk_INC_sim.c -o $(OUTDIR)/mk_INC_sim $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# 
# mk_ILC_sim	: 
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/mk_ILC_sim.c -o $(OUTDIR)/mk_ILC_sim $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# 
# map_decompose_WMAP	: 
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/map_decompose_WMAP.c -o $(OUTDIR)/map_decompose_WMAP $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# 
# calculate_alms_WMAP	: 
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/calculate_alms_WMAP.c -o $(OUTDIR)/calculate_alms_WMAP $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# 
# calculate_crossCl_WMAP : 
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/calculate_crossCl_WMAP.c -o $(OUTDIR)/calculate_crossCl_WMAP $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# 
# extract_Cl_WMAP : 
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/extract_Cl_WMAP.c -o $(OUTDIR)/extract_Cl_WMAP $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# 
# recalibrate_alms : 
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/recalibrate_alms.c -o $(OUTDIR)/recalibrate_alms $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# 
# calculate_minkowski_functionals_WMAP : 
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/calculate_minkowski_functionals_WMAP.c -o $(OUTDIR)/calculate_minkowski_functionals_WMAP $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# 
# calculate_minkowski_functionals_sim : 
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/calculate_minkowski_functionals_sim.c -o $(OUTDIR)/calculate_minkowski_functionals_sim $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# #	$(CC) $(OFLAGS) -c  calculate_minkowski_functionals_sim.c $(IFLAGS) 
# #	$(CC) $(OFLAGS)   calculate_minkowski_functionals_sim.o -o $(OUTDIR)calculate_minkowski_functionals_sim $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# 
# 
# plot_minkowski_functionals_sim : 
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/plot_minkowski_functionals_sim.c -o $(OUTDIR)/plot_minkowski_functionals_sim $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# 
# plot_minkowski_v_sim : 
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/plot_minkowski_v_sim.c -o $(OUTDIR)/plot_minkowski_v_sim $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# 
# calculate_chisq_simWMAP : 
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/calculate_chisq_simWMAP.c -o $(OUTDIR)/calculate_chisq_simWMAP $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# 
# calculate_chisq_simWMAP_SK : 
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/calculate_chisq_simWMAP_SK.c -o $(OUTDIR)/calculate_chisq_simWMAP_SK $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# 
# SKregstat : 
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/SKregstat.c -o $(OUTDIR)/SKregstat $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# SKregstat2 : 
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/SKregstat2.c -o $(OUTDIR)/SKregstat2 $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# 
# SKregstat2_s : 
# 	$(CC) -DVARONLY $(OFLAGS)   $(DISTR_DIRcpp)/SKregstat2.c -o $(OUTDIR)/SKregstat2_s $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# 
# SKregstat2_sr : 
# 	$(CC) -DVARRATIOONLY $(OFLAGS)   $(DISTR_DIRcpp)/SKregstat2.c -o $(OUTDIR)/SKregstat2_sr $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# 
# SKregstat2low : 
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/SKregstat2.c -o $(OUTDIR)/SKregstat2-lowres $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# SKregstat_convert : 
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/SKregstat_convert.c -o $(OUTDIR)/SKregstat_convert $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# SKregstat_convert64 : 
# 	$(CC) -DGO64BIT $(OFLAGS)   $(DISTR_DIRcpp)/SKregstat_convert.c -o $(OUTDIR)/SKregstat_convert64 $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# 
# SKregstat_convert_s : 
# 	$(CC) -DVARONLY $(OFLAGS)   $(DISTR_DIRcpp)/SKregstat_convert.c -o $(OUTDIR)/SKregstat_convert_s $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# 
# SKregstat_convert64_s : 
# 	$(CC) -DVARONLY -DGO64BIT $(OFLAGS)   $(DISTR_DIRcpp)/SKregstat_convert.c -o $(OUTDIR)/SKregstat_convert_s $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# 
# SKregstat_cov_stripes : 
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/SKregstat_cov_stripes.c -o $(OUTDIR)/SKregstat_cov_stripes $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# SKregstat--mk_diff_map_quick:
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/SKregstat--mk_diff_map_quick.c -o $(OUTDIR)/SKregstat--mk_diff_map_quick $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# 
# 
# statSKregstat : 
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/statSKregstat.c -o $(OUTDIR)/statSKregstat $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# statSKregstat2 : 
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/statSKregstat.c -o $(OUTDIR)/statSKregstat2 $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# statSKregstat3 : 
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/statSKregstat.c -o $(OUTDIR)/statSKregstat3 $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# statSKregstat--mk_det_maps : 
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/statSKregstat--mk_det_maps.c -o $(OUTDIR)/statSKregstat--mk_det_maps $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# statSKregstat_mk_localNG_map : 
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/statSKregstat_mk_localNG_map.c -o $(OUTDIR)/statSKregstat_mk_localNG_map $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# 
# statSKregstattxt : 
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/statSKregstat.c -o $(OUTDIR)/statSKregstattxt $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# 
# SKregstat_drv : 
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/SKregstat_drv.c -o $(OUTDIR)/SKregstat_drv $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# 
# lagplot : 
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/lagplot.c -o $(OUTDIR)/lagplot $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# 
# mfstatmpi : 
# 	$(CCMPI) $(OFLAGS) $(MPIFLAG)    $(DISTR_DIRcpp)/MFstat.c -o $(OUTDIR)/MFstat $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# 
# 
# 
# #-------------------------------------------
# # help programs
# 
# extract_spikes	: 
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/extract_spikes.cpp -o $(BIN_DIR)/extract_spikes $(IFLAGS) $(LFLAGS) $(LLFLAGSMIN)  
# 
# 
# #-------------------------------------------
# # test programs
# 
# sim_check	: 
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/sim_check.c -o $(OUTDIR)/sim_check $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# 
# cut-sky-test1 : 
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/cut-sky-test1.c -o $(OUTDIR)/cut-sky-test1 $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# 
# test-map : 
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/test-map.c -o $(BIN_DIR)/test-map $(IFLAGS) $(LFLAGS) $(LLFLAGS) 
# 
# test-healpix-th_vals:
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/$@.cpp -o $(BIN_DIR)/$@ $(IFLAGS) $(LFLAGS) $(LLFLAGS) 
# 
# sim17 : 
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/sim17.c -o $(OUTDIR)/sim17 $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# 
# plot_matrix_exp :	#Mscs-common.h Mscs-global-defs.h
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/plot_matrix_exp.c -o $(OUTDIR)/plot_matrix_exp $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# 
# plot_matrix_exp2 :	#Mscs-common.h Mscs-global-defs.h
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/plot_matrix_exp2.c -o $(OUTDIR)/plot_matrix_exp2 $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# 
# check-pixelization :	#Mscs-common.h Mscs-global-defs.h
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/check-pixelization.c -o $(OUTDIR)/check-pixelization $(IFLAGS) $(LLFLAGS) $(LFLAGS) 
# 
# ang2pix-test :
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/ang2pix-test.c -o $(OUTDIR)/ang2pix-test $(IFLAGS) $(LLFLAGS) $(LFLAGS) 	
# 
# almstest:
# 	$(CC) $(OFLAGS)   $(DISTR_DIRcpp)/test-alms.cpp -o $(OUTDIR)/test-alms $(IFLAGS) $(LLFLAGS) $(LFLAGS) 	
# 
# mapplot:
# 	$(CC) $(OFLAGS)  $(DISTR_DIRcpp)/test-plot-map.cpp -o $(OUTDIR)/test-plot-map $(IFLAGS) $(LFLAGS) $(LLFLAGS)
# 
# smoothingTest2d:
# 	$(CC) $(OFLAGS)  $(DISTR_DIRcpp)/smoothingTest2d.cpp -o $(OUTDIR)/smoothingTest2d $(IFLAGS) $(LFLAGS) $(LLFLAGS)

################################################################
################################################################
################################################################
# Define locations of libraries required for linking
#

# provided external libraries
set (VELKB_LIB velKB)
set (SLAREFR_LIB slaRefr)
set (NOVAS_LIB novas)
set (CCSHT3_LIB ccSHT3)

    


#	
# look for OMP library and set compiler flags
#
#find_library (OMP_LIBRARY
#	NAMES Mscsfn
#	PATHS $ENV{MSCS_PROGRAM_DIR}/lib
#	)

################################################################
################################################################
################################################################
# Set extra compiler flags for openmp
#
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

SET( CMAKE_CXX_FLAGS  
	"${CMAKE_CXX_FLAGS} -fopenmp -fPIC"  
	)

set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra") # -fPIC
set(CMAKE_CXX_FLAGS_RELEASE "-O3") # -fPIC

add_definitions(
	-DGO64BIT 
	-D_LARGEFILE_SOURCE 
	-D_FILE_OFFSET_BITS=64 
	-DGENERATE_CVS_VERSION_STRING 
#	-DDEBUG_USING_FFTW2AND3 
	-DHAVE_OPENMP 
#    -DDEBUG_MCMC_PDF
#    -DDEBUG_MCMC_PDF2
#    -Wall
	)
# add -DCMAKE_CXX_FLAGS=-pg when calling cmake to generate gprof-ready binaries
#	-DDEBUG_SUBDOMAIN
#	-DDEBUG_SUBDOMAIN2

################################################################
################################################################
################################################################
# Define list of libraries required for linking
#

set (CPEMS_LIBS
	CPEDS
	Mscsfn
	Mscscore
	Mscsplot
	MscsWMAP
	)

set (CPEMS_HEADERS
	${CPEDSLIB_HEADERS}
	${MSCSCORELIB_HEADERS}
	${MSCSFNLIB_HEADERS}
	${MSCSPLOTLIB_HEADERS}
	${MSCSWMAPLIB_HEADERS}
	)

set (LIBRARIES 
	${CPEMS_LIBS}
	${NOVAS_LIB}
	${CCSHT3_LIB}
	${VELKB_LIB}
	${SLAREFR_LIB}
	${NOVA_LIB}
	${GSL_LIBRARIES}
	${CGAL_LIBRARIES} CGAL
	QtCore
	${FFTW_LIBRARIES}
	${PGPLOT_LIBRARIES}
	${HDF5_LIBRARIES}
	${CFITSIO_LIBRARIES}
	${NOVA_LIBRARIES}
	${PROJ_LIBRARY}
	${ARMADILLO_LIBRARIES}
	)

################################################################
################################################################
################################################################
# Link CPEDS targets
#
target_link_libraries(cosmocalc ${LIBRARIES})
target_link_libraries(integrate ${LIBRARIES})
target_link_libraries(interpolate ${LIBRARIES})
target_link_libraries(compute_CL_ranges ${LIBRARIES})
target_link_libraries(bin_function ${LIBRARIES})
target_link_libraries(coco ${LIBRARIES})


# Link Mscs targets
#
target_link_libraries(select-map-pixels ${LIBRARIES})
target_link_libraries(matrix_oper ${LIBRARIES})
target_link_libraries(calcCMBflux ${LIBRARIES})
target_link_libraries(mkCMBfield ${LIBRARIES})
target_link_libraries(smoothPoints ${LIBRARIES})
target_link_libraries(draw_maps ${LIBRARIES})
target_link_libraries(LMdataFit ${LIBRARIES})
target_link_libraries(function_fft ${LIBRARIES})
target_link_libraries(mkFractal ${LIBRARIES})
target_link_libraries(mk_gauss_map ${LIBRARIES})
target_link_libraries(map_oper ${LIBRARIES} ${LIBRARIES})
target_link_libraries(calculate_Cth ${LIBRARIES} ${LIBRARIES})
target_link_libraries(calculate_minkowski_v ${LIBRARIES} ${LIBRARIES})

# Link Mscs test targets
#
target_link_libraries(test-hdf-openMP2 ${LIBRARIES})
target_link_libraries(test-hdf-openMP ${LIBRARIES})
#target_link_libraries(testfn-openMP ${LIBRARIES})
target_link_libraries(test-mscsfn ${LIBRARIES})
target_link_libraries(test-function-finter ${LIBRARIES})
target_link_libraries(test-mscsfn-hdf5 ${LIBRARIES})
target_link_libraries(test-domain-decomposition ${LIBRARIES})
target_link_libraries(test-function-fft ${LIBRARIES})
target_link_libraries(test-function-derivative ${LIBRARIES})
target_link_libraries(test-fits-read ${LIBRARIES})
target_link_libraries(test-cpedsMCMC-parabola_2param ${LIBRARIES})


################################################################
################################################################
################################################################
# Install targets
#
install (TARGETS 
	bin_function 
	coco 
	compute_CL_ranges 
	cosmocalc 
	integrate 
	interpolate 
#
	calcCMBflux
	draw_maps
	function_fft
	map_oper
	matrix_oper
	mk_gauss_map
	mkCMBfield
	mkFractal
	LMdataFit
	smoothPoints
	calculate_Cth
	calculate_minkowski_v
	DESTINATION bin)

install(TARGETS
    test-hdf-openMP2
    test-hdf-openMP
    test-mscsfn
    test-function-finter
    test-mscsfn-hdf5
    test-domain-decomposition
    test-function-fft
    test-function-derivative
    test-fits-read
    test-cpedsMCMC-parabola_2param
    DESTINATION cpems/tests)
    

install (TARGETS 
	${CPEMS_LIBS}
	DESTINATION lib)


install (FILES
    ${PROJECT_BINARY_DIR}/external_packages/velKB/libvelKB.a
    ${PROJECT_BINARY_DIR}/external_packages/sla-refr/libslaRefr.a
    ${PROJECT_BINARY_DIR}/external_packages/novas/libnovas.a
    ${PROJECT_BINARY_DIR}/external_packages/ccSHT3/libccSHT3.a
    
    DESTINATION lib)


install (FILES 
	${CPEMS_HEADERS}
	DESTINATION include/cpems)

#--prefix=${CMAKE_INSTALL_PREFIX} 
install (CODE "execute_process(COMMAND sudo python setup.py install -f
WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/pyth)")
 


# fetch this into setup.py too
install (FILES 
    ${PROJECT_SOURCE_DIR}/pyth/plot_likelihood.py
    ${PROJECT_SOURCE_DIR}/pyth/plot_function.py
    ${PROJECT_SOURCE_DIR}/pyth/plot_matrix.py
    ${PROJECT_SOURCE_DIR}/pyth/mat_oper.py
    ${PROJECT_SOURCE_DIR}/pyth/mat_oper.py
    
    PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                    GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ
	DESTINATION bin)
#    ${PROJECT_SOURCE_DIR}/pyth/pyCPEDScommonFunctions/cal2jd.py

#install(DIRECTORY DESTINATION $ENV{HOME}/.cpems)
    
#install (TARGETS 
#	${CPEMS_LIBS}
#	DESTINATION lib)

################################################################
################################################################
################################################################
# Make binaries from targets
#
# run "cpack --config CPackConfig.cmake" to make

# build a CPack driven installer package
include (InstallRequiredSystemLibraries)
set (CPACK_RESOURCE_FILE_LICENSE  
     "${CMAKE_CURRENT_SOURCE_DIR}/License.txt")
set (CPACK_PACKAGE_VERSION_MAJOR "1")
set (CPACK_PACKAGE_VERSION_MINOR "0")
include (CPack)
