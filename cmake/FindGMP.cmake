# Try to find the GNU Multiple Precision Arithmetic Library (GMP)
# See http://gmplib.org/

if (GMP_INCLUDES AND GMP_LIBRARIES)
  set(GMP_FIND_QUIETLY TRUE)
endif (GMP_INCLUDES AND GMP_LIBRARIES)

find_path(GMP_INCLUDES
	NAMES
	gmp.h
	PATHS
	$ENV{gmp_root}
	$ENV{GMPDIR}
	${INCLUDE_INSTALL_DIR}
	PATH_SUFFIXES "include"
)

find_library(GMP_LIBRARIES gmp 
	PATHS 
	$ENV{GMPDIR} 
	$ENV{gmp_root}
	${LIB_INSTALL_DIR}
	PATH_SUFFIXES lib lib64
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMP DEFAULT_MSG
                                  GMP_INCLUDES GMP_LIBRARIES)
mark_as_advanced(GMP_INCLUDES GMP_LIBRARIES)