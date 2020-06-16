cmake_minimum_required (VERSION 2.6)
set (CMAKE_CXX_STANDARD 17)
project (cpems)

#add_executable(cpeds_get_healpix_nested_pixels "cpeds_get_healpix_nested_pixels.cpp")

SET (LIB_BOOST
	boost_filesystem
	boost_system
)

#target_link_libraries(cpeds_get_healpix_nested_pixels ${LIBRARIES} )

add_subdirectory(unit-tests)
