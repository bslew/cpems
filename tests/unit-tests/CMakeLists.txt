cmake_minimum_required (VERSION 2.6)
set (CMAKE_CXX_STANDARD 17)
project (cpems)

add_executable(cpeds_get_healpix_nested_pixels "cpeds_get_healpix_nested_pixels.cpp")
add_executable(test-smoothing "test-smoothing.cpp")

SET (LIB_BOOST
	boost_program_options
	boost_filesystem
	boost_system
)

target_link_libraries(cpeds_get_healpix_nested_pixels ${LIBRARIES} )
target_link_libraries(test-smoothing 
	${LIBRARIES} 
	${LIB_BOOST}
	)