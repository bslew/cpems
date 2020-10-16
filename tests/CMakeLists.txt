
add_executable(testDirections "testDirections.cpp")

SET (LIB_BOOST
	boost_filesystem
	boost_system
)

target_link_libraries(testDirections ${LIBRARIES} )

add_subdirectory(unit-tests)
