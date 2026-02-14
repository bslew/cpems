cmake_minimum_required (VERSION 2.6)
set (CMAKE_CXX_STANDARD 17)
project (cpems)

# add_executable(func_test_example "func_test_example.cpp")
add_executable(cpems_arithmetic "cpems_arithmetic.cpp")
add_executable(test-cpeds-directions "test-cpeds-directions.cpp")

SET (LIB_UNIT_TEST
	gtest
	gtest_main
)

SET (LIB_UNIT_TEST_MIN
	gtest
)

target_link_libraries(cpems_arithmetic 
	${LIBRARIES} 
	${LIB_UNIT_TEST}
)

target_link_libraries(test-cpeds-directions 	${LIBRARIES} 	${LIB_UNIT_TEST})
