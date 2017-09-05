################################################################
################################################################
################################################################
# Set project name
#
project (cpems CXX)


################################################################
################################################################
################################################################
# Define source files
#
#set (VELKBLIB_SOURCES
#	velKB velKB.f
#	)

################################################################
################################################################
################################################################
# Add library targets
#
#add_library(velKB STATIC ${VELKBLIB_SOURCES} )



#ADD_CUSTOM_COMMAND(
#        WORKING_DIRECTORY external_packages/velKB
#        OUTPUT external_packages/velKB
#        COMMAND env CXXFLAGS=${ROOTS_COMPILER_FLAGS} ./configure ${TARGET_FLAG} ${ROOTS_HOST_OPTION}
#        COMMAND make
#        COMMENT "Building velKB"
#        VERBATIM
#    )
