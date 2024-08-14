cmake_minimum_required(VERSION 3.1...3.15)

project(ALC_TRAJECTORY VERSION 1.4 
        LANGUAGES Fortran)
set(AUTHOR "Ivan Scivetti")
set(AUTHOR_DETAILS "ivan.scivetti@stfc.ac.uk")
set(DESCRIPTION "A code for analysis of trajectories")
set(target_name alc_trajectory.x)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")

set(BUILDER "" CACHE STRING "Individual who has built the code")
cmake_host_system_information(RESULT AH QUERY FQDN)
set(HOST "${AH}" CACHE STRING "Name of the host machine used to build the code")

# wildcard for source
file(GLOB SOURCES source/*.F90)  # assgin all the *.F90 to the variable SOURCES, this is global variable.

# path to the folder with the executable
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
# path to the folder with Fortran modules
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# Name executable
add_executable(alc_trajectory.x ${SOURCES})

#Directives for buling the code automatically
if (NOT FLAGS_SET)
  if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")

    if (${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER "6.5" )
      set(CMAKE_Fortran_FLAGS_DEBUG "-g -Wextra -Wuse-without-only -frecursive -frecord-gcc-switches -O0 -std=f2008 -pedantic -fbacktrace -fcheck=all -finit-integer=2147483647 -finit-real=snan -finit-logical=true -finit-character=42 -finit-derived -ffpe-trap=invalid,zero,overflow -fdump-core -fstack-protector-all -Wall -pipe" CACHE STRING "Flags used by the GNU-Fortran compiler for DEBUG option." FORCE)
    else()
      set(CMAKE_Fortran_FLAGS_DEBUG "-g -Wextra -Wuse-without-only -frecord-gcc-switches -O0 -std=f2008 -pedantic -fbacktrace -fcheck=all -finit-integer=2147483647 -finit-real=snan -finit-logical=true -finit-character=42  -ffpe-trap=invalid,zero,overflow -fdump-core -fstack-protector-all -Wall -pipe" CACHE STRING "Flags used by the GNU-Fortran compiler for the DEBUG option." FORCE)
    endif()
    set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast -ftree-vectorize -funroll-loops -ffast-math" CACHE STRING "Flags used by the GNU-Fortran compiler for the RELEASE option." FORCE)

    if (${CMAKE_Fortran_COMPILER_VERSION} VERSION_LESS "5.4.0")
      message(FATAL_ERROR "***ERROR: Available GNU compiler version was not tested against Fortran2008 standards. Recommended minimum version: 5.4.0. The user can try building the code manually... at its own risk!")
    endif()

  elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")

   if (${CMAKE_Fortran_COMPILER_VERSION} VERSION_LESS "16.0.1")
     message(FATAL_ERROR "***ERROR: Available Intel compiler version was not tested against Fortran2008 standards. Recommended minimum version: 16.0.1. The user could try building the code manually... at its own risk!")
   endif()

   set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -stand f08 -traceback -C -check all,nouninit -ftrapuv -init=snan -init=arrays" CACHE STRING "Flags used by the Intel-Fortran compiler for the DEBUG option." FORCE)
   set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -C -check nouninit" CACHE STRING "Flags used by the Intel-Fortran compiler for the RELEASE option." FORCE)

  endif()
  set(FLAGS_SET 1 CACHE INTERNAL "Flags have been defined")
endif()


if(WITH_TESTING)
  include(CTest)
  message(STATUS "Building with implemented tests...")

  # Define a macro
  macro (do_test arg)
    # Look for the .tar file to perform the testing
    if (EXISTS ${CMAKE_SOURCE_DIR}/CI-tests/${arg}.tar)
      message(STATUS "Test file ${arg}.tar unpacked")
    else()
      message(FATAL_ERROR "Missing ${arg}.tar test")
    endif()

    # Some filing operartion 
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/Testing/new/${arg})
    set(TEST_NEW_PATH "${CMAKE_BINARY_DIR}/Testing/new/${arg}")
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/Testing/reference/${arg})
    set(TEST_REF_PATH "${CMAKE_BINARY_DIR}/Testing/reference/${arg}")
    file(COPY ${CMAKE_SOURCE_DIR}/CI-tests/${arg}.tar  DESTINATION ${TEST_NEW_PATH})
    file(COPY ${CMAKE_SOURCE_DIR}/CI-tests/${arg}.tar  DESTINATION ${TEST_REF_PATH})
    execute_process(COMMAND ${CMAKE_COMMAND} -E tar xvf "${arg}.tar"  WORKING_DIRECTORY ${TEST_NEW_PATH} OUTPUT_QUIET)
    execute_process(COMMAND ${CMAKE_COMMAND} -E tar xvf "${arg}.tar"  WORKING_DIRECTORY ${TEST_REF_PATH} OUTPUT_QUIET)
    execute_process(COMMAND rm "${arg}.tar"  WORKING_DIRECTORY ${TEST_NEW_PATH} OUTPUT_QUIET)
    execute_process(COMMAND sh ${CMAKE_SOURCE_DIR}/scripts/clean_data.sh  WORKING_DIRECTORY ${TEST_NEW_PATH} OUTPUT_QUIET)
    execute_process(COMMAND rm "${arg}.tar"  WORKING_DIRECTORY ${TEST_REF_PATH} OUTPUT_QUIET)
    execute_process(COMMAND bash ${CMAKE_SOURCE_DIR}/scripts/clean_data.sh  WORKING_DIRECTORY ${TEST_NEW_PATH} OUTPUT_QUIET)
    # Define the test
    add_test(NAME ${arg} COMMAND bash ${CMAKE_SOURCE_DIR}/tools/launch-test.sh "${alc_trajectory}" ${TEST_REF_PATH} WORKING_DIRECTORY ${TEST_NEW_PATH})
  endmacro (do_test)

  # Name executable + whole path for execution
  set(alc_trajectory "${CMAKE_BINARY_DIR}/bin/${target_name}")

  foreach(test RANGE 1 11 1)
    do_test("test${test}")
  endforeach()

else()
  message(STATUS "Build without testing")
endif()
