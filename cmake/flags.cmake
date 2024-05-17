# Directive for compilation
if (NOT FLAGS_SET)
  if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")

    if (${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER "6.5" )
      set(CMAKE_Fortran_FLAGS_DEBUG "-g -Wextra -Wuse-without-only -frecord-gcc-switches -O0 -std=f2008 -pedantic -fbacktrace -fcheck=all -finit-integer=2147483647 -finit-real=snan -finit-logical=true -finit-character=42 -finit-derived -ffpe-trap=invalid,zero,overflow -fdump-core -fstack-protector-all -Wall -pipe" CACHE STRING "Flags used by the GNU-Fortran compiler for DEBUG option." FORCE)
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


    set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -stand f08 -traceback -C -fp-stack-check -ftrapuv -init=snan -init=arrays" CACHE STRING "Flags used by the Intel-Fortran compiler for the DEBUG option." FORCE)
    set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast" CACHE STRING "Flags used by the Intel-Fortran compiler for the RELEASE option." FORCE)

  endif()
  set(FLAGS_SET 1 CACHE INTERNAL "Flags have been defined")
endif()
