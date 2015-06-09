CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

###########################################################################
#         Project information
###########################################################################
PROJECT(LUSOL Fortran C)

# global needed variables
SET(APPLICATION_NAME ${PROJECT_NAME})

SET(APPLICATION_VERSION_MAJOR "1")
SET(APPLICATION_VERSION_MINOR "0")
SET(APPLICATION_VERSION_PATCH "0")

SET(APPLICATION_VERSION "${APPLICATION_VERSION_MAJOR}.${APPLICATION_VERSION_MINOR}.${APPLICATION_VERSION_PATCH}")
SET(PACKAGE_VERSION ${APPLICATION_VERSION})
SET(PACKAGE_URL "https://github.com/nwh/lusol")


###########################################################################
#         Compilation flags
###########################################################################
IF("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  SET(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -ggdb")
  SET(CMAKE_Fortran_FLAGS_RELEASE "-O3 -g -ggdb")

  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-integer-8") # is this really needed?
ENDIF()


###########################################################################
#         General installation information
###########################################################################

SET(CMAKE_INSTALL_PREFIX "/usr" CACHE PATH "Default install path prefix")
SET(INSTALL_LIBRARY_LOCATION "${CMAKE_INSTALL_PREFIX}/lib")
SET(INSTALL_RUNTIME_LOCATION "${CMAKE_INSTALL_PREFIX}/bin")
SET(INSTALL_EXPORT_LOCATION "${CMAKE_INSTALL_PREFIX}/share/lusol/cmake")
SET(INSTALL_INCLUDE_LOCATION "${CMAKE_INSTALL_PREFIX}/include/lusol")
#SET(INSTALL_LIBRARY_PKG_LOCATION "${INSTALL_LIBRARY_LOCATION}/pkgconfig")
SET(INSTALL_SOURCE_LOCATION "${CMAKE_INSTALL_PREFIX}/share/lusol/src")


###########################################################################
#         Dependencies
###########################################################################

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

#
#             BLAS
#
SET(LUSOL_BLAS_VENDOR "Generic" CACHE STRING "BLAS library")
SET_PROPERTY(CACHE LUSOL_BLAS_VENDOR PROPERTY STRINGS 
             Generic OpenBLAS Goto ATLAS PhiPACK CXML DXML SunPerf SCSL SGIMATH IBMESSL ACML ACML_MP ACML_GPU Intel)
IF(LUSOL_BLAS_VENDOR STREQUAL "")
  SET(BLA_VENDOR CACHE "Generic")
ELSE()
  SET(BLA_VENDOR ${LUSOL_BLAS_VENDOR})
ENDIF()

FIND_PACKAGE(BLAS REQUIRED)
IF(${BLAS_FOUND})
  MESSAGE(STATUS "BLAS: ${BLAS_LIBRARIES}")
ENDIF()

#
#             Matlab
#

SET(MATLAB_FIND_DEBUG "YES")

FIND_PACKAGE(MATLAB COMPONENTS MAIN_PROGRAM)

IF(NOT Matlab_FOUND)
  MESSAGE(WARNING "Matlab not found!\n Please define variable Matlab_ROOT_DIR if you would like to build LUSOL for Matlab.")
  #SET(Matlab_ROOT_DIR "" CACHE PATH "Path to the root instalation directory of Matlab")
ELSEIF(NOT Matlab_MAIN_PROGRAM_FOUND)
  MESSAGE(WARNING "Matlab executable not found!")
ENDIF()

###########################################################################
#         Sources
###########################################################################

# list of Fortran code files
SET(F90_FILES "src/lusol.f90"
              "src/lusol8b.f90"
              "src/lusol7b.f90"
              "src/lusol6b.f90"
              "src/lusol_precision.f90")

# list of interface specification files
SET(INTERFACE_FILES "gen/interface.py"
                    "gen/interface_files.org"
                    "gen/lu1fac.org"
                    "gen/lu6mul.org"
                    "gen/lu6sol.org"
                    "gen/lu8adc.org"
                    "gen/lu8adr.org"
                    "gen/lu8dlc.org"
                    "gen/lu8dlr.org"
                    "gen/lu8mod.org"
                    "gen/lu8rpc.org"
                    "gen/lu8rpr.org")

##
#        Create C source and header
##
ADD_CUSTOM_COMMAND(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/clusol.c"
                   COMMAND ./gen/interface.py ARGS -i gen/interface_files.org -o "${CMAKE_CURRENT_BINARY_DIR}/clusol.c" -t source
                   WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
                   COMMENT "Creating C source")
                   
ADD_CUSTOM_COMMAND(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/clusol.h"
                   COMMAND ./gen/interface.py ARGS -i gen/interface_files.org -o "${CMAKE_CURRENT_BINARY_DIR}/clusol.h" -t header
                   WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
                   COMMENT "Creating C header")

ADD_CUSTOM_TARGET(clusol_header
                  DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/clusol.h")

##
#        Create lusol shared library
##
ADD_LIBRARY(clusol SHARED
            "${CMAKE_CURRENT_BINARY_DIR}/clusol.c"
            ${F90_FILES})

TARGET_LINK_LIBRARIES(clusol ${BLAS_LIBRARIES})
ADD_DEPENDENCIES(clusol clusol_header)

# dynamic library version
SET_TARGET_PROPERTIES(clusol
                      PROPERTIES SOVERSION ${APPLICATION_VERSION})

# includes for source compilation with generated C sources
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

##
#        Matlab
##

#SET(Matlab_MAIN_PROGRAM "~/MATLAB/R2015a/bin/matlab" CACHE PATH "Path to the root instalation directory of Matlab")
MESSAGE(WARNING "Matlab_MAIN_PROGRAM=${Matlab_MAIN_PROGRAM}")
IF(Matlab_MAIN_PROGRAM_FOUND)
  ADD_SUBDIRECTORY(matlab)
ENDIF()
                    
###########################################################################
#         export some targets for in tree building
###########################################################################

EXPORT(TARGETS clusol
       FILE "${CMAKE_BINARY_DIR}/lusol.cmake")

###########################################################################
#         Installation
###########################################################################
INSTALL(TARGETS clusol EXPORT clusol
        LIBRARY DESTINATION ${INSTALL_LIBRARY_LOCATION})
        
INSTALL(EXPORT clusol
        DESTINATION "${INSTALL_EXPORT_LOCATION}")

INSTALL(FILES "${CMAKE_CURRENT_BINARY_DIR}/clusol.h"
        DESTINATION "${INSTALL_INCLUDE_LOCATION}")

INSTALL(DIRECTORY "${CMAKE_SOURCE_DIR}/src/"
        DESTINATION "${INSTALL_INCLUDE_LOCATION}"
        FILES_MATCHING PATTERN "*.f")

# also add the fortran files since they have good documentation
INSTALL(DIRECTORY "${CMAKE_SOURCE_DIR}/src/"
        DESTINATION "${INSTALL_INCLUDE_LOCATION}"
        FILES_MATCHING PATTERN "*.f90")