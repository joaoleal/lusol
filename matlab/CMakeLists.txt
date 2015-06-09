##
#        Create matlab files
##
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  # settings for mac os x
  SET(MLARCH maci64)
ELSE()  
  # settins for linux
  SET(MLARCH glnxa64)
ENDIF()

###############################################################################

GET_PROPERTY(CLUSOL_LIB_PATH TARGET clusol PROPERTY LOCATION)
GET_FILENAME_COMPONENT(CLUSOL_LIB ${CLUSOL_LIB_PATH} NAME)

IF("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  SET(INSOURCEBUILD "YES")
ELSE()
  SET(INSOURCEBUILD "NO")
ENDIF()

###############################################################################
#                               MATLAB BUILD
###############################################################################

ADD_CUSTOM_COMMAND(OUTPUT ${CLUSOL_LIB}
                          "clusol.h"
                   COMMAND ${CMAKE_COMMAND} -E copy ${CLUSOL_LIB_PATH} ${CLUSOL_LIB}
                   COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_BINARY_DIR}/../clusol.h" "clusol.h"
                   DEPENDS clusol clusol_header)

IF(NOT ${INSOURCEBUILD})
  SET(LUSOL_BUILD "lusol_build.m")
  ADD_CUSTOM_COMMAND(OUTPUT lusol_build.m
                     COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/lusol_build.m" "lusol_build.m")
ENDIF()
                   
ADD_CUSTOM_COMMAND(OUTPUT "libclusol_thunk_${MLARCH}${CMAKE_SHARED_LIBRARY_SUFFIX}"
                          "libclusol_proto_${MLARCH}.m"
                   COMMAND ${Matlab_MAIN_PROGRAM} ARGS -nojvm -nodisplay -nosplash -nodesktop -r \"lusol_build\; exit\;\"
                   DEPENDS ${CLUSOL_LIB} "clusol.h" ${LUSOL_BUILD}
                   COMMENT "Creating matlab files")

ADD_CUSTOM_TARGET(matlab ALL
                  DEPENDS "libclusol_thunk_${MLARCH}${CMAKE_SHARED_LIBRARY_SUFFIX}")

###############################################################################
#                                TESTING
###############################################################################

IF(NOT ${INSOURCEBUILD})
  SET(LUSOL_BUILD "lusol_build.m")
  ADD_CUSTOM_COMMAND(OUTPUT lusol_build.m
                     COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/lusol_build.m" "lusol_build.m")

  SET(LUSOL_TEST "lusol_test.m")
  ADD_CUSTOM_COMMAND(OUTPUT lusol_test.m
                     COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/lusol_test.m" "lusol_test.m")

  SET(LUSOL_OBJ "lusol_obj.m")
  ADD_CUSTOM_COMMAND(OUTPUT lusol_obj.m
                     COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/lusol_obj.m" "lusol_obj.m")

  FILE(COPY "${CMAKE_CURRENT_SOURCE_DIR}/test_data" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
ENDIF()

ADD_CUSTOM_TARGET(test_lusol ALL
                  COMMAND ${Matlab_MAIN_PROGRAM} ARGS -nojvm -nodisplay -nosplash -nodesktop -r \"lusol_test\; exit\;\"
                  DEPENDS ${CLUSOL_LIB} ${LUSOL_BUILD} ${LUSOL_TEST} ${LUSOL_OBJ} ${MATRIX_FILES}
                  COMMENT "Testing LUSOL")
