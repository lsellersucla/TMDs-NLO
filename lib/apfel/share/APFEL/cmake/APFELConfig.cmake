
####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was APFELConfig.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

SET(APFEL_VERSION )
SET(APFEL_VERSION_MAJOR  3)
SET(APFEL_VERSION_MINOR  1)
SET(APFEL_VERSION_PATCH  1)

set_and_check(APFEL_INCLUDE_DIR ${PACKAGE_PREFIX_DIR}/include)
find_library(APFEL_LIB NAMES APFEL HINTS ${PACKAGE_PREFIX_DIR}/lib)
find_library(APFEL_EVOL_LIB NAMES APFELevol HINTS ${PACKAGE_PREFIX_DIR}/lib)

set(APFEL_LIBRARIES ${APFEL_LIB} ${APFEL_EVOL_LIB})

include(${CMAKE_CURRENT_LIST_DIR}/APFELTargets.cmake)

