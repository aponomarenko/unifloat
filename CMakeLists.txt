project (libunifloat C)
cmake_minimum_required (VERSION 2.6)

set(INC_DIR "include/unifloat")

set(INCLUDE_BASIC ${INC_DIR}/constants.h ${INC_DIR}/unifloat.h ${INC_DIR}/cstring.h ${INC_DIR}/libunifloat.h ${INC_DIR}/debug.h ${PROJECT_BINARY_DIR}/${INC_DIR}/config.h)
set(INCLUDE_MATH ${INC_DIR}/trig.h ${INC_DIR}/exp.h ${INC_DIR}/bessel.h ${INC_DIR}/gamma.h)
set(INCLUDE_MATH_COMPLEX ${INC_DIR}/unifloat_complex.h ${INC_DIR}/ctrig.h ${INC_DIR}/cexp.h)

set(SRC_BASIC src/unifloat.c src/cstring.c src/debug.c)
set(SRC_MATH src/trig.c src/exp.c src/math.c src/bessel.c src/gamma.c)
set(SRC_MATH_COMPLEX src/unifloat_complex.c src/ctrig.c src/cexp.c)

include_directories ("${PROJECT_SOURCE_DIR}/include")
include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${PROJECT_BINARY_DIR}/include") # to include config.h

set(APPLICATION_NAME ${PROJECT_NAME})

# the version number
set(APPLICATION_VERSION_MAJOR "1")
set(APPLICATION_VERSION_MINOR "0")
set(APPLICATION_VERSION_PATCH "0")
set(APPLICATION_VERSION "${APPLICATION_VERSION_MAJOR}.${APPLICATION_VERSION_MINOR}.${APPLICATION_VERSION_PATCH}")

set(LIBRARY_VERSION "1.0.0")
set(LIBRARY_SOVERSION "1")

# where to look first for cmake modules, before ${CMAKE_ROOT}/Modules/ is checked
set(CMAKE_MODULE_PATH
  ${CMAKE_SOURCE_DIR}/cmake/Modules
)
# add definitions
include(DefineCMakeDefaults)
include(DefineCompilerFlags)

# do not complain about missing prototypes in math.c
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-missing-prototypes")

# disallow in-source build
include(MacroEnsureOutOfSourceBuild)
macro_ensure_out_of_source_build("${PROJECT_NAME} requires an out of source build. Please create a separate build directory and run 'cmake /path/to/${PROJECT_NAME} [options]' there.")

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
OPTION( DEBUG_ON "Enable debuging" ON )  
else()
OPTION( DEBUG_ON "Enable debuging" OFF )  
endif()

configure_file(config.h.cmake ${INC_DIR}/config.h)

SET(BUILD_SHARED_LIBS ON)
add_library (unifloat ${SRC_BASIC} ${SRC_MATH} ${SRC_MATH_COMPLEX})

set_target_properties(
  unifloat
    PROPERTIES
      VERSION
        ${LIBRARY_VERSION}
      SOVERSION
        ${LIBRARY_SOVERSION}
      OUTPUT_NAME
        unifloat
)

# building tests
add_subdirectory (tests)

if (UNIX)
configure_file(libunifloat.pc.cmake ${PROJECT_BINARY_DIR}/libunifloat.pc)
install(
  FILES
    ${PROJECT_BINARY_DIR}/libunifloat.pc
  DESTINATION
    ${CMAKE_INSTALL_PREFIX}/lib/pkgconfig
  COMPONENT
    pkgconfig
)
endif (UNIX)

install (TARGETS unifloat DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")
install (FILES  ${INCLUDE_BASIC} ${INCLUDE_MATH} ${INCLUDE_MATH_COMPLEX} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/unifloat-1.0/unifloat")
