Name: ${APPLICATION_NAME}
Description: A C library for fast high-precision floating-point computations
Version: ${APPLICATION_VERSION}
Libs: -L${CMAKE_INSTALL_PREFIX}/lib -lunifloat
Cflags: -fno-builtin -I${CMAKE_INSTALL_PREFIX}/include/unifloat-1.0
