cmake_minimum_required (VERSION 2.6)
project (tests)

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-builtin -lm")
add_executable(tests EXCLUDE_FROM_ALL main.c)

target_link_libraries (tests unifloat) 
