cmake_minimum_required(VERSION 3.5)
project(ionocell)

set(Y_SDK "$ENV{HOME}/program/yttrium/sdk")
message( STATUS "Y_SDK=${Y_SDK}" )

include(${Y_SDK}/etc/y.cmake)

include_directories(${Y_SDK}/include)
link_directories(${Y_SDK}/lib)

add_executable(ph-vs-time ph-vs-time/main.cpp)
target_link_libraries(ph-vs-time y-mkl y)

