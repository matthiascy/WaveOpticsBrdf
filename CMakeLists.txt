cmake_minimum_required(VERSION 3.30)
project(WaveOpticsBRDF)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -march=native -fopenmp")

find_package(OpenEXR REQUIRED)

add_executable(gen_brdf
        exrimage.cpp
        exrimage.h
        gaborkernel.cpp
        gaborkernel.h
        genBrdf.cpp
        helpers.h
        spectrum.cpp
        spectrum.h
        waveBrdf.cpp
        waveBrdf.h)

target_include_directories(
        gen_brdf
        PUBLIC
#        /usr/include/OpenEXR
        /usr/include/eigen3
#        /usr/include/Imath
)

target_link_libraries(gen_brdf OpenEXR::OpenEXR)
