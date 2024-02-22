cmake_minimum_required(VERSION 3.20)
project(IntFFT)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

add_executable(IntFFT src/main.cpp)
