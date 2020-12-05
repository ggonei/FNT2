# CMakeList.txt : CMake project for FNT2, include source and define project specific logic here.
cmake_minimum_required (VERSION 3.8)

file(GLOB SOURCES
	include/*.h
	src/*.cpp
)
add_executable (FNT2 ${SOURCES})