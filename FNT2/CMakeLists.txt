# CMakeList.txt : CMake project for FNT2, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

file(GLOB SOURCES
	include/*.h
	src/*.cpp
)
# Add source to this project's executable.
add_executable (FNT2 ${SOURCES})

# TODO: Add tests and install targets if needed.
