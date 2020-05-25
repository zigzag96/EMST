# Set the minimum version of CMake that can be used
# To find the cmake version run
# $ cmake --version
cmake_minimum_required(VERSION 3.5)
set(CMAKE_CXX_STANDARD 11)

# Set the project name
project (parallel_EMST)

FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
message("OPENMP FOUND")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Create a sources variable with a link to all cpp files to compile
set(SOURCES
    src/main.cpp
    src/iofile.cpp
    src/Morton.cpp
    src/wspd.cpp
)

add_executable(parallel_EMST ${SOURCES})
target_include_directories(parallel_EMST
    PRIVATE 
        ${PROJECT_SOURCE_DIR}/include
)