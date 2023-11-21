cmake_minimum_required(VERSION 3.26)
project(Projet_Systemes_Distribues)

set(CMAKE_CXX_STANDARD 23)

find_package(MPI REQUIRED)

add_executable(Projet_Systemes_Distribues main.cpp)

target_link_libraries(Projet_Systemes_Distribues MPI::MPI_CXX)