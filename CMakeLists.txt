cmake_minimum_required(VERSION 3.26)

#partie calcul
project(Projet_Systemes_Distribues_calcul)
set(CMAKE_CXX_STANDARD 23)
find_package(MPI REQUIRED)
add_executable(Projet_Systemes_Distribues_calcul calcul/calcul.cpp)
target_link_libraries(Projet_Systemes_Distribues_calcul MPI::MPI_CXX)

#partie visualisation
project(Projet_Systemes_Distribues_visualisation)
set(CMAKE_CXX_STANDARD 23)
find_package(GLUT REQUIRED)
find_package(OpenGL REQUIRED)
add_executable(Projet_Systemes_Distribues_visualisation visualisation/visualisation.cpp)
target_link_libraries(Projet_Systemes_Distribues_visualisation GLUT::GLUT OpenGL::GLU OpenGL::GL)