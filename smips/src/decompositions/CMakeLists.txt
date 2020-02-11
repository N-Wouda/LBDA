cmake_minimum_required(VERSION 3.2)

add_subdirectory(gomory)
add_subdirectory(lagrangian)
add_subdirectory(relaxation)

add_library(decompositions empty.cpp)
target_link_libraries(decompositions
		gomory
		lagrangian
		relaxation)
