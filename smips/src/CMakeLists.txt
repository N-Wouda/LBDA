cmake_minimum_required(VERSION 3.2)

find_package(Gurobi REQUIRED)
include_directories(${GUROBI_INCLUDE_DIRS})

add_subdirectory(cuts)
add_subdirectory(data)
add_subdirectory(deterministicequivalent)
add_subdirectory(masterproblem)
add_subdirectory(problem)
add_subdirectory(solve)
add_subdirectory(subproblem)

add_executable(smips main.cpp)
target_link_libraries(smips
		cuts
		deterministicequivalent
		problem
		solve
		data
		${GUROBI_LIBRARIES})
