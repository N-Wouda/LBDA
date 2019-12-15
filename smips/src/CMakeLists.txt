cmake_minimum_required(VERSION 3.2)

add_subdirectory(benders)
add_subdirectory(data)
add_subdirectory(decompositions)
add_subdirectory(deqform)
add_subdirectory(master)
add_subdirectory(problem)
add_subdirectory(sub)

add_executable(smips
		main.cpp)

target_link_libraries(smips
		benders
		decompositions
		deqform
		problem
		data
		${GUROBI_C_LIBRARY}
		${GUROBI_CXX_LIBRARY})
