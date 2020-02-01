cmake_minimum_required(VERSION 3.2)

set(SOURCES
		copyconstructor.cpp
		destructor.cpp
		solve.cpp
		gominfo.cpp
		multipliers.cpp
		subproblem.cpp
		update.cpp)

add_library(sub ${SOURCES})
