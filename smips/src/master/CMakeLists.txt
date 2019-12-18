cmake_minimum_required(VERSION 3.2)

set(SOURCES
		addcut.cpp
		copy_ctor.cpp
		destructor.cpp
		master.cpp
		solve.cpp)

add_library(master ${SOURCES})
