cmake_minimum_required(VERSION 3.2)

set(SOURCES
		solve.cpp)

add_library(solve ${HEADERS} ${SOURCES})
target_link_libraries(solve
		master
		sub)
