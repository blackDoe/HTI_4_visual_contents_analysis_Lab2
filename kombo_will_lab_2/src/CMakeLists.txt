cmake_minimum_required(VERSION 3.7)
project(HTI4_visual_contents_analysis_Lab2)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        globals.h
        pdefilter.c
        pdefilter.h
        permal.c
        utilities.c
        utilities.h)

add_executable(HTI4_visual_contents_analysis_Lab2 ${SOURCE_FILES})