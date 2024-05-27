cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)
project(hpc LANGUAGES C CXX)

set(SOURCES src/hpc.c)

add_library(${PROJECT_NAME})
target_sources(${PROJECT_NAME} PRIVATE ${SOURCES})
target_include_directories(${PROJECT_NAME} PUBLIC include)

if(WITH_HPC_SINGLE_KX)
    target_compile_definitions(${PROJECT_NAME} PRIVATE HPC_USE_SINGLE_KX=${WITH_HPC_SINGLE_KX})
endif()
