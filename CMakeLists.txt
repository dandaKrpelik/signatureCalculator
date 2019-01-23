project(pMH4ISL)

cmake_minimum_required(VERSION 2.8)

find_package(OpenMP REQUIRED)

OPTION( WITH_PYPMH "Build pyIPSIG ?" ON)



set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -O2 -DCMAKE_C_FLAGS=-fopenmp -DCMAKE_EXE_LINKER_FLAGS=-fopenmp")
add_definitions("-std=c++11 -fopenmp")
include_directories("${CMAKE_SOURCE_DIR}/include")

if (WITH_PYPMH)

    add_subdirectory(pybind11)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/pylibs)
    pybind11_add_module(pyIPSIG ${CMAKE_SOURCE_DIR}/pyIPSIG.cpp)

    add_custom_command(OUTPUT "${PROJECT_BINARY_DIR}/pylibs/test.py"
                COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/pylibs/test.py ${PROJECT_BINARY_DIR}/pylibs/test.py
                DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/pylibs/test.py"
                COMMENT "copying test.py"
    )

    set(dep ${PROJECT_BINARY_DIR}/pylibs/test.py)

    add_custom_target(pylib_copy_module ALL DEPENDS ${dep})
 
endif()

#add_definitions("-std=c++11 -fopenmp")
add_executable(main main.cpp)
