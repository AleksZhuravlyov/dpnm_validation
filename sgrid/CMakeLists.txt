cmake_minimum_required(VERSION 3.17)

set(CMAKE_CXX_STANDARD 17)

project(sgrid)

find_package(pybind11 REQUIRED)

find_package(Eigen3 REQUIRED)

find_package(VTK REQUIRED)

find_package(pugixml)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

pybind11_add_module(${PROJECT_NAME} wrapper.cpp
        Sgrid.cpp SgridVtk.cpp SgridGetSet.cpp)

target_include_directories(${PROJECT_NAME} PUBLIC ${EIGEN3_INCLUDE_DIRS})

target_link_libraries(${PROJECT_NAME} PUBLIC ${VTK_LIBRARIES}
        pugixml)

set_target_properties(${PROJECT_NAME} PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/../..)

