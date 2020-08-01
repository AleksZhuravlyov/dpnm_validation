cmake_minimum_required(VERSION 3.17)

set(CMAKE_CXX_STANDARD 17)

project(diffusion)

find_package(pybind11 REQUIRED)

find_package(Eigen3 REQUIRED)

pybind11_add_module(${PROJECT_NAME} wrapper.cpp
        math/Props.cpp
        math/Local.cpp
        math/Convective.cpp
        Equation.cpp)

include_directories(.)

target_link_libraries(${PROJECT_NAME} ${EIGEN3_LIBRARIES})

target_include_directories(${PROJECT_NAME} PUBLIC ${EIGEN3_INCLUDE_DIRS})

set_target_properties(${PROJECT_NAME} PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/../..)