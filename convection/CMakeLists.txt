cmake_minimum_required(VERSION 3.17)

set(CMAKE_CXX_STANDARD 17)

project(convection)

find_package(pybind11 REQUIRED)

set(CMAKE_EXE_LINKER_FLAGS -Wl,-execute,-undefined,dynamic_lookup,-w)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(lnInclude)
include_directories($ENV{FOAM_SRC}/finiteVolume/lnInclude)
include_directories($ENV{FOAM_SRC}/OpenFOAM/lnInclude)
include_directories($ENV{FOAM_SRC}/OSspecific/POSIX/lnInclude)
include_directories($ENV{FOAM_SRC}/meshTools/lnInclude/)

add_definitions(
        -D$ENV{WM_ARCH}
        -DWM_ARCH_OPTION=$ENV{WM_ARCH_OPTION}
        -DWM_$ENV{WM_PRECISION_OPTION}
        -DWM_LABEL_SIZE=$ENV{WM_LABEL_SIZE}
        -DNoRepository -m64 -fPIC
)

pybind11_add_module(${PROJECT_NAME} wrapper.cpp convection.cpp)

target_link_libraries(${PROJECT_NAME} PRIVATE
        $ENV{FOAM_LIBBIN}/libfiniteVolume.dylib
        $ENV{FOAM_LIBBIN}/libOpenFOAM.dylib
        $ENV{FOAM_LIBBIN}/libmeshTools.dylib dl m)

set_target_properties(${PROJECT_NAME} PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/../..)
