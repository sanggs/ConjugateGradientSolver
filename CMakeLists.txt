cmake_minimum_required (VERSION 2.8.11)
set (CMAKE_CXX_STANDARD 11)

project (ConjugateGradient)

set (CMAKE_CXX_FLAGS "-Wno-deprecated ${CMAKE_CXX_FLAGS}")

add_definitions(-DBUILD_COMPONENT_SRC_PREFIX="" -DBUILD_OPTLEVEL_DEV)

set(USD_ROOT_DIR /opt/local/USD)
set(USD_LIB_DIR ${USD_ROOT_DIR}/lib/)
set(USD_INC_DIR ${USD_ROOT_DIR}/include/)

add_executable (ConjugateGradient
        src/main.cpp
        src/AnimatedMesh.h
        src/LatticeMesh.h
        src/CGVectorWrapper.h
        src/Laplacian.h
        src/ConjugateGradientSolver.h
        )

find_package(PythonLibs)

set(EIGEN3_INC_DIR /usr/local/include/eigen3)

target_include_directories(ConjugateGradient
        PUBLIC
            ${USD_INC_DIR}
            ${EIGEN3_INC_DIR}
            ${PYTHON_INCLUDE_PATH}
        )

target_link_libraries(
        ConjugateGradient
                ${USD_LIB_DIR}/libgf.dylib
                ${USD_LIB_DIR}/libsdf.dylib
                ${USD_LIB_DIR}/libtf.dylib
                ${USD_LIB_DIR}/libusd.dylib
                ${USD_LIB_DIR}/libusdGeom.dylib
                ${USD_LIB_DIR}/libvt.dylib
                ${USD_LIB_DIR}/libboost_python.dylib # todo find library
                ${PYTHON_LIBRARY}
        )
