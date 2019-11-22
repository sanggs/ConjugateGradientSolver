cmake_minimum_required (VERSION 2.8.11)
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_PREFIX_PATH "/usr/local/include/libtorch/share/cmake/Torch")

set (CMAKE_CXX_FLAGS "-Wno-deprecated ${CMAKE_CXX_FLAGS}")
set (LDFLAGS "-L/usr/local/opt/llvm/lib -Wl,-rpath,/usr/local/opt/llvm/lib")
set (GLIBCXX_USE_CXX11_ABI 0)
set (CMAKE_EXE_LINKER_FLAGS ${LDFLAGS} ${GLIBCXX_USE_CXX11_ABI})

add_definitions(-DBUILD_COMPONENT_SRC_PREFIX="" -DBUILD_OPTLEVEL_DEV)

set(USD_ROOT_DIR /opt/local/USD)
set(USD_LIB_DIR ${USD_ROOT_DIR}/lib/)
set(USD_INC_DIR ${USD_ROOT_DIR}/include/)

project (ConjugateGradient)

add_executable (ConjugateGradient
        src/main.cpp
        src/AnimatedMesh.h
        src/LatticeMesh.h
        src/CGVectorWrapper.h
        src/Laplacian.h
        src/ConjugateGradientSolver.h
        src/JacobiSolver.h
        )

find_package(PythonLibs)
find_package(Torch REQUIRED)
set(Torch_DIR /usr/local/include/libtorch)
set(EIGEN3_INC_DIR /usr/local/include/eigen3)
#set(TORCH_INCLUDE_PATH /usr/local/include/libtorch/include)
#set(TORCH_LIBRARIES /usr/local/include/libtorch/lib)
set(PYTHON_INCLUDE_PATH /Library/Frameworks/Python.framework/Versions/3.7/include/python3.7m)

target_include_directories(ConjugateGradient
        PUBLIC
            ${USD_INC_DIR}
            ${EIGEN3_INC_DIR}
            ${PYTHON_INCLUDE_PATH}
            ${TORCH_INCLUDE_PATH}
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
                ${TORCH_LIBRARIES}
        )
