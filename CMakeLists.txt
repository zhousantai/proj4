cmake_minimum_required(VERSION 3.23)
add_definitions(-DWITH_NEON)
set(CMAKE_C_COMPILER /opt/homebrew/bin/gcc-12)
add_compile_options(-O3)
#add_definitions(-DWITH_NEON)
project(proj4)
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -L/opt/homebrew/opt/libomp/lib" )
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I/opt/homebrew/opt/libomp/include" )
set(CMAKE_C_STANDARD 99)

set(LDFLAGS "-L/opt/homebrew/opt/openblas/lib")
set(USE_BLAS OpenBLAS)
set(BLAS_LIBRARY /opt/homebrew/opt/openblas/lib/libblas.dylib)




#OPTION (USE_OpenMP"Use OpenMP" ON)
#IF(USE_OpenMP)
#    FIND_PACKAGE(OpenMP)
#    IF(OPENMP_FOUND)
#        SET(CMAKE_C_FLAGS"${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    ENDIF()
#ENDIF()
#
##target_link_libraries(proj4 OpenMP::OpenMP_C)
#target_link_libraries(proj4 OpenMP::OpenMP_C)
#target_link_libraries(proj4 ${BLAS_LIBRARY})

set(OpenBlas_DIR /opt/homebrew/opt/openblas/lib)
set(OpenBlas_INC /opt/homebrew/opt/openblas/include)

# 寻找OpenCV库
find_package(OpenCV REQUIRED)
# 添加头文件
include_directories( ${OpenCV_INCLUDE_DIRS} )
include_directories(.)
include_directories(${OpenBlas_INC})
link_directories(${OpenBlas_DIR})
add_executable(proj4 main.c )
target_link_libraries(proj4 libopenblas.a)

set(OpenCV_LIBS opencv_core opencv_imgproc opencv_highgui opencv_imgcodecs)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message("OpenMP found.")
    target_link_libraries(proj4 OpenMP::OpenMP_C)
endif()
    target_link_libraries(proj4 ${OpenCV_LIBS})