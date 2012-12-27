
SET(NAME IntelMKL)
PROJECT (${NAME} Fortran)

cmake_minimum_required(VERSION 2.6)

find_package(MKL)

SET(SRCS blastest.f90 Timing.f90 ${MKL_INCLUDE_DIR}/blas.f90
 )


ADD_EXECUTABLE(${NAME}  ${SRCS})

IF(UNIX)
  SET(CMAKE_Fortran_FLAGS "-u")
  SET(CMAKE_Fortran_FLAGS_RELEASE " -O3 -xHOST -ipo -no-prec-div -static -align all -fno-alias")
  include_directories(${MKL_INCLUDE_DIR})
  TARGET_LINK_LIBRARIES(${NAME} mkl_blas95 mkl_intel mkl_intel_thread mkl_core iomp5 pthread )
ENDIF(UNIX)

