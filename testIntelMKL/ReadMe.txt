testIntelMKL-

The objective is to test Intel's MKL with Cmake and Fortran. You will need Intel's ifort compiler and the MKL library. The Cmake module FindMKL.cmake should assist with library path setup for compiling. On a debian machine, the FindMKL.cmake needs to be placed in the /usr/share/cmake-<version>/Modules/ path.

The blastest.f90 tests the mkl library versus native f77 and f90 style vector and matrix operations for both correctness and speed. On my Core 2 Duo machine, the MKL tradeoff for blas 2 (matrix) operations is about n=200 double precision arrays.

I hope this helps you compute faster and smarter. Thanks to Intel and the CMake guys for a nice scientific computing platform! Please send me an email if you have an improvement or comment.

Charles O'Neill
charles.oneill@gmail.com
7 Sept 2010
