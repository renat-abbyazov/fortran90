fortran90
=========

http://fortran90.org source code examples

fcython_mesh/ is for the illustration purpose of http://fortran90.org/src/best-practices.html#interfacing-with-python

fcython_qpck_1/ is for the same goal, but one need here to wrap a function written on python/cython and called by the quadpack fortran integrator
(also the code from http://people.sc.fsu.edu/~jburkardt/f_src/quadpack/quadpack.html was put here, which seems to be a more modern version than on
the http://www.netlib.org/quadpack/ one)

f90_vs_numpy_arrays_5.f90 is the version of the source code from http://people.sc.fsu.edu/~jburkardt/f_src/quadpack/quadpack.html rewritten
in to a modern fortran explicit interface style according to advices from http://fortran90.org/ (but gfortran still warns with debugging flags from
http://fortran90.org/src/faq.html#what-compiler-options-should-i-use-for-development)

fcython_qpck_2/ is for another C/Fortran wrapping style due to kleinert.pdf

testIntelMKL/ is for relative benchmarking of array and matrices operations, because fortran and python/cython integration is viewed as
very perspective especially on the array operations speedup http://technicaldiscovery.blogspot.ru/2011/06/speeding-up-python-numpy-cython-and.html
It is also interesing whether some speedup would take place in direct fortran usage in the way given in http://fortran90.org/src/best-practices.html#interfacing-with-python in comparison to f2py wrapping.
