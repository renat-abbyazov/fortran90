#cdef public double x
#cdef public double f(x)

cdef extern:
    void c_quad(double *c_f, double *a, double *b, double *epsabs, 
                double *epsrel, double *result, double *abserr, 
                int *neval, int *ier)
    double c_f(double *x)
    
    
def fcython_quad(double c_f, double a, double b):
    cdef double epsabs, epsrel, result, abserr
    cdef int neval, ier
    c_quad(&c_f, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier)
    return result

