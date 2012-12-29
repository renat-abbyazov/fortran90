cdef public double x
cdef public double f(x)

cdef extern:
    void c_quad(double *f, double *a, double *b, double *epsabs, 
                double *epsrel, double *result, double *abserr, 
                double *neval, double *ier)

def fcython_quad(double f(double x), double a, double b):
    c_quad(&f, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier)
    return result

