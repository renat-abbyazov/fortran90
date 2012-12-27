#include <stdio.h>
#include <math.h>

extern void dqags_(double f(double *x), double *a, double *b, double *epsabs,
		   double *epsrel, double *result, double *abserr, int *neval,
		   int *ier, int *limit, int *lenw, int *last, int iwork[],
		   double work[]);

double f2(double *x) {
return exp(*x);
}

#define LIMIT 5
#define LENW 4*LIMIT

int main(int argc, char **argv)
{
  int neval, ier, limit = LIMIT, lenw = LENW, last, iwork[LIMIT];
  double a = 0., b = 1.e0, epsabs = 1.e-10, epsrel = 0., result, abserr, work[LENW];
  double exact = exp(1.e0) - 1;
  
  dqags_(f2, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval,
	 &ier, &limit, &lenw, &last, iwork, work);
  if (ier > 0)
    {
      fprintf(stderr, "DQAGS returned IER = %d\n", ier);
    }
  printf ("C result = %20.15lf\n", result);
  printf ("estimated error = %20.15lf\n", abserr);
  printf ("actual error = %20.15lf\n",abs(exact - result));
  printf ("computed using %d function evaluations \n\n",neval);
  
  return 0;
}

