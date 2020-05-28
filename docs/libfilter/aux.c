#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include "aux.h"

// auxilliary functions

double
sinc(double x)
{
  double y;
  if(x<EPS) x=EPS;

  y=sin(PI*x) / (PI*x);
  return y;
}

double *
vec_alloc(int size)
{
	int i;
	double *v;
	v=(double *)malloc(size*sizeof(double));
	for(i=0;i<size;i++)
		v[i]=0.0;
	return v;
}

void
vec_free(double *v)
{
	free(v);
	return;
}

void 
vec_read(double *v, int size) 
{
  int i;
  for (i=0; i<size; i++) {
    scanf("%lf", &v[i]);
  }
}

void 
vec_write(double *v, int size)
{
	int i;

	for(i=0; i<size; i++){
		printf("%f ", v[i]);
	}
	printf("\n");
	return;
}

void 
vec_write_col(double *v, int size)
{
	int i;

	for(i=0; i<size; i++){
		printf("%f\n", v[i]);
	}
	printf("\n");
	return;
}

double 
median (int n, double *x)
{
  int i, j;
  int nr = n / 2;
  int nl = nr - 1;
  int even = 0;
  /* hi & lo are position limits encompassing the median. */
  int lo = 0;
  int hi = n-1;

  if (n==2*nr) even = 1;
  if (n<3)
  { if (n<1) return 0.;
    if (n == 1) return x[0];
    return 0.5*(x[0]+x[1]);
  }

  /* Find median of 1st, middle & last values. */
  do
  { int loop;
    int mid = (lo + hi)/2;
    double result = x[mid];
    double xlo = x[lo];
    double xhi = x[hi];
    if (xhi<xlo)
    { double temp = xlo;
      xlo = xhi;
      xhi = temp;
    }
    if (result>xhi) result = xhi;
    else if (result<xlo) result = xlo;
    /* The basic quicksort algorithm to move all values <= the sort key (XMED)
     * to the left-hand end, and all higher values to the other end.
     */
    i = lo;
    j = hi;
    do
    { while (x[i]<result) i++;
      while (x[j]>result) j--;
      loop = 0;
      if (i<j)
      { double temp = x[i];
        x[i] = x[j];
        x[j] = temp;
        i++;
        j--;
        if (i<=j) loop = 1;
      }
    } while (loop); /* Decide which half the median is in. */

    if (even)
    { if (j==nl && i==nr)
        /* Special case, n even, j = n/2 & i = j + 1, so the median is
         * between the two halves of the series.   Find max. of the first
         * half & min. of the second half, then average.
         */
        { int k;
          double xmax = x[0];
          double xmin = x[n-1];
          for (k = lo; k <= j; k++) xmax = max(xmax,x[k]);
          for (k = i; k <= hi; k++) xmin = min(xmin,x[k]);
          return 0.5*(xmin + xmax);
        }
      if (j<nl) lo = i;
      if (i>nr) hi = j;
      if (i==j)
      { if (i==nl) lo = nl;
        if (j==nr) hi = nr;
      }
    }
    else
    { if (j<nr) lo = i;
      if (i>nr) hi = j;
      /* Test whether median has been isolated. */
      if (i==j && i==nr) return result;
    }
  }
  while (lo<hi-1);

  if (even) return (0.5*(x[nl]+x[nr]));
  if (x[lo]>x[hi])
  { double temp = x[lo];
    x[lo] = x[hi];
    x[hi] = temp;
  }
  return x[nr];
}

int
kronecker_delta(int n) 
{
  int d;
  if(n==0) d=1; else d=0;
  return d;
}

void
convolve(double *x, double *y, int n, double *b, int width)
{
  int i, j;
  for(j=0; j<n; j++) y[j] = 0.0; 

  for(j=width-1; j<n; j++) {
    y[j] = 0.0; 
    for(i=0; i<width; i++) {
	y[j] += b[i]*x[j-i];
    }
  }
}

void
imp_resp(double *x, double *b, double *h, int size) 
{
  int i, j;
  for(j=0; j<size; j++) {
    h[j] = 0.0; 
    for(i=0; i<size; i++) {
      h[j]+=b[i]*kronecker_delta(j-i);
    }
  }
}


