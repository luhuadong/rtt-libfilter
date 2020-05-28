/* 
Discrete signal filtering functions.
++pac (Peter A. Cejchan) 2010 <tyapca7@gmail.com>
GPL3
Some inspiration taken from CIV Toolkit, and R:signal package.
*/

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include "aux.h"
#include "filter.h"



// filters

void
lowpass(int n, double *x, double *y, double alpha)
{
	int i;
	y[0] = x[0];
	for(i=1; i<n; i++) 
		y[i] = alpha * x[i] + (1.0-alpha) * y[i-1];
   	return;
}


void
sma(int n, double *x, double *y, int m)
{
	int i;
	double sum;
	for(i=0; i<m; i++) {
		y[i] = x[i];
		sum += x[i];
	}
	y[m] = sum/m;
	for(i=m+1; i<n; i++) {
		y[i] = y[i-1] - (x[i]-m)/m + x[i]/m;
	}
   	return;
}


void
smm(int n, double *x, double *y, int m)
{
	int i, j;
	double *w, sum;
	// allocate memory for window
	w=vec_alloc(m);
	// first m points are unchanged
	for(i=0; i<m; i++) {
		y[i] = x[i];
		sum += x[i];
	}
	// following points are the median
	for(i=m; i<n; i++) {
		// window
		for(j=i-m+1; j<=i; j++)
			w[j]=x[i];
		y[i] = median(m, w);
		
	}
   	return;
}

void
hanning3(int n, double *x, double *y)
{
  int i;

	//weights
	double w1=0.25;
	double w2=0.50;
	double w3=0.25;

	// first 2 points are unchanged

	// following points are convoluted by the weights w1, w2, w3
	for(i=2; i<n; i++) {
		// sliding window
		y[i] = (x[i-1]*w1 + x[i]*w2 +  x[i-1]*w3)/3;
	}
   	return;
}

// filter windows
// n --> i
// N --> size
void 
win_get_hamming(double *dest, unsigned int size) 
{
    unsigned int i = 0;
    
    size--;
    while (i <= size) {
        *dest++ = HAMMING_SUB - HAMMING_MUL * cos((2 * PI * i++) / size);
    } // while
}

void 
win_get_hann(double * dest, unsigned int size) 
{
    unsigned int i = 0;
    
    size--;
    while (i <= size) {
        *dest++ = 0.5*(1.0 - cos((2 * PI * i++) / size));
    } // while
}

void 
win_get_cos(double * dest, unsigned int size) 
{
    unsigned int i = 0;
    
    size--;
    while (i <= size) {
        *dest++ = sin((PI * i++) / size);
    } // while
}
 
void 
win_get_lanczos(double *dest, unsigned int size) 
{
    unsigned int i = 0;
    
    size--;
    while (i <= size) {
        *dest++ = sinc(((2.0 * i++) / (double)size) - 1.0);
    } // while
}

void 
win_get_bartlett_zero(double * dest, unsigned int size) 
{
    unsigned int i = 0;
    
    while (i < size) {
        *dest++ = 2.0/(size-1.0) * ((size-1.0)/2.0 - fabs(i++ -	(size-1.0)/2.0));
    } // while
}

void 
win_get_bartlett_nonzero(double * dest, unsigned int size) 
{
    unsigned int i = 0;
    
    while (i < size) {
        *dest++ = 2.0/(size) * ((size)/2.0 - fabs(i++ -	(size-1.0)/2.0));
    } // while
}


void 
win_get_gaussian(double * dest, unsigned int size, double sigma) 
{
    unsigned int i = 0;
    if(sigma > 0.5) sigma = 0.5; if (sigma == 0.0) sigma = EPS;	// adjust 0 < sigma <= 0.5
    while (i < size) {
        *dest++ = exp(0.5*pow(((i++ - (size-1)/2))/(sigma*(size-1.0)/2), 2.0));
    } // while
}

void 
win_get_bartlett_hann(double * dest, unsigned int size, double sigma) 
{
  unsigned int i = 0;
  double a0 = 0.62;
  double a1 = 0.48;
  double a2 = 0.38;

  while (i < size) {
    *dest++ = a0 - a1 * fabs(i/(size-1.0) - 0.5) - a2 * cos(2*PI*i/(size-1.0));
    i++;
  } // while
}


void 
win_get_blackman(double * dest, unsigned int size, double alpha) 
{
  unsigned int i = 0;
  double a0;
  double a1;
  double a2;

  if(alpha==0) alpha = 0.16;

  a0=1.0 - alpha/2.0;
  a1=0.5;
  a2=alpha/2.0;

  while (i < size) {
    *dest++ = a0 - a1 * cos(2*PI*i/(size-1.0)) + a2 * cos(4*PI*i/(size-1.0));
    i++;
  } // while
}

void 
win_get_nuttall(double * dest, unsigned int size) 
{
  unsigned int i = 0;
  double a0=0.355768;
  double a1=0.487396;
  double a2=0.144232;
  double a3=0.012604;

  while (i < size) {
    *dest++ = a0 - a1 * cos(2*PI*i/(size-1.0)) + a2 * cos(4*PI*i/(size-1.0)) - a3 * cos(6*PI*i/(size-1.0));
    i++;
  } // while
}
 
void 
win_get_blackman_harris(double * dest, unsigned int size) 
{
  unsigned int i = 0;
  double a0=0.355768;
  double a1=0.487396;
  double a2=0.144232;
  double a3=0.012604;

  while (i < size) {
    *dest++ = a0 - a1 * cos(2*PI*i/(size-1.0)) + a2 * cos(4*PI*i/(size-1.0)) - a3 * cos(6*PI*i/(size-1.0));
    i++;
  } // while
}

void 
win_get_blackman_nuttall(double * dest, unsigned int size)
{
  unsigned int i = 0;
  double a0=0.3635819;
  double a1=0.4891775;
  double a2=0.1365995;
  double a3=0.0106411;

  while (i < size) {
    *dest++ = a0 - a1 * cos(2*PI*i/(size-1.0)) + a2 * cos(4*PI*i/(size-1.0)) - a3 * cos(6*PI*i/(size-1.0));
    i++;
  } // while
}

void 
win_get_flattop(double * dest, unsigned int size)
{
  unsigned int i = 0;
  double a0=1.0;
  double a1=1.93;
  double a2=1.29;
  double a3=0.388;
  double a4=0.032;

  while (i < size) {
    *dest++ = a0 - a1 * cos(2*PI*i/(size-1.0)) + a2 * cos(4*PI*i/(size-1.0)) - a3 * cos(6*PI*i/(size-1.0)) + a4 * cos(8*PI*i/(size-1.0));
    i++;
  } // while
}

void 
win_get_rect(double *dest, unsigned int size) 
{
// same as boxcar

    while (size) {
        size -= 1;
        *dest++ = 1.0f;
    } // while
}

void 
win_get_boxcar(double *dest, unsigned int size) 
{
  unsigned int i;
//  if(! (size > 0)) exits("size must be an integer > 0");

  for(i=0; i<size; i++) 
    dest[i] = 1.0f;
    
  return;
}

void 
win_get_hanning(double *dest, unsigned int size) 
{
/*
hanning  <- function(n)  {

  if (! (length(n) == 1 && (n == round(n)) && (n > 0)))
    stop("hanning: n has to be an integer > 0")

  if (n == 1)
    c = 1
  else {
    n = n - 1
    c = 0.5 - 0.5 * cos(2 * pi * (0:n) / n)
  }
  c
} */
  unsigned int i;
//  if(! (size > 0)) exits("hanning: 'size' has to be an integer > 0");
  if (size == 1) dest[0] = 1.0;
  else {
    for(i=0; i<=size; i++) // sic
      dest[i] = 0.5 - 0.5 * cos(2 * PI * i / size);
  } 
  return;
} 

void 
win_get_kaiser(double *dest, unsigned int size, double beta) 
{
/*
kaiser <- function(n, beta)  { 
  
  if ( !(length(n) == 1 && (n == round(n)) && (n > 0))) 
    stop("kaiser:  n has to be a positive integer")

  if ( !(length(beta) == 1 && (beta == as.real(beta))))
    stop("kaiser:  beta has to be a real scalar")
  
  if (n == 1)
    w = 1
  else {
    m = n - 1
    k = 0:m
    k = 2 * beta / m * sqrt (k * (m - k))
    w = besselI(k, 0) / besselI(beta, 0)
  }
  w
}
  unsigned int i, m;
  if(! (size > 0)) exits("kaiser: 'size' has to be an integer > 0");
  if (size == 1)
    dest[0] = 1;
  else {
    m = n - 1;
##### TO BE CONTINUED...
  }
*/
}

