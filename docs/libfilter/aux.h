//Definitions
#ifndef max
#define max(a, b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
#define min(a, b) (((a) < (b)) ? (a) : (b))
#endif

#ifndef PI
	#define PI 3.141592654f
#endif
#define EPS 10e-16

//Functions
double *vec_alloc(int size);
// allocates memory for a vector

void vec_free(double *v);
// frees memory of the vector

void vec_write(double *v, int size);
// write out a vector as a row

void vec_read(double *v, int size);
// read a vector from ASCII stream

void vec_write_col(double *v, int size);
// write out a vector as a column

double sinc(double x);
/*
In digital signal processing and information theory, the normalized sinc function is commonly defined by
sinc(x) = sin(PI*x) / (PI*x)
It is qualified as normalized because its integral over all x is one. The Fourier transform of the normalized sinc function is the rectangular function with no scaling. This function is fundamental in the concept of reconstructing the original continuous bandlimited signal from uniformly spaced samples of that signal.
http://en.wikipedia.org/wiki/Sinc_function
*/

double median (int n, double *x);
/*
Find the median of X(1), ... , X(N), using as much of the quicksort
algorithm as is needed to isolate it.
N.B. On exit, the array X is partially ordered. If it is not desired, act on a copy of data vector.
Based on Alan J. Miller's median.f90 routine.
Copyright: (C) 2002 Michiel Jan Laurens de Hoon: mdehoon 'AT' gsc.riken.jp
Stolen from BioPython. See the license there.
*/

int kronecker_delta(int n); 
/*
Kronecker's delta, named after Leopold Kronecker (1823-1891), is a function of two variables, usually integers, which is 1 if they are equal, and 0 otherwise. In digital signal processing, the same concept is represented as a function on Z (the integers):
http://en.wikipedia.org/wiki/Kronecker_delta
*/

void convolve(double *x, double *y, int n, double *b, int width);

/*
convolution of the discrete data with filter
	x[j] is the input signal,
	y[j] is the output signal (convolution with the filter function),
	n is number of samples, or, equivalently, size of x and y
	b[i] are the filter coefficients, and
	width is the filter order – an nth-order filter has (n + 1) terms on the right-hand side; these are commonly referred to as taps.
*/

void imp_resp(double *x, double *b, double *h, int size);
/* Impulse response function
The impulse response h[j] can be calculated if we set x[n] = delta[n] in the above FIR relation, where delta[n] is the Kronecker delta impulse. The impulse response for an FIR filter then becomes the set of coefficients b[n].
http://en.wikipedia.org/wiki/Finite_impulse_response
	x[j] is the input signal,
	y[j] is the output signal,
	b[i] are the filter coefficients, and
	size is both number of samples, and the filter order
*/



