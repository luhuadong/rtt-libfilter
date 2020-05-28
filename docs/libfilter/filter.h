//Constants for Hamming window

static const double HAMMING_MUL = 0.46164f;
static const double HAMMING_SUB = 0.53836f;

//Filters
void lowpass(int n, double *x, double *y, double alpha);
/*
Returns RC low-pass filter output samples, given input samples, and smoothing factor alpha.
http://en.wikipedia.org/wiki/Low-pass_filter
The algorithm simulates the effect of a low-pass filter on a series of digital samples.
	n	number of samples
	x	time series data
	y	filtered time series
	alpha	smoothing factor, 0 <= alpha <= 1; low alpha means heavy smoothing; alpha==1.0 means no smoothing
*/

void sma(int n, double *x, double *y, int m);
/* simple moving average
http://en.wikipedia.org/wiki/Moving_average
A simple moving average (SMA) is the unweighted mean of the previous m data points.
	n	number of samples
	x	time series data
	y	filtered time series
	m	window width, 1 < m << n
*/

void smm(int n, double *x, double *y, int m);
/* simple moving median
http://en.wikipedia.org/wiki/Moving_average
A simple moving median (SMM) is the unweighted median of the previous m data points.
From a statistical point of view, the moving average, when used to estimate the underlying trend in a time series, is susceptible to rare events such as rapid shocks or other anomalies. A more robust estimate of the trend is the simple moving median over n time points:
	n	number of samples
	x	time series data
	y	filtered time series
	m	window width, 1 < m << n
*/

void hanning3(int n, double *x, double *y);
/* 
simple Hanning 3-point filter
see: Weedon G 2003: Time -Series Analysis and Cyclostratigraphy, p.97
	n	number of samples
	x	time series data
	y	filtered time series
*/

//Windows for convolution
void win_get_rect(double * dest, unsigned int size);
/*
Gives rectangular (Dirichlet) 0-1 window
	dest    Window destination address
	size    width of the filter, number of points
Destination array is filled by ones. It must be size floats wide. 
http://en.wikipedia.org/wiki/Window_function#Rectangular_window
*/

void win_get_hamming(double * dest, unsigned int size);
/*
Gives Hamming 0-1 window
	dest    Window destination address
	size    width of the filter, number of points
Destination array is filled by coefficients from 0 to 1, which forms Hamming window. 
It must be size floats wide.
http://en.wikipedia.org/wiki/Window_function#Hamming_window
*/

void win_get_hann(double * dest, unsigned int size);
/*
Gives Hann 0-1 window
	dest    Window destination address
	size    width of the filter, number of points
Destination array is filled by coefficients from 0 to 1, which forms Hann window. 
It must be size floats wide.
http://en.wikipedia.org/wiki/Window_function#Hann_window
*/

void win_get_cos(double * dest, unsigned int size);
/*
Gives Cosine 0-1 window
	dest    Window destination address
	size    width of the filter, number of points
Destination array is filled by coefficients from 0 to 1, which forms Cosine window. 
It must be size floats wide.
http://en.wikipedia.org/wiki/Window_function#Cosine_window
*/

void win_get_lanczos(double * dest, unsigned int size);
/*
Gives Lanczos 0-1 window
	dest    Window destination address
	size    width of the filter, number of points
Destination array is filled by coefficients from 0 to 1, which forms Lanczos window. 
It must be size floats wide.
http://en.wikipedia.org/wiki/Window_function#Lanczos_window
*/

void win_get_bartlett_zero(double * dest, unsigned int size);
/*
Gives Bartlett 0-1 window window with zero ends
	dest    Window destination address
	size    width of the filter, number of points
Destination array is filled by coefficients from 0 to 1, which forms Bartlett window. 
It must be size floats wide.
http://en.wikipedia.org/wiki/Window_function#Triangular_windows
*/

void win_get_bartlett_nonzero(double * dest, unsigned int size);
/*
Gives Bartlett non-zero 0-1 window with non-zero ends
	dest    Window destination address
	size    width of the filter, number of points
Destination array is filled by coefficients from 0 to 1, which forms Bartlett window. 
It must be size floats wide.
http://en.wikipedia.org/wiki/Window_function#Triangular_windows
*/

void win_get_gaussian(double * dest, unsigned int size, double sigma);
/*
Gives  Gaussian 0-1 window
	dest    Window destination address
	size    width of the filter, number of points
Destination array is filled by coefficients from 0 to 1, which forms  Gaussian window. 
It must be size floats wide.
http://en.wikipedia.org/wiki/Window_function#Gaussian_windows
*/

void win_get_bartlett_hann(double * dest, unsigned int size, double sigma);
/*
Gives  Bartlett–Hann 0-1 window
	dest    Window destination address
	size    width of the filter, number of points
Destination array is filled by coefficients from 0 to 1, which forms  Bartlett–Hann window. 
It must be size floats wide.
http://en.wikipedia.org/wiki/Window_function#Bartlett.E2.80.93Hann_window
*/

void win_get_blackman(double * dest, unsigned int size, double alpha);
/*
Gives   Blackman 0-1 window
	dest    Window destination address
	size    width of the filter, number of points
	alpha	filter parameter, usually when unqualified == 0.16
Destination array is filled by coefficients from 0 to 1, which forms  Blackman window. 
It must be size floats wide.
http://en.wikipedia.org/wiki/Window_function#Blackman_windows
*/

void win_get_nuttall(double * dest, unsigned int size);
/*
Nuttall window, continuous first derivative
	dest    Window destination address
	size    width of the filter, number of points
	alpha	filter parameter, usually when unqualified == 0.16
Destination array is filled by coefficients from 0 to 1, which forms  Nuttall window. 
It must be size floats wide.
http://en.wikipedia.org/wiki/Window_function#Nuttall_window.2C_continuous_first_derivative
*/

void win_get_blackman_harris(double * dest, unsigned int size);
/*
Blackman–Harris window
A generalization of the Hamming family, produced by adding more shifted sinc functions, meant to minimize side-lobe levels.
	dest    Window destination address
	size    width of the filter, number of points
Destination array is filled by coefficients from 0 to 1, which forms  Blackman–Harris window. 
It must be size floats wide.
http://en.wikipedia.org/wiki/Window_function#Blackman.E2.80.93Harris_window
*/

void win_get_blackman_nuttall(double * dest, unsigned int size);
/*
Blackman–Nuttall window
	dest    Window destination address
	size    width of the filter, number of points
Destination array is filled by coefficients from 0 to 1, which forms  Blackman–Nuttall window. 
It must be size floats wide.
http://en.wikipedia.org/wiki/Window_function#Blackman.E2.80.93Nuttall_window
*/

void win_get_flattop(double * dest, unsigned int size);
/*
Flat top window
	dest    Window destination address
	size    width of the filter, number of points
Destination array is filled by coefficients from 0 to 1, which forms  Flat top window. 
It must be size floats wide.
http://en.wikipedia.org/wiki/Window_function#Flat_top_window
*/


