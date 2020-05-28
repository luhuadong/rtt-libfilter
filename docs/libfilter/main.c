/*
main.c
apply a digital filter to the time-series
*/

#include <stdio.h>
#include <stdlib.h>
#include "aux.h"
#include "filter.h"

// convolution filters
#define HAMMING 0
#define DIRICHLET 1
#define HANN 2
#define COS 3
#define LANCZOS 4
#define BARTLETT0 5
#define BARTLETTN 6
#define GAUSSIAN 7
#define BARTLETT_HANN 8
#define BLACKMAN 9
#define NUTTALL 10
#define BLACKMAN_HARRIS 11
#define BLACKMAN_NUTTALL 12
#define FLATTOP 13

// non-convolution filters start at 91
#define LOWPASS 91
#define SMA 92
#define SMM 93
#define HANNING3 94

double 
*x, 		// data
*y, 		// results
*b, 		// window
alpha;		// filter coefficient

int 
  filter,	// filter type
  n, 		// number of samples == data points in x, y
  width;	// width of the filter (window)

void
usage(void){
  printf("filter filter_type filter_width [alpha] < data");
}


int
main(int argc, char *argv[])
{
  /* read-in command-line arguments */
  if ((argc>4) || (argc<2)) {usage(); return(-1);};
  if (argc== 2) {
    //default settings:
    n = atoi(argv[1]);
    filter = HAMMING;
    width = 5;
  } else if (argc== 3) {
    n = atoi(argv[1]);
    filter = atoi(argv[2]);
    width = atoi(argv[3]);
    alpha = 0.0;

  } else {
    n = atoi(argv[1]);
    filter = atoi(argv[2]);
    width = atoi(argv[3]);
    alpha = atof(argv[4]);
  }



  //read length
//  scanf("%d", &n);
  // allocate memory for the vectors
  x= vec_alloc(n);
  y= vec_alloc(n);
  b= vec_alloc(width);

  // read a vector from stdin ASCII stream (length x[1] x[2]...)
  vec_read(x, n);
//  vec_write(x, n);


  // if window not defined
  if (filter > 90) {
      switch (filter) {
        case SMA:
          sma(n, x, y, width);

        case LOWPASS:
        default:
          // apply lowpass filter
          lowpass(n, x, y, alpha);
	  break;
      }

  } else {
      switch (filter) {
        case LANCZOS:
          win_get_lanczos(b, width);
	  break;
        case HAMMING:
        default:
          // create a Hamming window
          win_get_hamming(b, width);
//  vec_write(b, width);
	  break;
      }
      // convolve the data with the filter
      convolve(x, y, n, b, width);
   }

  // write out a filtered vector
  vec_write_col(y, n);

  // free memory
  vec_free(x);
  vec_free(y);
  vec_free(b);
  return(0);
} //main

