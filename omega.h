
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <gmp.h>


/* Definition of constants. */

enum { 
  NFAC=1000
};


/* Function "headers" */

void omega(long int n, long int m, double tau, long int q, 
			  long int k1, long int k2, mpf_t *factoriales, mpf_t pi);

void readfac (mpf_t *factoriales);
void readpi(mpf_t pi);

long int min(long int i, long int j);
long int max(long int i, long int j);

void usage();

