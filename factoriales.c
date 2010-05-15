
/*
 * Program t obtain the coefficients of Omega
 */


#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <gmp.h>

void omega(long int n, long int m, double tau, long int q, 
			  long int k1, long int k2);



int
main() {

  mpf_t factorial;
  mpz_t intfac;

  FILE *fichero;
  char *filename;
  int i;

  mpf_set_default_prec( (unsigned long int) 100000);

  mpf_init(factorial);
  mpz_init(intfac);

  
  mpz_fac_ui(intfac, 10000);
  mpf_set_z(factorial, intfac);


  filename = (char *) malloc(50);
  for (i=1; i<=1000; i++) {
	 sprintf(filename, "factoriales/%4.4d.dat", i);
	 fichero = fopen(filename, "w");
	 printf("%s", filename);

	 mpz_fac_ui(intfac, i);
	 mpz_out_str(fichero, 10, intfac);
	 printf("->Hecho\n");
  }


  /*
  mpz_out_str(stdout, 10, intfac);
  printf("\n\n\n");
  mpf_out_str(stdout, 10, 0, factorial);
  printf("\n");
  */
  exit(1);
}

void 
omega(long int n, long int m, double tau, long int q, 
			  long int k1, long int k2) {

  

}
