
/*
 * Program to obtain the coefficients of Omega
 */


#include "omega.h"

/* Program to obtain the Omega coefficients. */
/* Alberto Ramos <alberto@martin.ft.uam.es> */


int
main(int argc, char **argv) {

  mpf_t *factoriales;
  mpf_t pi;
  long n, m, k1, k2; 
  int c;

  int nmin = 0, nmax = 100, kmin= 0, kmax= 10, q = 1, PREC = 1000;
  double tau = 1.0;


  while ( (c = getopt(argc, argv, "hn:N:k:K:p:t:q:")) != -1) {
	 switch (c) {
	 case 'h':
		usage();
		exit(0);
		
	 case 'n':
		nmin = atoi(optarg);
		break;

	 case 'N':
		nmax = atoi(optarg);
		break;

	 case 'k':
		kmin = atoi(optarg);
		break;

	 case 'K':
		kmax = atoi(optarg);
		break;

	 case 'p':
		PREC = atoi(optarg);
		break;

	 case 't':
		tau = strtod(optarg, NULL);
		break;

	 case 'q':
		q = atoi(optarg);
		break;
		
	 default:
		usage();
		exit(0);

	 }
  }
  
  argc -= optind;
  argv += optind;


  mpf_set_default_prec( (unsigned long int) PREC);

  if ( (factoriales = (mpf_t *) malloc((NFAC+1)*sizeof(mpf_t))) == NULL)
	 perror("Allocating space for factoriales");
  readfac(factoriales);
  readpi(pi);

  for (n=nmin; n<=nmax; n++) {
	 for (m=n; m<=nmax; m++) {
		for (k1=kmin; k1<=kmax; k1++) {
		  for (k2=k1; k2<=kmax; k2++) {
			 omega(n, m, tau, q, k1, k2, factoriales, pi);
		  }
		}
	 }
  }
  

  
  printf("\n");

  exit(1);
}


void 
readfac (mpf_t *factoriales) {

  FILE *fichero;
  char *filename;
  int i;
  
  mpz_t intfac;


  filename = (char *) malloc(50*sizeof(char));
  mpz_init(intfac);

  mpf_init_set_ui(factoriales[0], 1);
  for (i=1; i<=NFAC; i++) {
	 sprintf(filename, "factoriales/%4.4d.dat", i);
	 fichero = fopen(filename, "r");

	 mpz_inp_str(intfac, fichero, 10);
	 /* printf("Ahi va: %d\n", i);
	  * mpz_out_str(stdout, 10, intfac);
	  * printf("\n\n");
	  */

	 mpf_init(factoriales[i]);
	 mpf_set_z(factoriales[i], intfac);

	 fclose(fichero);
  }
  
  
}

void 
omega(long int n, long int m, double tau, long int q, 
			  long int k1, long int k2, mpf_t *factoriales, mpf_t pi) {

  mpf_t aux, aux2, aux3, sqrf, acum, Ltau, z, zelev;
  int j;
  unsigned long int qsqr;

  /* 
	* Set n = min(n,m), and m = max(n,m)
	*/
  j = n;
  n = min(n,m);
  m = max(j,m);
	 
  /*
	* The sqrt(n!m!) pre-factor
	*/
  mpf_init(aux);
  mpf_init(sqrf);
  mpf_mul(aux, factoriales[n], factoriales[m]);
  mpf_sqrt(sqrf, aux);
  
  /*
	* Calculus of tau
	*/
  mpf_init(aux2);
  mpf_init_set_d(Ltau, tau);
  mpf_init(z);
  mpf_init(zelev);
  qsqr = (unsigned long int) q;
  
  
  mpf_set_d(aux, (double) pow((double) k1, (double) 2));
  mpf_mul(aux, aux, Ltau);
  
  mpf_set_d(aux2, (double) pow((double) k2, (double) 2));
  mpf_div(aux, aux, Ltau);

  mpf_add(aux, aux, aux2);
  mpf_div_ui(aux, aux, qsqr);

  mpf_mul(z, aux, pi);

  

  if ( ((m-n)%2) == 0)
	 mpf_pow_ui(zelev, z, (unsigned long int) (m-n)/2);
  else {
	 mpf_pow_ui(zelev, z, m-n);
	 mpf_sqrt(zelev, zelev);
  } 
	 

  /*  mpf_pow_ui(zelev, z, m-n);
	*mpf_pow_ui(z, z, 2);
	*/
  /*  mpf_out_str(stdout, 10, 20, z);
		printf("\n");*/
  /* 
	* The loop
	*/
  mpf_init(acum);
  mpf_init(aux3);
  mpf_set_ui(acum, (unsigned long int) 0); 

  for (j=0; j <= n; j++) {
	 mpf_mul(aux, factoriales[j], factoriales[n-j]);
	 mpf_mul(aux2, aux, factoriales[j+m-n]);

	 mpf_pow_ui(aux3, z, j);
	 mpf_div(aux, aux3, aux2);
	 if ((j%2) == 0) 
		mpf_set(aux2, aux);
	 else
		mpf_neg(aux2, aux);
	 
	 mpf_add(acum, acum, aux2);

  }
  


  mpf_mul(aux, acum, sqrf);
  mpf_mul(aux, aux, zelev);

  gmp_printf("%4d %4d %4d %4d %20.20Fe\n", n, m, k1, k2, aux);
  /*mpf_out_str(stdout, 10, 20, aux);*/

}

void 
readpi(mpf_t pi) {

  char *filename;
  FILE *file;

  if ( (filename = (char *) malloc(50*sizeof(char))) == NULL )
	 perror("readpi");

  sprintf(filename, "factoriales/pi.dat");
  file = fopen(filename, "r");

  mpf_init(pi);
  mpf_inp_str(pi, file, 10);
  fclose(file);
  

}

long int 
min(long int i, long int j) {

  if (i<j) 
	 return i;
  else
	 return j;

}

long int 
max(long int i, long int j) {

  if (i>j) 
	 return i;
  else
	 return j;

}

void
usage() {

  printf("\n");
  printf("Usage: \n");
  printf("  omega.out <options>\n");
  printf("\n");
  printf("Options:\n");
  printf("  -n nmin: Set min value of (n,m) to nmin (default 0)\n");
  printf("  -N nmax: Set max value of (n,m) to nmax (default 100)\n");	
  printf("  -k kmin: Set min value of (k1,k2) to kmin (default 0)\n");
  printf("  -K kmax: Set max value of (k1,k2) to kmax (default 10)\n");
  printf("  -p prec: Set precision to prec significant bits (default 1000)\n");
  printf("  -t tau: Set the aspect ratio l2/l1 to tau (default 1.0)\n");
  printf("  -q q: Set the value of the flux to q units (default 1)\n");

  printf("\n");

}

