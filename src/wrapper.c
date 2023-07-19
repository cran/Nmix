#include <R.h>
void F77_SUB(rndstart)(void) { GetRNGstate(); }
void F77_SUB(rndend)(void) { PutRNGstate(); }
double F77_SUB(normrnd)(void) { return norm_rand(); }
double F77_SUB(unifrnd)(void) { return unif_rand(); }
double F77_SUB(exprnd)(void) { return exp_rand(); }
void F77_SUB(countdown)(int *k) {Rprintf("%d ", *k); if(*k==0) Rprintf("\n");}

