#include <R.h>
#include <Rmath.h>
#include <R_ext/Random.h>
#include <R_ext/Utils.h>

/* call before generating any random variates */
void F77_SUB(rndstart)(void) {GetRNGstate(); }

/* call after done generating all random variates */
void F77_SUB(rndend)(void) {PutRNGstate(); }

/* call to generate one uniform random variate */
double F77_SUB(rnd)(double *alpha1, double *alpha2)
{
    return runif(alpha1[0], alpha2[0]);
}

/* These are exposed in Utils.h and are misguidedly in the API */
void F77_SUB(qsort4)(double *v, int *indx, int *ii, int *jj)
{
  R_qsort_I(v, indx, *ii, *jj);
}