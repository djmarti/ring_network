#include <math.h>
#include <stdio.h>
#include "cubature.h"

/* compile: gcc -o my_test curbature.o -lm */

int count = 0;

void my_integrand(unsigned dim, const double *x, void *data_,
	    unsigned func_dim, double *retval);

int main()
{
    double reqAbsError, reqRelError;
    double val, err;
    unsigned i, maxEval;
    unsigned dim, integrand_fdim;

    dim = 2;
    reqAbsError = 0;    /* required Absolute Error. 0: ignore */
    reqRelError = 1e-5; /* required Relative Error            */ 
    maxEval = 0;        /* no limit in number of interations  */
    integrand_fdim = 1;

    /* C99 */
    double xmin[dim];
    double xmax[dim];
    for (i = 0; i < dim; ++i) {
        xmin[i] = 0;
        xmax[i] = 1;
    }

    printf("%u-dim integral, tolerance = %g\n", dim, reqRelError);
    adapt_integrate(integrand_fdim, 
            my_integrand,
            NULL,  /* fdata, additional info (pars,...) */
            dim, xmin, xmax, 
            maxEval, 
            reqAbsError, 
            reqRelError, 
            &val, &err); /* Output vars */

    printf("integral = %g, est err = %g, true err = %g\n", 
            val, err, fabs(val - 0.484999));
    printf("#evals = %d\n", count);

    return 0;
}


void my_integrand(unsigned dim, const double *x, void *data_,
	    unsigned func_dim, double *retval)
{
     double tmp = 0;
     count++;
     unsigned i;
     (void) data_; /* not used */
     if (func_dim > 1)
         printf("Only functions on R\n");
     for (i = 0; i < dim; ++i)
         tmp += x[i] * x[i];
     *retval = exp(-sqrt(tmp));
}

