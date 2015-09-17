#include "utils.h"

double distance_on_torus(double x_1, double y_1, double x_2, double y_2)
{
    /* We assume the torus has dimensions 
     *   [-pi/2, +pi/2] x [-pi/2, +pi/2] */
    double t1, t2;
    t1 = fabs(x_1 - x_2);
    t2 = fabs(y_1 - y_2);
    t1 = fmin(t1, M_PI - t1);
    t2 = fmin(t2, M_PI - t2);
    return sqrt(t1 * t1 + t2 * t2);
}


/* This is the integrand function in the form required by 
 * the cubature library */
void connectivity_kernel_function(unsigned dim, const double *x,
                                  void *data_, unsigned func_dim __attribute__ ((unused)),
                                  double *retval)
{
    double *m = data_;          /* concentration parameter */
    double d = 0;
    for (unsigned i = 0; i < dim; ++i)
        d += x[i] * x[i];
    *retval = exp(*m * cos(2 * sqrt(d)));
}

void Fourier_coefficient_integrand(unsigned dim, const double *x,
                                   void *data_, unsigned func_dim __attribute__ ((unused)),
                                   double *retval)
{
    struct params *p = data_;   /* parameters */
    double d = 0;
    for (unsigned i = 0; i < dim; ++i)
        d += x[i] * x[i];
    *retval =
        exp(p->m * cos(2 * sqrt(d))) * cos(2 * (p->k1 * x[0] + p->k2 * x[1]));
}

double sum_on_torus(double m, int N1, int N2)
{
    double sum = 0;
    double dist;

    for (int i = 0; i < N1; i++) {
        for (int j = 0; j < N1; j++) {
            dist = distance_on_torus(0, 0, 
                                     i * M_PI / (double) N1, 
                                     j * M_PI / (double) N2);
            sum += exp(m * cos(2 * dist));
        }
    }
    return sum;
}

double integral_on_torus(double m)
{
    unsigned dim = 2;
    unsigned integrand_fdim = 1;

    double reqAbsError = 0;     /* required Absolute Error. 0: ignore */
    double reqRelError = 1e-5;  /* required Relative Error            */
    unsigned maxEval = 0;       /* no limit in number of interations  */

    double val, err;

    /* C99 */
    double xmin[dim];
    double xmax[dim];


    /* We assume the torus has dimensions 
     *  [-pi/2, +pi/2] x [-pi/2, +pi/2] 
     * Taking into account the symmetry in the distance function, 
     * the relevant * interval of integration for the connectivity 
     * kernel is [0,pi/2]  */

    for (unsigned i = 0; i < dim; ++i) {
        xmin[i] = 0.0;
        xmax[i] = M_PI / 2.0;
    }
    adapt_integrate(integrand_fdim, connectivity_kernel_function, &m,   /* fdata, additional info (pars,...) */
                    dim, xmin, xmax, maxEval, reqAbsError, reqRelError, &val, &err);    /* Output vars */
    return 4 * val;
}

double Fourier_integral_2D(double m, double k1, double k2)
{
    unsigned dim = 2;
    unsigned integrand_fdim = 1;

    double reqAbsError = 0;     /* required Absolute Error. 0: ignore */
    double reqRelError = 1e-5;  /* required Relative Error            */
    unsigned maxEval = 0;       /* no limit in number of interations  */

    double val, err;

    /* C99 */
    double xmin[dim];
    double xmax[dim];

    struct params p = { m, k1, k2 };

    /* We assume the torus has dimensions  */
      /* [-pi/2, +pi/2] x [-pi/2, +pi/2]  */
    /* Taking into account the symmetry in the distance function,  */
    /* the relevant * interval of integration for the connectivity  */
    /* kernel is [0,pi/2]  */


    for (unsigned i = 0; i < dim; ++i) {
        xmax[i] = M_PI / 2.0;
        xmin[i] = -xmax[i];
    }
    adapt_integrate(integrand_fdim, Fourier_coefficient_integrand, &p, dim, xmin, xmax, maxEval, reqAbsError, reqRelError, &val, &err); /* Output vars */
    return val;
}
