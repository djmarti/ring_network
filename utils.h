#ifndef _UTILS_H
#define _UTILS_H 1
#include <gsl/gsl_math.h>
#include "cubature-20101018/cubature.h"

/* Parameters to enter in Fourier_coefficient_integrand */
struct params {
    double m;
    double k1;
    double k2;
};

/* utils.c */
double distance_on_torus(double x_1, double y_1, double x_2, double y_2);
void connectivity_kernel_function(unsigned dim, const double *x, void *data_,
                                  unsigned func_dim, double *retval);
void Fourier_coefficient_integrand(unsigned dim, const double *x, void *data_,
                                   unsigned func_dim, double *retval);
double sum_on_torus(double m, int N1, int N2);
double integral_on_torus(double m);
double Fourier_integral_2D(double m, double k1, double k2);
#endif
