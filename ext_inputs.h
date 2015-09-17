#ifndef _EXT_INPUTS_H
#define _EXT_INPUTS_H 1
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include "utils.h"

/* A list containing the properties of
 * the external stimulation pulses */
struct pulse {
    double time_onset;
    double duration;
    double theta;
    double Imax;
    double concentration;
    double *precomputed_value;
    _Bool switched_on;          /* this switch avoids redundant checks */
    struct pulse *next;         /* next in the list */
};

/* 2D version. Pulses are assumed to be isotropic */
struct pulse_2D {
    double time_onset;
    double duration;
    double theta1;
    double theta2;
    double Imax;
    double concentration;
    double *precomputed_value;
    _Bool switched_on;          /* this switch avoids redundant checks */
    struct pulse_2D *next;      /* next in the list */
};


/* ext_inputs.c */
struct pulse *new_pulse(int N, double t, double d, double q, double input,
                        double m);
struct pulse_2D *new_pulse_2D(int N1, int N2, double t, double d, double q1,
                              double q2, double input, double m);
double *precompute_vector_pulses(int N, struct pulse *plse);
double *precompute_vector_pulses_2D(int N1, int N2, struct pulse_2D *plse);
struct pulse *append_pulse(struct pulse *plse, struct pulse *newplse);
struct pulse_2D *append_pulse_2D(struct pulse_2D *plse,
                                 struct pulse_2D *newplse);
void print_pulses(struct pulse *plse);
void print_pulses_2D(struct pulse_2D *plse);
double max_amplitude_pulse(struct pulse *plse);
double max_amplitude_pulse_2D(struct pulse_2D *plse);
void free_extinputs(struct pulse *plse);
void free_extinputs_2D(struct pulse_2D *plse);
#endif
