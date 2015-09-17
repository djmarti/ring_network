#ifndef _PARAMETERS_H
#define _PARAMETERS_H 1
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#define MAX_CHARS 200
/* const int MAX_CHARS = 200; I'd like this but doesn't work */

typedef struct Parameters Parameters;
typedef struct Sgmoid_pars Sgmoid_pars;
typedef struct Connectivity_pars Connectivity_pars;
typedef struct Input_pars Input_pars;

struct Sgmoid_pars {
    double max_rate;
    double beta;
    double thr;
};

struct Connectivity_pars {
    double J_E;
    double J_I;
    double m_E;
    double m_I;
};

struct NegFeedback_pars {
    /* Adaptation variables */
    double J_A; /* This var acts as a flag as well */
    double tau_A;
    /* Short term synaptic plasticity */
    bool facilitation_f; /* flag */
    bool depression_f;   /* flag */
    double U; /* Utilization factor */
    double tau_F; /* Facilitation timescale */
    double tau_D; /* Depression timescale */
};

struct Parameters {
    char output_rate_file[MAX_CHARS];
    char output_adaptation_file[MAX_CHARS];
    char output_depression_file[MAX_CHARS];
    char output_facilitation_file[MAX_CHARS];
    char extinputs_file[MAX_CHARS];
    char conf_file[MAX_CHARS];
    _Bool normalize;
    _Bool verbose;
    int N1;
    int N2;
    int N;
    int N_layers;
    double T;                   /* simulation time */
    double dt;
    struct Sgmoid_pars sgmoid;
    struct NegFeedback_pars negfbck;
    struct Connectivity_pars kernel;    /* Lateral connectivity */
    struct Connectivity_pars ff_kernel; /* Forward connectivity */
};

int handler(void *, const char *, const char *);

void init_parameters(struct Parameters *);

void show_parameters(void *);

#endif
