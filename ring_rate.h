#ifndef _RING_RATE_H
#define _RING_RATE_H 1
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdarg.h>
#include <complex.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_vector.h>
#include "parameters.h"
#include "conf_file.h"
#include "colormap.h"
#include "circ_buffer.h"
#include "ext_inputs.h"
#include "utils.h"

struct Layer {
    double *rate;
    double *rate_old;
    double *adaptation; /* adaptation current */
    double *adaptation_old;
    double *depression;     /* depression variable */
    double *depression_old;
    double *facilitation;    /* facilitation variable */
    double *facilitation_old;
    double *kernel;             /* kernel indpendent of location (homogeneous networks) */
    double **kernel_mat;        /* kernel matrix, for inhomogeneous networks. Built from the kernel vector */
    struct Sgmoid_pars *spars;
};


struct Network {
    double I_firstlayer;
    double I_innerlayer;
    double *I_stim;
    double *I_stim_ahead;
    double I_max;
    struct Layer *layer;
    double *ff_kernel;
    
    /* Convenient vars */
    double max_evoked_rate;     /* this is the maximum rate achieved throughout the simulation */
    double min_evoked_rate;     /* and this is the minimum */
    double max_adaptation;
    double min_adaptation;
    double max_facilitation;
    double min_facilitation;
    double max_depression;
    double min_depression;
};

struct State {
    struct Parameters pconfig;  /* connectivity and 'neuronal' parameters  */
    struct pulse *pulse_list;   /* Linked list of pulses of external input */
    struct pulse_2D *pulse_list_2D;     /* Linked list of 2D pulses of external input */
    struct Network nw;
    struct circ_buffer *rate_images;
    struct circ_buffer *adaptation_images;
    struct circ_buffer *depression_images;
    struct circ_buffer *facilitation_images;
    struct circ_buffer input_image;
    bool is_2D;
};

/* Function Prototypes */

/* ring_rate.c */
/* Primary functions */
void read_network_parameters(int argc, char *argv[], struct State *S);
void read_input_parameters(struct State *S);
void generate_pulses(struct State *S, int N_categories, 
        double *prototype_location, double dispersion, double duration, 
        double intensity, double concentration_pulse);
void initialize_network(struct State *S);
void initialize_network_2D(struct State *S);
void reset_rates_and_inputs(struct State *S);
void initialize_rate_profile(struct State *S, double r_max, double m);
void initialize_circular_buffer(struct State *S);
void simulate_single_trajectory(struct State *S);
void simulate_single_trajectory_self_stop(struct State *S, double *time_formation, double *angle_location);
void simulate_single_trajectory_2D(struct State *S);
double get_angle_firing_highest(struct State *S, double threshold);
// void save_firing_activity(struct State *S, const char *output_fname);
void save_dynamic_variables(struct State *S);
void save_dynamic_variables_in_text(struct State *S);
void save_dynamic_variables_in_bitmaps(struct State *S);
void reset_all_bounds(struct State *S);
void determine_bounds(struct State *S);
void save_ppm(const struct circ_buffer *buff, const char *fname, double
        lobound, double hibound);
void save_snapshots_firing_activity_2D(struct State *S);
void save_inputs(struct State *S);
void free_network_variables(struct State *S);
void free_input_variables(struct State *S);

/* Secondary_functions */
void decide_if_network_is_2D(struct State *S);
double sigmoid(double x, const struct Sgmoid_pars *spars);
double sigmoid_inv(double r, const struct Sgmoid_pars *spars);
double tune_I_firstlayer(double r, const struct Parameters *pconfig);
double tune_I_innerlayer(double r, const struct Parameters *pconfig);
void set_all_values(double *v, int N, double x);
void initialize_kernel_vector(double *krnl, int N,
                              const struct Connectivity_pars *cp);
void initialize_kernel_vector_2D(double *krnl, int N1, int N2,
                                 const struct Connectivity_pars *cp);
void initialize_modulation (double *v, int N, double lambda, double concentration);
void initialize_kernel_matrix(double **krnl_mat, const double *krnl, int N,
                              double lambda, double concentration);
void initialize_kernel_matrix_with_noise(double **krnl_mat, const double *krnl, int N,
                              double epsilon);
void initialize_sigmoid_pars(struct Sgmoid_pars *, int N,
                             const struct Sgmoid_pars *sp);
void show_fourier_coeffs(const struct Parameters *pconfig);
void show_fourier_coeffs_2D(const struct Parameters *pconfig);
void initialize_rng(void);
struct pulse *inputs_update(double t, struct State *S);
void white_noise_input_update(struct State *S);
struct pulse_2D *inputs_update_2D(double t, struct State *S);
void network_update_Euler(struct State *S);
void network_update_RK4(struct State *S);
void advance_inputs_one_step(struct State *S);
void network_update_2D(struct State *S);
void update_bounds_if_necessary(double v, double *min, double *max);
void set_filename(char *filename, const char *basename, int i, const char *ext);
void save(const struct circ_buffer *buff, const char *fname);
void save_population_vector(const struct circ_buffer *buff, const char *fname);
void save_pgm(const struct circ_buffer *buff, const char *fname, double lo,
              double hi, int N);
void save_circ_buffer_ppm(const struct circ_buffer *buff, const char *fname,
                          double lo, double hi);
void save_2D(const struct circ_buffer *buff, const char *fname, int N1, int N2);
void save_2D_ppm(const struct circ_buffer *buff, const char *fname, double lo,
                 double hi, int N1, int N2);
void replace_substring(char *p, const char *str, const char *orig, const char *rep);
#endif
