#include "ring_rate.h"

/* External (global) variables are not necessarily bad */

const int TIME_WINDOW_WIDTH = 500;     /* in time units */
const int N_SNAPSHOTS = 20;     /* number of frames in a 2D movie (per trial) */
const double R = 0.1;           /* firing rate stationary uniform state */
const int MAXVAL = 255;         /* maximum value for color codes */

gsl_rng *Rng;                   /* Global random number generator */
const gsl_rng_type *Rng_T;      /* Type of rng */

/* Primary functions (those appearing in main) */

void read_network_parameters(int argc, char *argv[], struct State *S)
{
    struct Parameters *pc = &S->pconfig;
    init_parameters(pc);
    /* Check if the user specified a config file in the command line. If true, update
     * the config_file member. */
    process_config_in_options(pc, argc, argv);

    /* Parse the config file */
    if (conf_parse(pc->conf_file, handler, pc) < 0)
        printf("'%s' cannot be loaded. Using default values and options.",
               pc->conf_file);

    /* Modify the parameters given in the config file with the options
     * specified in the command line */
    process_options(pc, argc, argv);
    if (pc->verbose)
        show_parameters(pc);

    decide_if_network_is_2D(S);
}

void decide_if_network_is_2D(struct State *S)
{
    if (S->pconfig.N2 != 1 && S->pconfig.N1 != 1) {
        S->is_2D = true;
    }
    /* if (S->pconfig.N2 == 1 || S->pconfig.N1 == 1) { */
        /* printf ("The network is actually 1D. "); */
        /* printf ("Use 'simulate_trajectory' instead, please.\n"); */
    /* } else { */
        /* S->is_2D = true; */
    /* } */
}

void read_input_parameters(struct State *S)
{
    struct Parameters *pc = &S->pconfig;

    if (!S->is_2D) {
        if (extinputs_parse(pc->extinputs_file, &S->pulse_list, pc->N) != 0) {
            printf("'%s' for a 1D network cannot be parsed.",
                   pc->extinputs_file);
            printf("No external input applied.\n");
        }
        if (pc->verbose) {
            print_pulses(S->pulse_list);
        }
    } else {
        if (extinputs_parse_2D
            (pc->extinputs_file, &S->pulse_list_2D, pc->N1, pc->N2) != 0) {
            printf("'%s' for a 2D network cannot be parsed.",
                   pc->extinputs_file);
            printf("No external input applied.\n");
        }
        if (pc->verbose) {
            print_pulses_2D(S->pulse_list_2D);
        }
    }
}

void generate_pulses(struct State *S, int N_categories, 
        double *prototype_location, double dispersion, double duration, 
        double intensity, double concentration_pulse)
{
    S->pulse_list = NULL;
    S->pulse_list_2D = NULL;
    struct Parameters *pc = &S->pconfig;
    
    /* inter-pulse separation */
    /* d_0 = 1e-2 when d = 0.1 (10%). To keep the rate nu constant we
                          should rescale d_0 as well */
    double sep0 = 1e-2; 
    double nu = 0.1 / (0.1 + sep0);
    double sep = duration / nu - duration;
    int total_number_pulses = (int) (pc->T / (duration + sep)); /* this covers the whole range */

    struct pulse **pulse_list = &S->pulse_list;

    /* Sample total_number_pulses items from the container */
    /* This works when all categories are equiprobable */
    int container[N_categories];
    int chosen_categories[total_number_pulses];
    for (int i = 0; i < N_categories; i++)
        container[i] = i;
    gsl_ran_sample(Rng, chosen_categories, total_number_pulses, 
            container, N_categories, sizeof(int));

    double mu;
    double time;
    /* The following works only when N_categories=2 */
    /* if (N_categories != 2) { */
        /* printf("The code assumes there are 2 categories only. Recode, please.\n"); */
        /* exit(1); */
    /* } */
    /* unsigned i_cat; */
    /* double p = 2./3.; [> probability of choosing first category <] */
    for (int i = 0; i < total_number_pulses; i++) {
        /* i_cat = 1 - gsl_ran_bernoulli(Rng, p); */
        mu = (180 / M_PI) * (prototype_location[chosen_categories[i]] + gsl_ran_gaussian(Rng, 1.0/sqrt(2.0*dispersion)));
        /* printf("%5.2f ", mu); */
        time = 0.01 + i * (duration + sep);
        *pulse_list = append_pulse(*pulse_list, 
                new_pulse(pc->N, time, duration, -mu, intensity, concentration_pulse));
    }
}

void initialize_network(struct State *S)
{
    /* Make the code a bit more readable */
    struct Parameters *pc = &S->pconfig;
    int N = pc->N;
    int N_layers = pc->N_layers;

    /* allocate space */
    /* ============== */
    S->nw.layer = (struct Layer *) malloc(N_layers * sizeof(struct Layer));
    for (int i = 0; i < N_layers; i++) {
        S->nw.layer[i].spars = malloc(N * sizeof(struct Sgmoid_pars));
        S->nw.layer[i].rate = malloc(N * sizeof(double));
        S->nw.layer[i].rate_old = malloc(N * sizeof(double));
        S->nw.layer[i].adaptation = malloc(N * sizeof(double));
        S->nw.layer[i].adaptation_old = malloc(N * sizeof(double));
        S->nw.layer[i].depression = malloc(N * sizeof(double));
        S->nw.layer[i].depression_old = malloc(N * sizeof(double));
        S->nw.layer[i].facilitation = malloc(N * sizeof(double));
        S->nw.layer[i].facilitation_old = malloc(N * sizeof(double));
        S->nw.layer[i].kernel = malloc(N * sizeof(double));
        S->nw.layer[i].kernel_mat = malloc(N * sizeof(double *));
        for (int j = 0; j < N; j++)
            S->nw.layer[i].kernel_mat[j] = malloc(N * sizeof(double));
    } 
    S->nw.I_firstlayer = tune_I_firstlayer(R, pc);
    S->nw.I_innerlayer = tune_I_innerlayer(R, pc);
    S->nw.ff_kernel = malloc(N * sizeof(double));
    S->nw.I_stim = malloc(N * sizeof(double));
    S->nw.I_stim_ahead = malloc(N * sizeof(double));

    /* and then, assign */
    struct Sgmoid_pars *sp = &pc->sgmoid;
    struct Connectivity_pars *cp = &pc->kernel;
    struct Connectivity_pars *cpff = &pc->ff_kernel;

    /* inhomogeneity parameters */
    /* double inhomog_lambda = 0.0;  */
    /* double inhomog_kappa = 0.0; */
    double u_steady;
    if (pc->negfbck.facilitation_f) {
        u_steady = pc->negfbck.U * ( 1 + pc->negfbck.tau_F * R );
        u_steady = u_steady / (1.0 + pc->negfbck.U * R * pc->negfbck.tau_F );
    } else {
        u_steady = 1.0;
    }
        
    float epsilon = 0.0;
    for (int i = 0; i < N_layers; i++) {
        initialize_sigmoid_pars(S->nw.layer[i].spars, N, sp);
        initialize_kernel_vector(S->nw.layer[i].kernel, N, cp);
        /* initialize_kernel_matrix(S->nw.layer[i].kernel_mat, S->nw.layer[i].kernel, N, inhomog_lambda, inhomog_kappa);  */
        initialize_kernel_matrix_with_noise(S->nw.layer[i].kernel_mat, S->nw.layer[i].kernel, N, epsilon); 
        set_all_values(S->nw.layer[i].rate, N, R);
        set_all_values(S->nw.layer[i].adaptation, N, 0);
        set_all_values(S->nw.layer[i].depression, N, 1);
        set_all_values(S->nw.layer[i].facilitation, N, u_steady);
        /* initialize_modulation(S->nw.layer[i].facilitation, N, 0.2, 10); */
        /* the _old counterparts are assigned at update_network */
    }
    set_all_values(S->nw.I_stim, N, 0.0);
    set_all_values(S->nw.I_stim_ahead, N, 0.0);

    /* Initialize lists of inputs */
    S->pulse_list = NULL;
    S->pulse_list_2D = NULL;

    if (N_layers > 1)
        initialize_kernel_vector(S->nw.ff_kernel, N, cpff);

    reset_all_bounds(S);

    if (pc->verbose) {
        show_fourier_coeffs(pc);
    }
}

void reset_rates_and_inputs(struct State *S)
{
    int N = S->pconfig.N;
    struct Parameters *pc = &S->pconfig;

    double U0;
    if (pc->negfbck.facilitation_f) {
        U0 = S->pconfig.negfbck.U;
    } else {
        U0 = 1.0;
    }
    for (int i = 0; i < S->pconfig.N_layers; i++) {
        set_all_values(S->nw.layer[i].rate, N, R);
        set_all_values(S->nw.layer[i].rate_old, N, R);
        set_all_values(S->nw.layer[i].adaptation, N, 0);
        set_all_values(S->nw.layer[i].adaptation_old, N, 0);
        set_all_values(S->nw.layer[i].facilitation, N, U0);
        set_all_values(S->nw.layer[i].facilitation_old, N, U0);
        set_all_values(S->nw.layer[i].depression, N, 1.0);
        set_all_values(S->nw.layer[i].depression_old, N, 1.0);
    }
    reset_all_bounds(S);
    set_all_values(S->nw.I_stim, N, 0.0);
    set_all_values(S->nw.I_stim_ahead, N, 0.0);
}

void initialize_rate_profile(struct State *S, double r_max, double m)
{
    int N = S->pconfig.N;
    int N_layers = S->pconfig.N_layers;
    double theta, tmp;
    for (int i = 0; i < N; i++) {
        theta = (-0.5 + (double) i / (double) N) * M_PI;
        /* Profile centered at Theta=0 */
        tmp = R + (r_max - R) * exp(m * (cos(2 * theta) - 1));
        for (int j = 0; j < N_layers; j++)
            S->nw.layer[j].rate[i] = tmp;
    }
}

void initialize_network_2D(struct State *S)
{
    /* Make the code a bit more readable */
    struct Parameters *pc = &S->pconfig;
    int N = pc->N;
    int N1 = pc->N1;
    int N2 = pc->N2;
    int N_layers = pc->N_layers;

    if (pc->N2 == 1 || pc->N1 == 1) {
        printf
            ("The network is actually 1D. Use 'simulate_trajectory' instead, please.\n");
        exit(1);
    } else {
        S->is_2D = true;
    }

    /* allocate space */
    /* ============== */
    S->nw.layer = (struct Layer *) malloc(N_layers * sizeof(struct Layer));
    for (int i = 0; i < N_layers; i++) {
        S->nw.layer[i].spars = malloc(N * sizeof(struct Sgmoid_pars));
        S->nw.layer[i].rate = malloc(N * sizeof(double));
        S->nw.layer[i].rate_old = malloc(N * sizeof(double));
        S->nw.layer[i].adaptation = malloc(N * sizeof(double));
        S->nw.layer[i].adaptation_old = malloc(N * sizeof(double));
        S->nw.layer[i].kernel = malloc(N * sizeof(double));
        /* S->nw.layer[i].kernel_mat = malloc(N * sizeof(double *)); */
        /* for (int j = 0; j < N; j++) */
            /* S->nw.layer[i].kernel_mat[j] = malloc(N * sizeof(double)); */
    }
    S->nw.I_firstlayer = tune_I_firstlayer(R, pc);
    S->nw.I_innerlayer = tune_I_innerlayer(R, pc);
    S->nw.ff_kernel = malloc(N * sizeof(double));
    S->nw.I_stim = malloc(N * sizeof(double));
    S->nw.I_stim_ahead = malloc(N * sizeof(double));

    /* and then, assign */
    struct Sgmoid_pars *sp = &pc->sgmoid;
    struct Connectivity_pars *cp = &pc->kernel;
    struct Connectivity_pars *cpff = &pc->ff_kernel;

    /* inhomogeneity parameters */
    /* double inhomog_lambda = 0.0; */
    /* double inhomog_kappa = 0.0; */
    for (int i = 0; i < N_layers; i++) {
        initialize_sigmoid_pars(S->nw.layer[i].spars, N, sp);
        initialize_kernel_vector_2D(S->nw.layer[i].kernel, N1, N2, cp);
        /* initialize_kernel_matrix(S->nw.layer[i].kernel_mat, S->nw.layer[i].kernel, N, inhomog_lambda, inhomog_kappa); */
        set_all_values(S->nw.layer[i].rate, N, R);
        set_all_values(S->nw.layer[i].rate_old, N, R);
        set_all_values(S->nw.layer[i].adaptation, N, 0);
        set_all_values(S->nw.layer[i].adaptation_old, N, 0);
    }
    set_all_values(S->nw.I_stim, N, 0.0);
    set_all_values(S->nw.I_stim_ahead, N, 0.0);
    S->pulse_list = NULL;
    S->pulse_list_2D = NULL;

    if (N_layers > 1)
        initialize_kernel_vector_2D(S->nw.ff_kernel, N1, N2, cpff);
    S->nw.I_max = max_amplitude_pulse_2D(S->pulse_list_2D);
    S->nw.max_evoked_rate = 0.0;
    S->nw.min_evoked_rate = pc->sgmoid.max_rate;

    if (pc->verbose) {
        show_fourier_coeffs_2D(pc);
    }
}

void initialize_circular_buffer(struct State *S)
{
    /* Make the code a bit more readable */
    struct Parameters *pc = &S->pconfig;
    int N = pc->N;
    int N_layers = pc->N_layers;

    /* Circular buffers storing the network activity and the inputs. */
    S->rate_images = malloc(N_layers * sizeof(struct circ_buffer));
    S->adaptation_images = malloc(N_layers * sizeof(struct circ_buffer));
    S->depression_images = malloc(N_layers * sizeof(struct circ_buffer));
    S->facilitation_images = malloc(N_layers * sizeof(struct circ_buffer));
    for (int i = 0; i < N_layers; i++) {
        S->rate_images[i] = circ_buffer_create(TIME_WINDOW_WIDTH, N);
        S->adaptation_images[i] = circ_buffer_create(TIME_WINDOW_WIDTH, N);
        S->depression_images[i] = circ_buffer_create(TIME_WINDOW_WIDTH, N);
        S->facilitation_images[i] = circ_buffer_create(TIME_WINDOW_WIDTH, N);
    }
    S->input_image = circ_buffer_create(TIME_WINDOW_WIDTH, N);
}

void simulate_single_trajectory(struct State *S)
{
    /* ----------------------------------------------- *
     *  Actual simulation
     * ----------------------------------------------- */
    double duration = S->pconfig.T;
    double dt = S->pconfig.dt;
    int N_layers = S->pconfig.N_layers;
    double inter_sample_interval =
        (double) duration / (double) TIME_WINDOW_WIDTH;
    if (inter_sample_interval <= dt) {
        printf("the inverse of the sampling frequency is shorter than dt.\n"
               "Change parameters conveniently\n");
        exit(1);
    }

    double time = 0.0;
    double next_sample_time = 0.0;
    bool after_perturbation = false;

    /* we need this to know the inputs for the following timestep */
    S->pulse_list = inputs_update(time, S);
    advance_inputs_one_step(S);
    
    while (time <= duration) {
        /* Look inputs ahead */
        S->pulse_list = inputs_update(time + dt, S);
        if (time < 10) { 
                white_noise_input_update(S); /* this is a temporary line */
        } else if (!after_perturbation) {
                set_all_values(S->nw.I_stim_ahead, S->pconfig.N, 0.0);
                after_perturbation = true;
        }
        network_update_RK4(S);
        if (time >= next_sample_time) {
            for (int i = 0; i < N_layers; i++) {
                circ_buffer_add_vector(&S->rate_images[i], S->nw.layer[i].rate);
                circ_buffer_add_vector(&S->adaptation_images[i], S->nw.layer[i].adaptation);
                circ_buffer_add_vector(&S->depression_images[i], S->nw.layer[i].depression);
                circ_buffer_add_vector(&S->facilitation_images[i], S->nw.layer[i].facilitation);
            }
            circ_buffer_add_vector(&S->input_image, S->nw.I_stim);
            next_sample_time += inter_sample_interval;
        }
        time += dt;
        advance_inputs_one_step(S);
    }
}

void advance_inputs_one_step(struct State *S)
{
    for (int i = 0; i < S->pconfig.N; i++) {
        S->nw.I_stim[i] = S->nw.I_stim_ahead[i];
    }
}

void simulate_single_trajectory_self_stop(struct State *S, double *time_formation, double *angle_location)
{
    /* ----------------------------------------------- *
     *  Actual simulation
     * ----------------------------------------------- */
    double duration = S->pconfig.T;
    double dt = S->pconfig.dt;
    double time = 0.0;
    double threshold = 0.9;

    advance_inputs_one_step(S);

    while (time <= duration) {
        S->pulse_list = inputs_update(time, S);
        network_update_Euler(S);
        if (S->nw.max_evoked_rate >= threshold) {
            break;
        }
        time += dt;
        advance_inputs_one_step(S);
    }
    *time_formation = time;
    *angle_location = get_angle_firing_highest(S, threshold);
}

void simulate_single_trajectory_2D(struct State *S)
{
    /* This function is like the 1D counterpart, except
     * that we use network_update_2D and that we use much
     * fewer snapshots (time samples).
     * ----------------------------------------------- */
    double duration = S->pconfig.T;
    double dt = S->pconfig.dt;
    int N_layers = S->pconfig.N_layers;
    /* We don't need many shots to figure out what's going on. */
    double inter_sample_interval = (double) duration / N_SNAPSHOTS;


    double time = 0.0;
    double next_sample_time = 0.0;

    while (time < duration) {
        S->pulse_list_2D = inputs_update_2D(time, S);
        network_update_2D(S);
        if (time >= next_sample_time) {
            for (int i = 0; i < N_layers; i++) {
                circ_buffer_add_vector(&S->rate_images[i], S->nw.layer[i].rate);
                circ_buffer_add_vector(&S->adaptation_images[i], S->nw.layer[i].adaptation);
                circ_buffer_add_vector(&S->depression_images[i], S->nw.layer[i].depression);
                circ_buffer_add_vector(&S->facilitation_images[i], S->nw.layer[i].facilitation);
            }
            circ_buffer_add_vector(&S->input_image, S->nw.I_stim);
            next_sample_time += inter_sample_interval;
        }
        time += dt;
    }
}

double get_angle_firing_highest(struct State *S, double threshold)
{
    /* Layer 0 only */
    double tmp_i = 0.0;
    int tmp_n = 0.0;
    for (int i = 0; i < S->pconfig.N; i++) {
        if (S->nw.layer[0].rate[i] > threshold) {
            tmp_i += i;
            tmp_n++;
        }
    }
    tmp_i = tmp_i / (double) tmp_n;
    return (-0.5 + (double) tmp_i / (double) S->pconfig.N) * M_PI;
}

void save_dynamic_variables(struct State *S)
{
    save_dynamic_variables_in_text(S);
    save_dynamic_variables_in_bitmaps(S);
}

void save_dynamic_variables_in_text(struct State *S)
{
    /* -----------------------------------------------
     *  Save firing rate activity across time. 
     * ----------------------------------------------- */
    /* Save data in text format */

    char filename[50];
    char filename_popv[50];
    struct Parameters *pc = &S->pconfig;

    for (int i = 0; i < pc->N_layers; i++) {
        /* Save rates */
        set_filename(filename, pc->output_rate_file, i, ".dat");
        save(&S->rate_images[i], filename);

        /* Save adaptation variable */
        set_filename(filename, pc->output_adaptation_file, i, ".dat");
        save(&S->adaptation_images[i], filename);

        /* Save depression variable */
        set_filename(filename, pc->output_depression_file, i, ".dat");
        save(&S->depression_images[i], filename);

        /* Save facilitation variable */
        set_filename(filename, pc->output_facilitation_file, i, ".dat");
        save(&S->facilitation_images[i], filename);

        /* Save population vectors */
        set_filename(filename_popv, "pop_vector", i, ".dat");
        save_population_vector(&S->rate_images[i], filename_popv);
    }
}

void save_dynamic_variables_in_bitmaps(struct State *S)
{
    char filename[50];
    struct Parameters *pc = &S->pconfig;
    struct Network *nw = &S->nw;

    /* Save data in image (colored) format */
    /* For B\W figures use save_pgm instead of save_circ_buffer_ppm
     * and modify the filename extension accordingly. You
     * may also want to use '-negate' as an option for 'convert' */
    /* determine_bounds(S); */
    for (int i = 0; i < pc->N_layers; i++) {
        set_filename(filename, pc->output_rate_file, i, ".ppm");
        save_ppm(&S->rate_images[i], filename, nw->min_evoked_rate, nw->max_evoked_rate);
        printf("Min rate = %6.3f,  Max rate = %6.3f\n", nw->min_evoked_rate, nw->max_evoked_rate);

        set_filename(filename, pc->output_adaptation_file, i, ".ppm");
        save_ppm(&S->adaptation_images[i], filename, nw->min_adaptation, nw->max_adaptation);

        set_filename(filename, pc->output_depression_file, i, ".ppm");
        save_ppm(&S->depression_images[i], filename, nw->min_depression, nw->max_depression);

        set_filename(filename, pc->output_facilitation_file, i, ".ppm");
        save_ppm(&S->facilitation_images[i], filename, nw->min_facilitation, nw->max_facilitation);
    }
}

void reset_all_bounds(struct State *S)
{
    /* Set max and min values to sensible initial values,
     * meaning that we make sure they will be updated during the simulation. */
    struct Parameters *pc = &S->pconfig;

    S->nw.I_max = max_amplitude_pulse(S->pulse_list);
    S->nw.max_evoked_rate = 0.0;
    S->nw.min_evoked_rate = pc->sgmoid.max_rate;
    if (pc->negfbck.J_A != 0.0) {
        S->nw.max_adaptation = 0.0;
        S->nw.min_adaptation = pc->sgmoid.max_rate;
    } else {
        S->nw.max_adaptation = 0.0;
        S->nw.min_adaptation = 0.0;
    }
    if (pc->negfbck.facilitation_f) {
        S->nw.max_facilitation = pc->negfbck.U;
        S->nw.min_facilitation = 1.0;
    } else {
        S->nw.max_facilitation = pc->negfbck.U;
        S->nw.min_facilitation = pc->negfbck.U;
    }
    if (pc->negfbck.depression_f) {
        S->nw.max_depression = 0.0;
        S->nw.min_depression = 1.0;
    } else {
        S->nw.max_depression = 1.0;
        S->nw.min_depression = 1.0;
    }
}

void determine_bounds(struct State *S)
{
    struct Parameters *pc = &S->pconfig;
    struct Network *nw = &S->nw;
    if (!pc->normalize) {
        nw->min_evoked_rate = 0.0;
        nw->max_evoked_rate = pc->sgmoid.max_rate;
        nw->min_adaptation = 0.0;
        nw->max_adaptation = pc->sgmoid.max_rate;
        nw->min_facilitation = pc->negfbck.U;
        nw->max_facilitation = 1.0;
        nw->min_depression = 0;
        nw->max_depression = 1.0;
    } else {
    /* Even if we ask for normalization, it makes no sense
     * to normalize constant variables */
        if (pc->negfbck.J_A == 0) {
            nw->min_adaptation = 0.0;
            nw->max_adaptation = pc->sgmoid.max_rate;
        }
        if (pc->negfbck.facilitation_f) {
            nw->min_facilitation = pc->negfbck.U;
            nw->max_facilitation = 1.0;
        }
        if (pc->negfbck.depression_f) {
            nw->min_depression = 0;
            nw->max_depression = 1.0;
        }
    }
}

void save_ppm(const struct circ_buffer *buff, const char *fname, double
        lobound, double hibound) 
{
    int status = 0;
    char filename[50];
    char command_str[100];

    save_circ_buffer_ppm(buff, fname, lobound, hibound);
    
    replace_substring(filename, fname, ".ppm", ".jpg");

    /* WARNING: This is extremely system dependent. I'm assuming
     * "imagemagick" is installed in your system.*/
    /* sprintf(command_str, "convert %s -transpose %s", fname, filename); */
    sprintf(command_str, "convert %s -rotate -90 %s", fname, filename); 
    status = system(command_str);
    if (status != 0)
        exit(status);
}

void save_snapshots_firing_activity_2D(struct State *S)
{
    /* -----------------------------------------------
     *  Save N_SNAPSHOTS frames of activity during
     *  a trial. Useful for 2D networks.
     * ----------------------------------------------- */

    /* Save data in text format */
    char filename[50];
    strcpy(filename, S->pconfig.output_rate_file);

    struct Parameters *pc = &S->pconfig;
    struct Network *nw = &S->nw;

    /* Save data in image (colored) format */
    char tmp[50];
    double lobound, hibound;
    for (int i = 0; i < pc->N_layers; i++) {
        sprintf(tmp, "%s_%02d", filename, i);
        if (pc->normalize) {
            lobound = nw->min_evoked_rate;
            hibound = nw->max_evoked_rate;
        } else {
            lobound = 0;
            hibound = pc->sgmoid.max_rate;
        }
        save_2D(&S->rate_images[i], tmp, pc->N1, pc->N2);
        save_2D_ppm(&S->rate_images[i], tmp, lobound, hibound, pc->N1, pc->N2);

    }
}


void save_inputs(struct State *S)
{
    struct Parameters *pc = &S->pconfig;
    /* TODO avoid using fix-sized arrays */
    char command_str[100];
    char filename[50];
    char tmp1[50];
    char tmp2[50];

    /* Save input vs time */
    strcpy(filename, pc->extinputs_file);
    sprintf(tmp1, "%s.dat", filename);
    save(&S->input_image, tmp1);

    /* Generate ppm format */
    sprintf(tmp1, "%s.ppm", filename);
    save_circ_buffer_ppm(&S->input_image, tmp1, 0, S->nw.I_max); // TODO use config info

    /* Generate image */
    int status = 0;
    sprintf(tmp2, "%s.jpg", filename);
    sprintf(command_str, "convert %s -rotate -90 %s", tmp1, tmp2);
    status = system(command_str);
    if (status != 0)
        exit(status);
}

void free_network_variables(struct State *S)
{
    struct Parameters *pc = &S->pconfig;
    struct Network *nw = &S->nw;
    int N_layers = pc->N_layers;
    int N = pc->N;

    for (int i = 0; i < N_layers; i++) {
        free(nw->layer[i].spars);
        free(nw->layer[i].rate);
        free(nw->layer[i].rate_old);
        free(nw->layer[i].adaptation);
        free(nw->layer[i].adaptation_old);
        free(nw->layer[i].kernel);
        for (int j = 0; j < N; j++)
            free(nw->layer[i].kernel_mat[j]);
        free(nw->layer[i].kernel_mat);
        circ_buffer_free(&S->rate_images[i]);
        circ_buffer_free(&S->adaptation_images[i]);
        circ_buffer_free(&S->depression_images[i]);
        circ_buffer_free(&S->facilitation_images[i]);
    }
    free(nw->layer);
    free(nw->I_stim);
    free(nw->I_stim_ahead);
    free(nw->ff_kernel);
    free(Rng);
    circ_buffer_free(&S->input_image);
}

void free_input_variables(struct State *S)
{
    free_extinputs(S->pulse_list);
    free_extinputs_2D(S->pulse_list_2D);
}

/* lower-level functions */

double sigmoid(double x, const struct Sgmoid_pars *spars)
{
    return spars->max_rate / (1.0 + exp(-spars->beta * (x - spars->thr)));
}

double sigmoid_inv(double r, const struct Sgmoid_pars *spars)
{
    /* This function returns the "input" necessary to have a SU at value r */
    return spars->thr - (1.0 / spars->beta) * log(spars->max_rate / r - 1.0);
}

double tune_I_firstlayer(double r, const struct Parameters *pconfig)
{
    double J0 = pconfig->kernel.J_E - pconfig->kernel.J_I;
    return sigmoid_inv(r, &pconfig->sgmoid) - J0 * r;
}

double tune_I_innerlayer(double r, const struct Parameters *pconfig)
{
    double J0_FF = pconfig->ff_kernel.J_E - pconfig->ff_kernel.J_I;
    return tune_I_firstlayer(r, pconfig) - J0_FF * r;
}

void set_all_values(double *v, int N, double x)
{
    for (int i = 0; i < N; i++) {
        v[i] = x;
    }
}

void initialize_kernel_vector(double *krnl, int N,
                              const struct Connectivity_pars *cp)
{
    /* Feedforward conn kernel is identical */
    double theta;
    /* Just for readability */
    double Je = cp->J_E;
    double Ji = cp->J_I;
    double me = cp->m_E;
    double mi = cp->m_I;

    for (int i = 0; i <= N / 2; i++) {
        theta = ((double) i / (double) N) * M_PI;
        /* Use symmetry */
        krnl[i] = Je * exp(me * cos(2 * theta)) / gsl_sf_bessel_I0(me)
            - Ji * exp(mi * cos(2 * theta)) / gsl_sf_bessel_I0(mi);
        if (i > 0 && i < N / 2)
            krnl[N - i] = krnl[i];
    }
}

void initialize_kernel_vector_2D(double *krnl, int N1, int N2,
                                 const struct Connectivity_pars *cp)
{
    double theta1[N1];
    double theta2[N2];

    /* Just for readability */
    double Je = cp->J_E;
    double Ji = cp->J_I;
    double me = cp->m_E;
    double mi = cp->m_I;

    for (int i = 0; i < N1; i++)
        theta1[i] = ((double) i / (double) N1) * M_PI;
    for (int i = 0; i < N2; i++)
        theta2[i] = ((double) i / (double) N2) * M_PI;


    /* Compute numerically the normalization constants
     * of the exc. and inh. footprints */
    /* double N_E = integral_on_torus(me) / (M_PI * M_PI); */
    /* double N_I = integral_on_torus(mi) / (M_PI * M_PI); */
    double N_E = sum_on_torus(me, N1, N2) / (N1 * N2);
    double N_I = sum_on_torus(mi, N1, N2) / (N1 * N2);

    double dist;
    for (int i = 0; i < N1; i++) {
        for (int j = 0; j < N2; j++) {
            dist = distance_on_torus(0, 0, theta1[i], theta2[j]);
            krnl[N2 * i + j] = Je * exp(me * cos(2 * dist)) / N_E
                - Ji * exp(mi * cos(2 * dist)) / N_I;
        }
    }
    /* Use symmetry */
    /* for (int i = 0; i <= N1 / 2; i++) { */
        /* for (int j = 0; j <= N2 / 2; j++) { */
            /* dist = sqrt(pow(theta1[i], 2) + pow(theta2[j], 2)); */
            /* [> upper left quarter of the whole matrix <] */
            /* krnl[N2 * i + j] = Je * exp(me * cos(2 * dist)) / N_E */
                /* - Ji * exp(mi * cos(2 * dist)) / N_I; */
            /* matrix divided in quarters, each of which is a mirrored
             * copy of one another. */
            /* [> lower left <] */
            /* if (i > 0 && i < N1 / 2) */
                /* krnl[(N1 - i) * N2 + j] = krnl[N2 * i + j]; */
            /* [> upper right <] */
            /* if (j > 0 && j < N2 / 2) */
                /* krnl[N2 * i + (N2 - j)] = krnl[N2 * i + j]; */
            /* [> lower right <] */
            /* if (i > 0 && j > 0 && i < N1 / 2 && j < N2 / 2) */
                /* krnl[(N1 - i) * N2 + (N2 - j)] = krnl[N2 * i + j]; */
        /* } */
    /* } */
}

void initialize_modulation (double *v, int N,
                              double lambda, double concentration)
{
    double theta;
    double favored_angle = 0.0; /* maximum local excitability at favored_angle */
    for (int i = 0; i < N; i++) {
        theta = (-0.5 + (double) i / (double) N) * M_PI;
        v[i] = 1.0 + lambda * exp(concentration * (cos(2 * (theta - favored_angle)) - 1.0));
    }
}

void initialize_kernel_matrix(double **krnl_mat, const double *krnl, int N,
                              double lambda, double concentration)
{
    /* "lambda" and "concentration" are, respectively, the amplification factor and the
     * concentration paramater of the modulation */

    double theta;
    double modulation;
    double favored_angle = 0.0; /* where the maximum local excitability is located */

    /* TODO get around the construction of this matrix
     * when the network is homogeneous */
    if (concentration == 0 && lambda == 0) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++)
                krnl_mat[i][j] = krnl[j];
        }
    } else {
        for (int i = 0; i < N; i++) {
            theta = (-0.5 + (double) i / (double) N) * M_PI;
            modulation =
                1.0 +
                lambda * exp(concentration *
                             (cos(2 * (theta - favored_angle)) - 1.0));
            for (int j = 0; j < N; j++)
                krnl_mat[i][j] = modulation * krnl[j];
        }
    }
}

void initialize_kernel_matrix_with_noise(double **krnl_mat, const double *krnl, int N,
                              double epsilon)
{
    /* "epsilon" is the std of the gaussian r.v. we add to the connectivity kernel */

    double sqrtN = sqrt(N);
    /* TODO get around the construction of this matrix
     * when the network is homogeneous */
    if (epsilon == 0) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++)
                krnl_mat[i][j] = krnl[j];
        }
    } else {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++)
                krnl_mat[i][j] = krnl[j] + sqrtN * gsl_ran_gaussian(Rng, epsilon);
        }
    }
}

void initialize_sigmoid_pars(struct Sgmoid_pars *spars, int N,
                             const struct Sgmoid_pars *sp)
{
    /* Inhomogeneities */
    /* double theta; */
    /* double extra = 3.0;  the maximum amount of gain/thr we add  */
    /* double modulation = 2; the concentration parameter of the modulation */

    for (int i = 0; i < N; i++) {
        /* theta = (-0.5 + (double) i / (double) N) * M_PI; */
        /* we modulate the beta parameter of the sigmoid */
        spars[i].beta = sp->beta;       /*  + extra * exp(modulation * ( cos(2 * theta) - 1 )) */
        spars[i].max_rate = sp->max_rate;       /* (1 + extra * exp(modulation * ( cos(2 * theta) - 1 )));  */
        spars[i].thr = sp->thr;
        /* spars[i].thr = sp->thr +  gsl_ran_gaussian(Rng, 1.0); */
    }
}

void show_fourier_coeffs(const struct Parameters *pconfig)
{
    printf("\nLinear stability analysis...\n");
    printf("  steady uniform rate, R=%4.2f\n", R);
    printf("  critical value Fourier mode, J_cr=%5.2f\n",
           1 / (pconfig->sgmoid.beta * R * (1 - R / pconfig->sgmoid.max_rate)));
    printf("  Fourier decomposition:\n");
    for (int i = 0; i < 10; i++) {
        printf("  \t%d\t%10.6f\n", i, pconfig->kernel.J_E *
               gsl_sf_bessel_In(i, pconfig->kernel.m_E) /
               gsl_sf_bessel_I0(pconfig->kernel.m_E)
               - pconfig->kernel.J_I *
               gsl_sf_bessel_In(i, pconfig->kernel.m_I) /
               gsl_sf_bessel_I0(pconfig->kernel.m_I));
    }
}

void show_fourier_coeffs_2D(const struct Parameters *pconfig)
{
    printf("\nLinear stability analysis for the 2D network...\n");
    printf("  steady uniform rate, R=%4.2f\n", R);
    printf("  critical value Fourier mode, J_cr=%5.2f\n",
           1 / (pconfig->sgmoid.beta * R * (1 - R / pconfig->sgmoid.max_rate)));
    printf("  Fourier decomposition (rows: k_1, columns: k_2)\n");
    printf("\n  \t    %d % 10d % 10d\n", 0, 1, 2);

    double m_E = pconfig->kernel.m_E;
    double m_I = pconfig->kernel.m_I;
    double J_E = pconfig->kernel.J_E;
    double J_I = pconfig->kernel.J_I;

    double N_E = integral_on_torus(m_E);
    double N_I = integral_on_torus(m_I);
    double J;
    for (int i = 0; i < 10; i++) {
        printf("  \t%d", i);
        for (int j = 0; j < 3; j++) {
            J = (J_E * Fourier_integral_2D(m_E, i, j) / N_E
                - J_I * Fourier_integral_2D(m_I, i, j) / N_I);
            printf(" % 10.6f", J);
        }
        printf("\n");
    }
}

void initialize_rng(void)
{
    /* Set up random generator */
    gsl_rng_env_setup();
    Rng_T = gsl_rng_default;
    Rng = gsl_rng_alloc(Rng_T);
}

void white_noise_input_update(struct State *S)
{
        int N = S->pconfig.N;
        for (int i = 0; i < N; i++) {
                S->nw.I_stim_ahead[i] = 1e-5 * gsl_ran_gaussian(Rng, 1.0);
                if (S->nw.I_stim_ahead[i] > S->nw.I_max)
                        S->nw.I_max = S->nw.I_stim_ahead[i];
        }
}

struct pulse *inputs_update(double t, struct State *S)
{
    struct pulse *p, *prev;
    int N = S->pconfig.N;

    /* Check if we have to switch any pulse on */
    prev = NULL;
    for (p = S->pulse_list; p != NULL; p = p->next) {
        /* This pulse is new */
        if (p->switched_on == 0 && t > p->time_onset) {
            p->switched_on = 1;
            /* Add the new pulse to the input vector */
            for (int i = 0; i < N; i++) {
                S->nw.I_stim_ahead[i] += p->precomputed_value[i];
                if (S->nw.I_stim_ahead[i] > S->nw.I_max)
                    S->nw.I_max = S->nw.I_stim_ahead[i];
            }
        }
        /* This pulse is stale */
        if (p->switched_on == 1 && t > p->time_onset + p->duration) {
            /* Remove pulse */
            for (int i = 0; i < N; i++)
                S->nw.I_stim_ahead[i] -= p->precomputed_value[i];
            /* Be nice and clean */
            /* Remove pulse from the list. See Kernighan & Pike p49 */
            if (prev == NULL)
                S->pulse_list = p->next;
            else
                prev->next = p->next;   /* bridge out the element */
            free(p->precomputed_value);
            free(p);
            continue;
        }
        prev = p;
    }
    return S->pulse_list;
}

struct pulse_2D *inputs_update_2D(double t, struct State *S)
{
    /* This is a verbatim copy of inputs_update, changing
     * pulse with pulse_2D */
    struct pulse_2D *p, *prev;
    int N = S->pconfig.N;

    /* Check if we have to switch any pulse on */
    prev = NULL;
    for (p = S->pulse_list_2D; p != NULL; p = p->next) {
        /* This pulse is new */
        if (p->switched_on == 0 && t > p->time_onset) {
            p->switched_on = 1;
            /* Add the new pulse to the input vector */
            for (int i = 0; i < N; i++) {
                S->nw.I_stim_ahead[i] += p->precomputed_value[i];
                if (S->nw.I_stim_ahead[i] > S->nw.I_max)
                    S->nw.I_max = S->nw.I_stim_ahead[i];
            }
        }
        /* This pulse is stale */
        if (p->switched_on == 1 && t > p->time_onset + p->duration) {
            /* Remove pulse */
            for (int i = 0; i < N; i++)
                S->nw.I_stim_ahead[i] -= p->precomputed_value[i];
            /* Be nice and clean */
            /* Remove pulse from the list. See Kernighan & Pike p49 */
            if (prev == NULL)
                S->pulse_list_2D = p->next;
            else
                prev->next = p->next;   /* bridge out the element */
            free(p->precomputed_value);
            free(p);
        }
        prev = p;
    }
    return S->pulse_list_2D;
}

void network_update_Euler(struct State *S)
{
    double rec_input, ff_input;
    int N = S->pconfig.N;
    int N_layers = S->pconfig.N_layers;
    double dt = S->pconfig.dt;
    struct Network *nw = &S->nw;
    struct NegFeedback_pars *negfbck = &S->pconfig.negfbck;
    struct Layer *l;

    /* rold_{t} = r_{t-1} */
    for (int i = 0; i < N_layers; i++) {
        for (int j = 0; j < N; j++) {
            nw->layer[i].rate_old[j] = nw->layer[i].rate[j];
            nw->layer[i].adaptation_old[j] = nw->layer[i].adaptation[j];
            nw->layer[i].depression_old[j] = nw->layer[i].depression[j];
            nw->layer[i].facilitation_old[j] = nw->layer[i].facilitation[j];
        }
    }
    for (int i = 0; i < N_layers; i++) {
        l = &nw->layer[i];
        for (int j = 0; j < N; j++) {
            /* Recurrent input: spatial convolution */
            rec_input = 0.0;
            for (int k = 0; k < N; k++) {
                /* Homogeneous network */
                /* rec_input += l->rate_old[k] * l->facilitation[k] * l->depression[k] * l->kernel[abs(j - k)]; */
                /* Heterogeneous network */
                rec_input += l->rate_old[k] * l->facilitation[k] * l->depression[k] * l->kernel_mat[j][abs(j - k)];
                /* rec_input += l->rate_old[k] * l->kernel_mat[j][abs(j - k)]; */
            }
            rec_input /= N;
            /* External input enters in the input (1st) layer */
            if (i == 0) {
                l->rate[j] +=
                    dt * (-l->rate_old[j] +
                          sigmoid(rec_input +nw->I_firstlayer + nw->I_stim[j], l->spars));
            } else {
                /* Feedforward input propagated to remaining layers */
                ff_input = 0.0;
                for (int k = 0; k < N; k++) {
                    ff_input +=
                        nw->layer[i -
                                  1].rate_old[k] * nw->ff_kernel[abs(j - k)];
                }
                ff_input = ff_input / N;
                l->rate[j] +=
                    dt * (-l->rate_old[j] +
                          sigmoid(rec_input + nw->I_innerlayer + ff_input,
                                  l->spars));
            }

            /* store lower and upper bounds of activity */
            if (l->rate[j] > nw->max_evoked_rate) {
                nw->max_evoked_rate = l->rate[j];
            } else if (l->rate[j] < nw->min_evoked_rate) {
                nw->min_evoked_rate = l->rate[j];
            }

            /* Update adaptation current */
            if (negfbck->J_A != 0) {
                l->adaptation[j] += dt * ( - l->adaptation_old[j]
                        + negfbck->J_A * ( l->rate_old[j] ) ) / negfbck->tau_A;

                if (l->adaptation[j] > nw->max_adaptation) {
                    nw->max_adaptation = l->adaptation[j];
                } else if (l->adaptation[j] < nw->min_adaptation) {
                    nw->min_adaptation = l->adaptation[j];
                }
            } else {
                l->adaptation[j] = l->adaptation_old[j];
            }

            /* Update depression (available resources) variable */
            if (negfbck->depression_f) {
                l->depression[j] += dt * ( (1.0 - l->depression_old[j]) / negfbck->tau_D
                        - l->depression_old[j] * l->facilitation_old[j] * l->rate_old[j] ) ;
                if (l->depression[j] > nw->max_depression) {
                    nw->max_depression = l->depression[j];
                } else if (l->depression[j] < nw->min_depression) {
                    nw->min_depression = l->depression[j];
                }
            } else {
                l->depression[j] = l->depression_old[j];
            }

            /* Update facilitation (fraction of resources used per spike) variable */
            if (negfbck->facilitation_f) {
                l->facilitation[j] += dt * ( - (l->facilitation_old[j] - negfbck->U) / negfbck->tau_F
                                            + negfbck->U * l->rate_old[j] * (1 - l->facilitation_old[j] ) ) ;
                if (l->facilitation[j] > nw->max_facilitation) {
                    nw->max_facilitation = l->facilitation[j];
                } else if (l->facilitation[j] < nw->min_facilitation) {
                    nw->min_facilitation = l->facilitation[j];
                }
            } else {
                l->facilitation[j] = l->facilitation_old[j];
            }

        }
    }
}

void network_update_RK4(struct State *S)
{
    double rec_input, ff_input;
    int N = S->pconfig.N;
    int N_layers = S->pconfig.N_layers;
    double dt = S->pconfig.dt;
    struct Network *nw = &S->nw;
    struct NegFeedback_pars *negfbck = &S->pconfig.negfbck;
    struct Layer *l;

    /* Aux vars for RK4. Four N-dimensional variables: rate, adaptation, facilitation, and depression */
    double k1[4][N]; /*[N]*/
    double k2[4][N];
    double k3[4][N];
    double k4[4][N];
    double r;
    double a;
    double f;
    double d;

    /* rold_{t} = r_{t-1} */
    for (int i = 0; i < N_layers; i++) {
        for (int j = 0; j < N; j++) {
            nw->layer[i].rate_old[j] = nw->layer[i].rate[j];
            nw->layer[i].adaptation_old[j] = nw->layer[i].adaptation[j];
            nw->layer[i].depression_old[j] = nw->layer[i].depression[j];
            nw->layer[i].facilitation_old[j] = nw->layer[i].facilitation[j];
        }
    }
    for (int i = 0; i < N_layers; i++) {
        l = &nw->layer[i];
        for (int j = 0; j < N; j++) {
            /* ===========================================================================*/
            /* 1st step RK4                                                               */
            /* ===========================================================================*/
            /* Recurrent input: spatial convolution */
            rec_input = 0.0;
            for (int k = 0; k < N; k++) {
                /* -- Homogeneous network -- */
                /* rec_input += l->rate_old[k] * l->facilitation_old[k] * l->depression_old[k] * l->kernel[abs(j - k)]; */
                /* -- Heterogeneous network -- */
                rec_input += l->rate_old[k] * l->facilitation_old[k] * l->depression_old[k] * l->kernel_mat[j][abs(j - k)];
            }
            rec_input /= N;
            if (i == 0) {
                k1[0][j] = dt * (-l->rate_old[j] + sigmoid(rec_input - l->adaptation_old[j]
                                        +nw->I_firstlayer + nw->I_stim[j], l->spars));
            } else {
                /* TODO Only layer 0 is computed with RK4. The remaining layers use Euler. Fix that. */
                /* Feedforward input propagated to remaining layers */
                ff_input = 0.0;
                for (int k = 0; k < N; k++) {
                    ff_input +=
                        nw->layer[i -
                                  1].rate_old[k] * nw->ff_kernel[abs(j - k)];
                }
                ff_input = ff_input / N;
                k1[0][j] = dt * (-l->rate_old[j] +
                          sigmoid(rec_input + nw->I_innerlayer + ff_input,
                                  l->spars));
            }

            /* Update adaptation current */
            if (negfbck->J_A != 0)
                k1[1][j] = dt * ( -l->adaptation_old[j] + negfbck->J_A * ( l->rate_old[j] ) ) / negfbck->tau_A;
            else
                k1[1][j] = 0;

            /* Update facilitation (fraction of resources used per spike) variable */
            if (negfbck->facilitation_f)
                k1[2][j] = dt * ( - (l->facilitation_old[j] - negfbck->U) / negfbck->tau_F
                                            + negfbck->U * l->rate_old[j] * (1 - l->facilitation_old[j] ) ) ;
            else
                k1[2][j] = 0;

            /* Update depression (available resources) variable */
            if (negfbck->depression_f)
                k1[3][j] = dt * ( (1.0 - l->depression_old[j]) / negfbck->tau_D
                        - l->depression_old[j] * l->facilitation_old[j] * l->rate_old[j] ) ;
            else
                k1[3][j] = 0;
        }

        /* ===========================================================================*/
        /* 2nd step RK4                                                               */
        /* ===========================================================================*/
        for (int j = 0; j < N; j++) {
            rec_input = 0.0;
            for (int k = 0; k < N; k++) {
                rec_input += (l->rate_old[k] + 0.5 * k1[0][k]) * (l->facilitation_old[k] + 0.5 * k1[2][k]) * (l->depression_old[k] + 0.5 * k1[3][k]) * l->kernel_mat[j][abs(j - k)];
            }
            rec_input /= N;

            r = l->rate_old[j] + 0.5 * k1[0][j];
            a = l->adaptation_old[j] + 0.5 * k1[1][j];
            f = l->facilitation_old[j] + 0.5 * k1[2][j];
            d = l->depression_old[j] + 0.5 * k1[3][j];
            /* External input enters in the input (1st) layer */
            if (i == 0) {
                k2[0][j] = dt * (- r + sigmoid(rec_input - a
                                  + nw->I_firstlayer + nw->I_stim[j], l->spars));
            } else {
                /* TODO Only layer 0 is computed with RK4. The remaining layers use Euler. Fix that. */
                /* Feedforward input propagated to remaining layers */
                ff_input = 0.0;
                for (int k = 0; k < N; k++) {
                    /* this is wrong; we are not updating properly */
                    ff_input += nw->layer[i - 1].rate_old[k] * nw->ff_kernel[abs(j - k)];
                }
                ff_input = ff_input / N;
                k2[0][j] = dt * (-r + sigmoid(rec_input + nw->I_innerlayer + ff_input,
                                  l->spars));
            }

            /* Update adaptation current */
            if (negfbck->J_A != 0)
                k2[1][j] = dt * ( -a + negfbck->J_A * r ) / negfbck->tau_A;
            else
                k2[1][j]  = 0 ;

            /* Update facilitation (fraction of resources used per spike) variable */
            if (negfbck->facilitation_f)
                k2[2][j] = dt * ( - (f - negfbck->U) / negfbck->tau_F + negfbck->U * r * (1 - f) ) ;
            else
                k2[2][j] = 0;

            /* Update depression (available resources) variable */
            if (negfbck->depression_f)
                k2[3][j] = dt * ( (1.0 - d) / negfbck->tau_D - d * f * r ) ;
            else
                k2[3][j] = 0;

        }

        /* ===========================================================================*/
        /* 3rd step RK4                                                               */
        /* ===========================================================================*/
        for (int j = 0; j < N; j++) {
            rec_input = 0.0;
            for (int k = 0; k < N; k++) {
                rec_input += (l->rate_old[k] + 0.5 * k2[0][k]) * (l->facilitation_old[k] + 0.5 * k2[2][k]) * (l->depression_old[k] + 0.5 * k2[3][k]) * l->kernel[abs(j - k)];
            }
            rec_input /= N;
            r = l->rate_old[j] + 0.5 * k2[0][j];
            a = l->adaptation_old[j] + 0.5 * k2[1][j];
            f = l->facilitation_old[j] + 0.5 * k2[2][j];
            d = l->depression_old[j] + 0.5 * k2[3][j];
            /* External input enters in the input (1st) layer */
            if (i == 0) {
                k3[0][j] = dt * (- r + sigmoid(rec_input - a
                                  + nw->I_firstlayer + nw->I_stim[j], l->spars));
            } else {
                /* TODO Only layer 0 is computed with RK4. The remaining layers use Euler. Fix that. */
                /* Feedforward input propagated to remaining layers */
                ff_input = 0.0;
                for (int k = 0; k < N; k++) {
                    /* this is wrong; we are not updating properly */
                    ff_input += nw->layer[i - 1].rate_old[k] * nw->ff_kernel[abs(j - k)];
                }
                ff_input = ff_input / N;
                k3[0][j] = dt * (-r + sigmoid(rec_input + nw->I_innerlayer + ff_input,
                                  l->spars));
            }

            /* Update adaptation current */
            if (negfbck->J_A != 0)
                k3[1][j] = dt * ( -a + negfbck->J_A * r ) / negfbck->tau_A;
            else
                k3[1][j]  = 0 ;

            /* Update facilitation (fraction of resources used per spike) variable */
            if (negfbck->facilitation_f)
                k3[2][j] = dt * ( - (f - negfbck->U) / negfbck->tau_F + negfbck->U * r * (1 - f)) ;
            else
                k3[2][j] = 0;

            /* Update depression (available resources) variable */
            if (negfbck->depression_f)
                k3[3][j] = dt * ( (1.0 - d) / negfbck->tau_D - d * f * r ) ;
            else 
                k3[3][j] = 0;
        }

        /* ===========================================================================*/
        /* 4th step RK4                                                               */
        /* ===========================================================================*/
        for (int j = 0; j < N; j++) {
            rec_input = 0.0;
            for (int k = 0; k < N; k++) {
                rec_input += (l->rate_old[k] + 0.5 * k3[0][k]) * (l->facilitation_old[k] + 0.5 * k3[2][k]) * (l->depression_old[k] + 0.5 * k3[3][k]) * l->kernel[abs(j - k)];
            }
            rec_input /= N;
            r = l->rate_old[j] + k3[0][j];
            a = l->adaptation_old[j] + k3[1][j];
            f = l->facilitation_old[j] + k3[2][j];
            d = l->depression_old[j] + k3[3][j];

            /* External input enters in the input (1st) layer */
            if (i == 0) {
                k4[0][j] = dt * (- r + sigmoid(rec_input - a
                                  + nw->I_firstlayer + nw->I_stim_ahead[j], l->spars));
            } else {
                /* TODO Only layer 0 is computed with RK4. The remaining layers use Euler. Fix that. */
                /* Feedforward input propagated to remaining layers */
                ff_input = 0.0;
                for (int k = 0; k < N; k++) {
                    /* this is wrong; we are not updating properly */
                    ff_input += nw->layer[i - 1].rate_old[k] * nw->ff_kernel[abs(j - k)];
                }
                ff_input = ff_input / N;
                k4[0][j] = dt * (-r + sigmoid(rec_input + nw->I_innerlayer + ff_input,
                                  l->spars));
            }

            /* Update adaptation current */
            if (negfbck->J_A != 0)
                k4[1][j] = dt * ( -a + negfbck->J_A * r ) / negfbck->tau_A;
            else
                k4[1][j]  = 0 ;

            /* Update facilitation (fraction of resources used per spike) variable */
            if (negfbck->facilitation_f)
                k4[2][j] = dt * ( - (f - negfbck->U) / negfbck->tau_F + negfbck->U * r * (1 - f)) ;
            else
                k4[2][j] = 0;

            /* Update depression (available resources) variable */
            if (negfbck->depression_f)
                k4[3][j] = dt * ( (1.0 - d) / negfbck->tau_D - d * f * r ) ;
            else 
                k4[3][j] = 0;

            l->rate[j] = l->rate_old[j] + (1.0/6.0) * (k1[0][j] + 2 * k2[0][j] + 2 * k3[0][j] + k4[0][j]);
            l->adaptation[j] = l->adaptation_old[j] + (1.0/6.0) * (k1[1][j] + 2 * k2[1][j] + 2 * k3[1][j] + k4[1][j]);
            l->facilitation[j] = l->facilitation_old[j] + (1.0/6.0) * (k1[2][j] + 2 * k2[2][j] + 2 * k3[2][j] + k4[2][j]);
            l->depression[j] = l->depression_old[j] + (1.0/6.0) * (k1[3][j] + 2 * k2[3][j] + 2 * k3[3][j] + k4[3][j]);

            /* store lower and upper bounds of activity */
            update_bounds_if_necessary(l->rate[j], &nw->min_evoked_rate, &nw->max_evoked_rate);
            if (negfbck->J_A != 0)
                update_bounds_if_necessary(l->adaptation[j], &nw->min_adaptation, &nw->max_adaptation);
            if (negfbck->facilitation_f)
                update_bounds_if_necessary(l->facilitation[j], &nw->min_facilitation, &nw->max_facilitation);
            if (negfbck->depression_f)
                update_bounds_if_necessary(l->depression[j], &nw->min_depression, &nw->max_depression);
        }
    }
}

void network_update_2D(struct State *S)
{
    double rec_input, ff_input;
    int N = S->pconfig.N;
    int N2 = S->pconfig.N2;
    int N_layers = S->pconfig.N_layers;
    double dt = S->pconfig.dt;
    struct Network *nw = &S->nw;


    /* rold_{t} = r_{t-1} */
    for (int i = 0; i < N_layers; i++) {
        for (int j = 0; j < N; j++)
            nw->layer[i].rate_old[j] = nw->layer[i].rate[j];
    }
    int post[2];    /* indices in the 2D representation */
    int pre[2]; 
    int id;
    for (int i = 0; i < N_layers; i++) {
        for (int j = 0; j < N; j++) {
            /* 2D indices of post-synaptic neuron */
            post[0] = j / N2;
            post[1] = j % N2;
            /* Recurrent input: spatial convolution */
            rec_input = 0.0;
            for (int k = 0; k < N; k++) {
                /* Homogeneous network */
                /* For a 2D network and for a 2d-array representation, the syn
                 * strength between two neurons (i,j) and (i', j') is
                 * kernel[|i - i'|][|j - j'|]. But we use a 1D representation: 
                 *              I  = N_2 * i  + j,
                 *              I' = N_2 * i' + j'                                     */
                pre[0] = k / N2;
                pre[1] = k % N2;
                id = abs(post[0] - pre[0]) * N2 + abs(post[1] - pre[1]);
                rec_input += nw->layer[i].rate_old[k] * nw->layer[i].kernel[id];
            }
            rec_input = rec_input / N;
            /* External input enters in the input (1st) layer */
            if (i == 0) {
                nw->layer[i].rate[j] +=
                    dt * (-nw->layer[i].rate_old[j] +
                          sigmoid(rec_input + nw->I_firstlayer + nw->I_stim[j],
                                  nw->layer[i].spars));
            } else {
                /* Feedforward input propagated to remaining layers */
                ff_input = 0.0;
                for (int k = 0; k < N; k++) {
                    pre[0] = k / N2;
                    pre[1] = k % N2;
                    id = abs(post[0] - pre[0]) * N2 + abs(post[1] - pre[1]);
                    ff_input +=
                        nw->layer[i - 1].rate_old[k] * nw->ff_kernel[id];
                }
                ff_input = ff_input / N;
                nw->layer[i].rate[j] +=
                    dt * (-nw->layer[i].rate_old[j] +
                          sigmoid(rec_input + nw->I_innerlayer + ff_input,
                                  nw->layer[i].spars));
            }
            /* store lower and upper bounds of activity */
            if (nw->layer[i].rate[j] > nw->max_evoked_rate) {
                nw->max_evoked_rate = nw->layer[i].rate[j];
            } else if (nw->layer[i].rate[j] < nw->min_evoked_rate) {
                nw->min_evoked_rate = nw->layer[i].rate[j];
            }
        }
    }
}

void update_bounds_if_necessary(double v, double *min, double *max) 
{
    if (v > *max) {
        *max = v;
    } else if (v < *min) {
        *min = v;
    }
}

void set_filename(char *filename, const char *basename, int i, const char *ext)
{
    // TODO Check sizes and all that
    strcpy(filename, basename);
    sprintf(filename, "%s_%02d%s", filename, i, ext);
}


void save(const struct circ_buffer *buff, const char *fname)
{
    FILE *fin = fopen(fname, "w");
    size_t n_rows;

    size_t offset = buff->start;
    if ( (buff->start - buff->next) % buff->size == 1 ) {
        /* the buffer has been filled */
        n_rows = buff->size - 1;
    } else {
        n_rows = buff->next;
    }
    for (size_t i = 0; i < n_rows - 1; i++) {
        for (size_t j = 0; j < buff->ncols; j++)
            fprintf(fin, "% 14.6e", buff->image[(i + offset) % buff->size][j]);
        fprintf(fin, "\n");
    }
    fclose(fin);
}

void save_population_vector(const struct circ_buffer *buff, const char *fname)
{
    FILE *fin = fopen(fname, "w");

    size_t offset = buff->start;
    size_t N = buff->ncols;

    double theta[N];            /* Thank god we are using C99! */
    double cos2theta[N];
    double sin2theta[N];
    /* precalculated values */
    for (size_t i = 0; i < N; i++) {
        theta[i] = (-0.5 + (double) i / (double) N) * M_PI;
        cos2theta[i] = cos(2 * theta[i]);
        sin2theta[i] = sin(2 * theta[i]);
    }

    double rate;
    double tmp_rate;
    double complex tmp_pop;
    for (size_t i = 0; i < buff->size - 1; i++) {
        tmp_rate = 0;
        tmp_pop = 0;
        for (size_t j = 0; j < N; j++) {
            rate = buff->image[(i + offset) % buff->size][j];
            tmp_pop += rate * (cos2theta[j] + I * sin2theta[j]);
            tmp_rate += rate;
        }
        tmp_pop = tmp_pop / tmp_rate;
        fprintf(fin, "% 10.4f % 10.4f", cabs(tmp_pop), carg(tmp_pop) / 2.0);
        fprintf(fin, "\n");
    }
    fclose(fin);
}

void save_pgm(const struct circ_buffer *buff, const char *fname, double lo,
              double hi, int N)
{
    /* Write it down */
    FILE *fin = fopen(fname, "w");
    /* Pgm format */
    fprintf(fin, "P2\n%d %d\n%d\n", N, TIME_WINDOW_WIDTH, MAXVAL);

    size_t offset = buff->start;
    for (size_t i = 0; i < buff->size - 1; i++) {
        for (size_t j = 0; j < buff->ncols; j++)
            fprintf(fin, "% 4d", (int) rint((double) MAXVAL *
                                            (buff->
                                             image[(i +
                                                    offset) % buff->size][j] -
                                             lo)
                                            / (hi - lo)));
        fprintf(fin, "\n");
    }
    fclose(fin);
}

void save_circ_buffer_ppm(const struct circ_buffer *buff, const char *fname,
                          double lo, double hi)
    /* Save data in ppm format (color) using the bone color scheme */
{
    /* Create the lookup table */
    /* unsigned int rgb_jet[MAXVAL + 1][3]; */
    /* fill_cm_lookuptable_jet(rgb_grey, MAXVAL); */
    unsigned int rgb_bone[MAXVAL + 1][3];
    fill_cm_lookuptable_bone(rgb_bone, MAXVAL);

    size_t N = buff->ncols;
    size_t n_rows;

    if ( (buff->start - buff->next) % buff->size == 1 ) {
        /* the buffer has been filled */
        n_rows = buff->size;
    } else {
        n_rows = buff->next;
    }

    /* Write data down */
    FILE *fin = fopen(fname, "w");
    /* Pgm format */
    fprintf(fin, "P3\n%d %d\n%d\n", (int)N, (int)n_rows, MAXVAL);

    size_t offset = buff->start;
    int pos;
    double normalization = hi - lo;
    if (hi == lo)
        normalization = 1.0;    /* Avoid division by zero */
    for (size_t i = 0; i < n_rows; i++) {
        for (size_t j = 0; j < buff->ncols; j++) {
            pos =
                (int) rint((double) MAXVAL *
                           (buff->image[(i + offset) % buff->size][j] -
                            lo) / normalization);
            fprintf(fin, "% 3d % 3d % 3d   ", rgb_bone[pos][0], rgb_bone[pos][1],
                    rgb_bone[pos][2]);
        }
        fprintf(fin, "\n");
    }
    fclose(fin);
}

void save_2D(const struct circ_buffer *buff, const char *fname, int N1, int N2)
/* Save data file with a snapshot of firing rates in the 2D network */
{
    char filename[50];
    strcpy(filename, fname);

    size_t offset;
    char tmp_base[50];
    char tmp_suff[50];

    for (int i = 0; i < N_SNAPSHOTS; i++) {
        sprintf(tmp_base, "%s_frame_%02d", filename, i);
        sprintf(tmp_suff, "%s.dat", tmp_base);
        /* Write data down */
        FILE *fin = fopen(tmp_suff, "w");
        offset = buff->start;
        for (int j = 0; j < N1; j++) {
            for (int k = 0; k < N2; k++) {
                fprintf(fin, "% 8.3f ",
                        buff->image[(i + offset) % buff->size][N2 * j + k]);
            }
            fprintf(fin, "\n");
        }
        fclose(fin);
    }
}

void save_2D_ppm(const struct circ_buffer *buff, const char *fname, double lo,
                 double hi, int N1, int N2)
/* Save .ppm and .jpg with a snapshot of firing rates in the 2D network */
{
    char filename[50];
    char command_str[100];
    int status = 0;

    strcpy(filename, fname);

    /* Create the lookup table */
    /* unsigned int rgb_jet[MAXVAL + 1][3]; */
    /* fill_cm_lookuptable_jet(rgb_jet, MAXVAL); */
    unsigned int rgb_bone[MAXVAL + 1][3];
    fill_cm_lookuptable_bone(rgb_bone, MAXVAL);

    size_t offset;
    int pos;
    double normalization;
    char tmp_base[50];
    char tmp_suff[50];

    for (int i = 0; i < N_SNAPSHOTS; i++) {
        sprintf(tmp_base, "%s_frame_%02d", filename, i);
        sprintf(tmp_suff, "%s.ppm", tmp_base);
        /* Write data down */
        FILE *fin = fopen(tmp_suff, "w");
        /* Pgm format */
        fprintf(fin, "P3\n%d %d\n%d\n", N1, N2, MAXVAL);
        offset = buff->start;
        normalization = hi - lo;
        if (normalization == 0.0)
            normalization = 1.0;        /* Avoid division by zero */
        for (int j = 0; j < N1; j++) {
            for (int k = 0; k < N2; k++) {
                pos = (int) rint((double) MAXVAL *
                                 (buff->
                                  image[(i + offset) % buff->size][N2 * j + k] -
                                  lo)
                                 / normalization);
                fprintf(fin, "% 3d % 3d % 3d   ", rgb_bone[pos][0],
                        rgb_bone[pos][1], rgb_bone[pos][2]);
            }
            fprintf(fin, "\n");
        }
        fclose(fin);

        /* WARNING: This is extremely system dependent. I'm assuming
         * "imagemagick" is installed in your system.*/
        sprintf(command_str, "convert %s -rotate -90 %s.jpg", tmp_suff,
                tmp_base);
        status = system(command_str);
        if (status != 0)
            exit(status);

    }
}

void replace_substring(char *pst, const char *str, const char *orig, const char *rep)
{
    char *p;
    static char buffer[128];

    if(!(p = strstr(str, orig)))  /* Is 'orig' even in 'str'? */
        strcpy(pst, str);

    strncpy(buffer, str, p - str); /* Copy characters from 'str' start to 'orig' */
    buffer[p - str] = '\0';

    sprintf(buffer + (p - str), "%s%s", rep, p + strlen(orig));
    strcpy(pst, buffer);
}
