#include "parameters.h"

/* Provide reasonable defaults */
void init_parameters(struct Parameters *pconfig)
{
    strcpy(pconfig->conf_file, "config");
    strcpy(pconfig->extinputs_file, "inputs");
    strcpy(pconfig->output_rate_file, "rates");
    strcpy(pconfig->output_adaptation_file, "adaptation");
    strcpy(pconfig->output_depression_file, "depression");
    strcpy(pconfig->output_facilitation_file, "facilitation");
    pconfig->verbose = 0;
    pconfig->normalize = 0;
    pconfig->N1 = 128;
    pconfig->N2 = 1;
    pconfig->N = pconfig->N1;
    pconfig->N_layers = 3;
    pconfig->T = 10.;
    pconfig->dt = 2e-3;

    pconfig->sgmoid.max_rate = 1.0;
    pconfig->sgmoid.beta = 3.0;
    pconfig->sgmoid.thr = 1.0;

    /* Adaptation */
    pconfig->negfbck.J_A = 1.0;
    pconfig->negfbck.tau_A = 5.0;
    /* Short term depression and facilitation */
    pconfig->negfbck.depression_f = false;
    pconfig->negfbck.facilitation_f = false;
    pconfig->negfbck.tau_F = 100.0;
    pconfig->negfbck.tau_D = 10.0;
    pconfig->negfbck.U = 1;
    pconfig->negfbck.tau_F = 100.0;
    pconfig->negfbck.tau_D = 10.0;

    pconfig->kernel.J_E = 1.0;
    pconfig->kernel.J_I = 1.0;
    pconfig->kernel.m_E = 1.0;
    pconfig->kernel.m_I = 1.0;

    pconfig->ff_kernel.J_E = 0.1;
    pconfig->ff_kernel.J_I = 0.1;
    pconfig->ff_kernel.m_E = 1.0;
    pconfig->ff_kernel.m_I = 0.0;
}

void show_parameters(void *usercfg)
{
    struct Parameters *pconfig = (struct Parameters *) usercfg;

    printf("Input files:\n");
    printf("  Config loaded from '%s':\n", pconfig->conf_file);
    printf(  "External inputs file: '%s'\n", pconfig->extinputs_file);
    printf("Output files:\n");
    printf("  Rate output filename: '%s'\n", pconfig->output_rate_file);
    printf("  Adaptation output filename: '%s'\n", pconfig->output_adaptation_file);
    printf("  Depression output filename: '%s'\n", pconfig->output_depression_file);
    printf(  "Facilitation output filename: '%s'\n", pconfig->output_facilitation_file);
    printf("\nNumber of neurons    = %d (%d x %d)\n",
           pconfig->N, pconfig->N1, pconfig->N2);
    printf("Number of layers     = %d\n", pconfig->N_layers);
    printf("Time step            = %.1e\n", pconfig->dt);
    printf("Normalize responses? = %d\n", pconfig->normalize);
    printf("Total simulated time = %.1f\n", pconfig->T);
    printf("\nMaximum rate = %.1f\n", pconfig->sgmoid.max_rate);
    printf("Beta         = %.1f\n", pconfig->sgmoid.beta);
    printf("Sigmoid bias = %.1f\n", pconfig->sgmoid.thr);
    printf("\nJ_A          = %.1f\n", pconfig->negfbck.J_A);
    printf("tau_A         = %.1f\n", pconfig->negfbck.tau_A);
    if(pconfig->negfbck.depression_f) {
        printf("Depression: ON\n");
        printf("  tau_D         = %.1f\n", pconfig->negfbck.tau_D);
    } else {
        printf("Depression: OFF\n");
    }
    if(pconfig->negfbck.facilitation_f) {
        printf("Facilitation: ON\n");
        printf("  Utilization   = %.1f\n", pconfig->negfbck.U);
        printf("  tau_F         = %.1f\n", pconfig->negfbck.tau_F);
    } else {
        printf("Facilitation: OFF\n");
    }
    printf("\nLateral connectivity:\n  J_E = %.1f\t", pconfig->kernel.J_E);
    printf("  J_I = %.1f\n", pconfig->kernel.J_I);
    printf("  m_E = %.1f\t", pconfig->kernel.m_E);
    printf("  m_I = %.1f\n", pconfig->kernel.m_I);
    printf("\nFeedforward connectivity:\n  J_E = %.1f\t",
           pconfig->ff_kernel.J_E);
    printf("  J_I = %.1f\n", pconfig->ff_kernel.J_I);
    printf("  m_E = %.1f\t", pconfig->ff_kernel.m_E);
    printf("  m_I = %.1f\n", pconfig->ff_kernel.m_I);
}

int handler(void *usercfg, const char *name, const char *value)
{
    struct Parameters *pconfig = (struct Parameters *) usercfg;

    if (strcasecmp(name, "extinputs_file") == 0) {
        strcpy(pconfig->extinputs_file, value);
    } else if (strcasecmp(name, "output_rate_file") == 0) {
        strcpy(pconfig->output_rate_file, value);
    } else if (strcasecmp(name, "output_adaptation_file") == 0) {
        strcpy(pconfig->output_adaptation_file, value);
    } else if (strcasecmp(name, "output_depression_file") == 0) {
        strcpy(pconfig->output_depression_file, value);
    } else if (strcasecmp(name, "output_facilitation_file") == 0) {
        strcpy(pconfig->output_facilitation_file, value);
    } else if (strcasecmp(name, "normalize") == 0) {
        pconfig->normalize = atoi(value);
    } else if (strcasecmp(name, "verbose") == 0) {
        pconfig->verbose = atoi(value);
    } else if (strcasecmp(name, "N") == 0) {
        pconfig->N = atoi(value);
        pconfig->N1 = pconfig->N;
    } else if (strcasecmp(name, "N1") == 0) {
        pconfig->N1 = atoi(value);
        pconfig->N = pconfig->N1 * pconfig->N2;
    } else if (strcasecmp(name, "N2") == 0) {
        pconfig->N2 = atoi(value);
        pconfig->N = pconfig->N1 * pconfig->N2;
        if (pconfig->N2 > 1)
            strcpy(pconfig->extinputs_file, "inputs_2D");
    } else if (strcasecmp(name, "N_layers") == 0) {
        pconfig->N_layers = atoi(value);
    } else if (strcasecmp(name, "dt") == 0) {
        pconfig->dt = atof(value);
    } else if (strcasecmp(name, "T") == 0) {
        pconfig->T = atof(value);
    } else if (strcasecmp(name, "max_rate") == 0) {
        pconfig->sgmoid.max_rate = atof(value);
    } else if (strcasecmp(name, "beta") == 0) {
        pconfig->sgmoid.beta = atof(value);
    } else if (strcasecmp(name, "thr") == 0) {
        pconfig->sgmoid.thr = atof(value);
    } else if (strcasecmp(name, "U") == 0) {
        pconfig->negfbck.U = atof(value);
    } else if (strcasecmp(name, "depression") == 0) {
        if (atoi(value) == 0) {
            pconfig->negfbck.depression_f = false;
        } else if (atoi(value) == 1) {
            pconfig->negfbck.depression_f = true;
        } else {
            printf("'depression' is a boolean flag. It takes values 0 or 1. Please check your config file\n");
            exit(1);
        }
    } else if (strcasecmp(name, "facilitation") == 0) {
        if (atoi(value) == 0) {
            pconfig->negfbck.facilitation_f = false;
        } else if (atoi(value) == 1) {
            pconfig->negfbck.facilitation_f = true;
        } else {
            printf("'facilitation' is a boolean flag. It takes values 0 or 1. Please check your config file\n");
            exit(1);
        }
    } else if (strcasecmp(name, "J_A") == 0) {
        pconfig->negfbck.J_A = atof(value);
    } else if (strcasecmp(name, "tau_A") == 0) {
        pconfig->negfbck.tau_A = atof(value);
    } else if (strcasecmp(name, "tau_F") == 0) {
        pconfig->negfbck.tau_F = atof(value);
    } else if (strcasecmp(name, "tau_D") == 0) {
        pconfig->negfbck.tau_D = atof(value);
    } else if (strcasecmp(name, "J_rec_E") == 0) {
        pconfig->kernel.J_E = atof(value);
    } else if (strcasecmp(name, "J_rec_I") == 0) {
        pconfig->kernel.J_I = atof(value);
    } else if (strcasecmp(name, "m_rec_E") == 0) {
        pconfig->kernel.m_E = atof(value);
    } else if (strcasecmp(name, "m_rec_I") == 0) {
        pconfig->kernel.m_I = atof(value);
    } else if (strcasecmp(name, "J_ff_E") == 0) {
        pconfig->ff_kernel.J_E = atof(value);
    } else if (strcasecmp(name, "J_ff_I") == 0) {
        pconfig->ff_kernel.J_I = atof(value);
    } else if (strcasecmp(name, "m_ff_E") == 0) {
        pconfig->ff_kernel.m_E = atof(value);
    } else if (strcasecmp(name, "m_ff_I") == 0) {
        pconfig->ff_kernel.m_I = atof(value);
    }
    return 0;
}
