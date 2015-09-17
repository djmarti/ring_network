#include "ext_inputs.h"

struct pulse *new_pulse(int N, double t, double d, double q, double input,
                        double m)
{
    struct pulse *p;

    p = (struct pulse *) malloc(sizeof(struct pulse));
    p->time_onset = t;
    p->duration = d;
    p->theta = q; /* q is given in degrees! */
    p->Imax = input;
    p->concentration = m;
    p->switched_on = 0;
    p->precomputed_value = precompute_vector_pulses(N, p);
    p->next = NULL;
    return p;
}

struct pulse_2D *new_pulse_2D(int N1, int N2, double t, double d, double q1,
                              double q2, double input, double m)
{
    struct pulse_2D *p;

    p = (struct pulse_2D *) malloc(sizeof(struct pulse_2D));
    p->time_onset = t;
    p->duration = d;
    p->theta1 = q1;
    p->theta2 = q2;
    p->Imax = input;
    p->concentration = m;
    p->switched_on = 0;
    p->precomputed_value = precompute_vector_pulses_2D(N1, N2, p);
    p->next = NULL;
    return p;
}

double *precompute_vector_pulses(int N, struct pulse *plse)
{
    double *v = malloc(N * sizeof(double));
    double angle[N];
    double theta_rad = plse->theta * M_PI / 180.;

    for (int j = 0; j < N; j++) {
        /* TODO avoid computing angle each time a new pulse is added */
        angle[j] = M_PI * (-1. / 2. + (double) j / (double) N);
        v[j] = plse->Imax * exp(plse->concentration
                                * (cos(2.0 * (angle[j] - theta_rad)) - 1));
    }
    return v;
}

double *precompute_vector_pulses_2D(int N1, int N2, struct pulse_2D *plse)
{
    double *v = malloc(N1 * N2 * sizeof(double));
    double angle1[N1];
    double angle2[N2];
    double theta1_rad = plse->theta1 * M_PI / 180.;
    double theta2_rad = plse->theta2 * M_PI / 180.;
    double dist;

    /* TODO avoid computing angle each time a new pulse is added */
    for (int i = 0; i < N1; i++)
        angle1[i] = M_PI * (-1. / 2. + (double) i / (double) N1);
    for (int i = 0; i < N2; i++)
        angle2[i] = M_PI * (-1. / 2. + (double) i / (double) N2);

    for (int i = 0; i < N1; i++) {
        for (int j = 0; j < N2; j++) {
            dist =
                distance_on_torus(angle1[i], angle2[j], theta1_rad, theta2_rad);
            v[N2 * i + j] =
                plse->Imax * exp(plse->concentration * (cos(2.0 * dist) - 1));
        }
    }
    return v;
}

struct pulse *append_pulse(struct pulse *plse, struct pulse *newplse)
{
    struct pulse *p;

    /* this is the first item in the list */
    if (plse == NULL) {
        return newplse;
    } else {
        /* walk to the end of the list */
        for (p = plse; p->next != NULL; p = p->next)
            /* do nothing */ ;
        p->next = newplse;
        return plse;
    }
}

struct pulse_2D *append_pulse_2D(struct pulse_2D *plse,
                                 struct pulse_2D *newplse)
{
    struct pulse_2D *p;

    /* this is the first item in the list */
    if (plse == NULL) {
        return newplse;
    } else {
        /* walk to the end of the list */
        for (p = plse; p->next != NULL; p = p->next)
            /* do nothing */ ;
        p->next = newplse;
        return plse;
    }
}

void print_pulses(struct pulse *plse)
{
    printf("\nExternal inputs...\n");
    for (; plse != NULL; plse = plse->next) {
        printf("% 8.2f % 8.2f % 8.2f % 8.2f % 8.2f\n",
               plse->time_onset, plse->duration, plse->theta,
               plse->Imax, plse->concentration);
    }
}

void print_pulses_2D(struct pulse_2D *plse)
{
    printf("\nExternal inputs...\n");
    for (; plse != NULL; plse = plse->next) {
        printf("% 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f\n",
               plse->time_onset, plse->duration,
               plse->theta1, plse->theta2, plse->Imax, plse->concentration);
    }
}


double max_amplitude_pulse(struct pulse *plse)
{
    double tmp = 0;
    for (; plse != NULL; plse = plse->next) {
        if (plse->Imax > tmp) {
            tmp = plse->Imax;
        }
    }
    return tmp;
}

double max_amplitude_pulse_2D(struct pulse_2D *plse)
{
    double tmp = 0;
    for (; plse != NULL; plse = plse->next) {
        if (plse->Imax > tmp) {
            tmp = plse->Imax;
        }
    }
    return tmp;
}

void free_extinputs(struct pulse *plse)
{
    struct pulse *next;
    for (; plse != NULL; plse = next) {
        next = plse->next;
        free(plse->precomputed_value);
        free(plse);
    }
}

void free_extinputs_2D(struct pulse_2D *plse)
{
    struct pulse_2D *next;
    for (; plse != NULL; plse = next) {
        next = plse->next;
        free(plse->precomputed_value);
        free(plse);
    }
}
