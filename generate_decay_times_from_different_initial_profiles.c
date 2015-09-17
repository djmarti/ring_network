#include "ring_rate.h"

int main(int argc, char *argv[])
{
    struct State S;
    S.is_2D = false;
    double m;
    char filename[50];

    read_network_parameters(argc, argv, &S);
    initialize_network(&S);
    initialize_circular_buffer(&S);
    double r_max = 0.2;
    for (int i = 0; i < 4; i++){
        m = (double) i * 2.0;
        initialize_rate_profile(&S, r_max, m);
        simulate_single_trajectory(&S);
        /* Save the data */
        sprintf(filename, "rate_decay_max%3.1f_%02d.dat", r_max, (int) m);
        printf("Save simulation for m=-%3.1f and r_max=%3.1f in %s.\n", 
                m, r_max, filename);
        save_dynamic_variables(&S);
    }
    free_network_variables(&S);
    free_input_variables(&S);
    return 0;
}
