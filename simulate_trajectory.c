#include "ring_rate.h"

int main(int argc, char *argv[])
{
    struct State S;
    S.is_2D = false;

    read_network_parameters(argc, argv, &S);
    initialize_rng();
    initialize_network(&S);
    initialize_circular_buffer(&S);
    read_input_parameters(&S);

    simulate_single_trajectory(&S);

    save_dynamic_variables(&S);
    save_inputs(&S);

    free_network_variables(&S);
    free_input_variables(&S);
    return 0;
}
