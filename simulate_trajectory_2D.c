#include "ring_rate.h"

int main(int argc, char *argv[])
{
    struct State S;

    read_network_parameters(argc, argv, &S);
    read_input_parameters(&S);
    initialize_network_2D(&S);
    initialize_circular_buffer(&S);

    simulate_single_trajectory_2D(&S);

    save_snapshots_firing_activity_2D(&S);
    /* save_inputs(&S); */

    free_network_variables(&S);
    free_input_variables(&S);
    return 0;
}
