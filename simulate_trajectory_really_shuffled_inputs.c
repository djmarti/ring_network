
#include "ring_rate.h"

int main(int argc, char *argv[])
{
    struct State S;
    S.is_2D = false;

    read_network_parameters(argc, argv, &S);
    initialize_network(&S);
    initialize_circular_buffer(&S);
    initialize_rng();

    /*const int N_categories = 2;*/
    double m_s = 10.0;
    double dispersion = 50.0;
    double duration = 0.1;
    double intensity = 0.2;
    double location[2] = {-M_PI/4, M_PI/4}; /* in radians */

    generate_pulses(&S, 2, location, dispersion, duration, intensity, m_s);
    simulate_single_trajectory(&S);
    save_firing_activity(&S, S.pconfig.output_file);
    save_inputs(&S);

    free_network_variables(&S);
    free_input_variables(&S);
    return 0;
}
