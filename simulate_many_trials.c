#include "ring_rate.h"

int main(int argc, char *argv[])
{
    struct State S;
    S.is_2D = false;
    int N_trials = 200;
    double time_formation;
    double angle_location;
    char filename[100];

    read_network_parameters(argc, argv, &S);
    initialize_rng();
    initialize_network(&S);
    initialize_circular_buffer(&S);

    double m_s = 5.0;
    double dispersion;
    double duration = 0.10;
    double intensity = 0.2;
    double location[1] = {0};

    FILE *fin;

    for (int i=0; i < 12; i++) {
        dispersion = 15 + 5 * i;
        sprintf(filename,
                "bump_formation_times_ms_%02d_duration_%4.2f_i_%3.1f_dispersion_%02d", 
                (int) m_s, duration, intensity, (int)dispersion);
        fin = fopen(filename, "w");
        for (int j=0; j < N_trials; j++) {
            printf("%d\n", j);
            generate_pulses(&S, 1, location, dispersion, duration, intensity, m_s);
            simulate_single_trajectory_self_stop(&S, &time_formation, &angle_location);
            fprintf(fin, "% 8.3f % 10.3f\n", time_formation, angle_location);
            printf("% 8.3f % 10.3f\n", time_formation, angle_location);
            reset_rates_and_inputs(&S);
            free_extinputs(S.pulse_list);
        }
        fclose(fin);
    }
    free_network_variables(&S);
    fclose(fin);
    return 0;
}
