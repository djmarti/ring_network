#!/usr/bin/env python
"""
This script generates an 'input' data file specifying the time onset,
duration, location, concentration, and intensity of each pulse in the input
stream. The data could have been entered by hand, but the script comes in handy
when the input stream contains more than 5 pulses.
"""

import numpy as np
from numpy import pi
from random import sample

# Peaks of stimulus features

mus = {
    'one_peak': [0],
    'two_peaks': [-pi/4, pi/4],
    'two_closer_peaks': [-pi/10, pi/10],
    'two_broad_peaks': [-pi/4, pi/4],
    'three_peaks': [-pi/3, 0, pi/3],
    'four_peaks': [-3*pi/8, -pi/8, pi/8, 3*pi/8]
}

dispersion = {
    'one_peak': 50,
    'two_peaks': 50,
    'two_closer_peaks': 50,
    'two_broad_peaks': 5,
    'three_peaks': 50,
    'four_peaks': 50
}

T_max = 20.
duration = 0.1  # 0.25 for three_random_peaks and four_random_peaks
gap_between_pulses = 1e-2
onset = 0.1     # 0.4 for  three_random_peaks and four_random_peaks

I_s = 0.2
concentration = 10
n_params = 5

for c in mus:
    n_peaks = len(mus[c])
    n_pulses_peak = int(
        (T_max - onset) /
        ((duration + gap_between_pulses) * n_peaks))
    n_rows = n_peaks * n_pulses_peak
    A = np.zeros((n_rows, n_params))

    # Fill the entries
    for i, mu in enumerate(mus[c]):
        init = i * n_pulses_peak
        end = (i + 1) * n_pulses_peak
        A[init:end, 2] = np.degrees(
            np.random.vonmises(mu=mu * 2,
                               kappa=dispersion[c],
                               size=(n_pulses_peak))) / 2
    indices = np.arange(n_rows).tolist()

    # Uncomment this if you want a really random process
    indices_shff = np.array(sample(indices, n_rows))

    # This will generate an alternating sequence:
    # indices_resorted = zeros_like(indices)
    # for i in arange(n_peaks):
    #     init = i * n_pulses_peak
    #     end = (i + 1) * n_pulses_peak
    #     indices_resorted[init:end] = arange(i, n_rows, n_peaks)

    # Common values for both peaks
    A[:, 0] = onset + indices_shff * (duration + gap_between_pulses)
    # A[:, 0] = onset + indices_resorted * (duration + gap_between_pulses)
    A[:, 1] = duration * np.ones(n_rows)
    A[:, 3] = I_s * np.ones(n_rows)
    A[:, 4] = concentration * np.ones(n_rows)

    A = A[np.argsort(A[:, 0])]

    ## ----------------------------------------------------------
    ## Apply a slow drift at the end of
    ## the sequence
    ## n_plses = 180
    ## B = zeros((n_plses, n_params))
    ## speed = 0.6 # degree / s
    ## t_drift_starts = A[-1,0] + A[-1,1] + 2 * gap_between_pulses
    ## A[-1,2] = mus[0] - 0.10
    ## duration = 0.2
    ## B[:,0] = t_drift_starts + arange(n_plses) * duration # (duration + gap_between_pulses)
    ## B[:,1] = duration * ones(n_plses)
    ## B[:,2] = mus[0] + speed * (B[:,0] - t_drift_starts)
    ## B[:,3] = I_s * ones(n_plses)
    ## B[:,4] = concentration * ones(n_plses)
    ##
    ## A = vstack([A,B])

    init_string = """#
# Configuration of external inputs
#
# Distribution of location parameter:
"""
    prob_peaks = 1.0 / len(mus[c])
    for i, m in enumerate(mus[c]):
        init_string += ("#   Peak %d (%4.2f): mu = %4.1f, dispersion = %4.1f\n"
                        % (i, prob_peaks, np.degrees(m), dispersion[c]))

    init_string += """#
# Single pulses
#   Each line specifies the characteristics of a particular pulse of stimulation,
#   and is structured in 5 different fields, separated by spaces. The order of the
#   fields is:
#
#       1. time onset of the pulse,
#       2. duration of the pulse
#       3. angle (in degrees) where the pulse is applied,
#       4. intensity of the pulse,
#       5. concentration of the pulse (in a von Mises distribution)
#
# time onset    duration      theta     I_s   concentration
# -----------------------------------------------------------
"""

    fname = "i_%s_really_shuffled" % c
    #fname = "i_%s_peak_with_drift" % numbers[n_peaks-1]
    #fname = "test"
    with open(fname, "w") as fin:
        fin.write(init_string)
        for i in range(len(A)):
            fin.write("% 12.2f % 11.2f % 10.2f % 8.2f % 15.2f\n" % (
                A[i, 0], A[i, 1], A[i, 2], A[i, 3], A[i, 4]))

        print(f"Input data saved in '{fname}'.")
