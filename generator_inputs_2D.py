#!/usr/bin/env python
'''This script generates an 'input' data file for the 2D network. It specifies
the time onset, duration, location (theta_1, theta_2), concentration (assumed
isotropic), and intensity of each pulse in the input stream. The data could
have been entered by hand, but the script comes in handy when the input stream
contains more than 5 pulses.'''

import numpy as np

# Peaks of stimulus features
#mus = [ (-45, 30), (45,-50) ]      # Location of each peak
mus = [ (1,-2) ]
n_peaks = len(mus)
dispersion = 20.    # concentration parameter of the distribution of input
                    # centers in a peak (assumed equal for all peaks)
n_pulses_peak = 20
Imax = 0.5
duration = 0.2
onset = 0.2
concentration = 4       # Concentration parameter of each single pulse
interpulse_sep = 0.5 * duration # silent period between pulses

A = np.zeros((n_peaks * n_pulses_peak, 6))

# Fill the entries
# A[:,2] = mus[0] * ones (2 * n_pulses_peak)
# Split first the pulses into peaks, we'll sort by time onset later
# Iterate over peaks
for i in range(n_peaks):
    A[:n_pulses_peak,2] = np.degrees(
        np.random.vonmises(mu=np.radians(mus[i][0] * 2), kappa=dispersion,
                           size=(n_pulses_peak))) / 2
    A[:n_pulses_peak,3] = np.degrees(
        np.random.vonmises(mu=np.radians(mus[i][1] * 2), kappa=dispersion,
                           size=(n_pulses_peak))) / 2
    A[:n_pulses_peak,0] = onset + i * (duration + interpulse_sep) + \
                    np.arange(0, n_peaks * (duration + interpulse_sep)
                              * n_pulses_peak,
                              n_peaks * (duration + interpulse_sep))
# Common values for both peaks
A[:,1] = duration * np.ones(n_peaks * n_pulses_peak)
A[:,4] = Imax * np.ones(n_peaks * n_pulses_peak)
A[:,5] = concentration * np.ones(n_peaks * n_pulses_peak)

A = A[np.argsort(A[:,0])]

init_string = '''#
# Configuration of external inputs
#
# Distribution of location parameter:
'''
prob_peaks = 1.0 / n_peaks
for i, m in enumerate(mus):
    init_string += '#   Peak %d (%4.2f): mu = (%4.1f, %4.1f), dispersion = %4.1f\n' %  \
        (i, prob_peaks, m[0], m[1], dispersion)

init_string += '''#
# Single pulses
#   Each line specifies the characteristics of a particular pulse of stimulation,
#   and is structured in 5 different fields separated by spaces. The order of the
#   fields is:
#       1. time onset of the pulse,
#       2. duration of the pulse
#       3. angle_x (in degrees) where the pulse is applied,
#       4. angle_y (in degrees)
#       5. intensity of the pulse,
#       6. concentration of the pulse (in a von Mises distribution)
#
# time onset  duration  theta_x  theta_y   Imax  concentration
# ------------------------------------------------------------
'''

numbers = ['one', 'two', 'three', 'four']
fname = 'i_2D_%s_peaks.txt' % numbers[n_peaks-1]
fin = open(fname, 'w')
fin.write(init_string)
for i in range(len(A)):
    fin.write('% 10.1f % 9.1f % 8.1f % 8.1f % 8.1f % 9.1f\n' % (A[i,0], A[i,1],
                                                                A[i,2], A[i,3],
                                                                A[i,4], A[i,5]))

print 'Input data saved in ' + fname + '.'
fin.close()
