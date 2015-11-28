# ring_network
Code for the network of rate units used in the article ['Dynamics of feature categorization'](http://www.mitpressjournals.org/doi/abs/10.1162/NECO_a_00383) by Daniel Martí and John Rinzel in Neural Computation, January 2013, vol. 25, No. 1, pp 1-45. The repository contains tools to stimulate the network with different types of spatio-temporal input patterns.

## Installation and compilation
Just get the code and compile the code with either `make` or [`scons`](http://www.scons.org) in the `ring_network` folder:
```shell
git clone https://github.com/djmarti/ring_network
cd ring_network
make  # or scons
```

## Basic usage
The main executable is `simulate_trajectory`. It simulates the time evolution of the activity of the ring network and stores the resulting activity in several datafiles, one per dynamical variable (firing rate and, if present, adaptation, facilitation, or depression). These activities are stored both in text format and in bitmap format.

A typical run of the executable would be
```shell
./simulate_trajectory -c config_pulse_sequences -x i_single_pulse
```
where the file `config_pulse_sequences` is a configuration file specifying the network and synaptic parameters (see below), and `i_single_pulse` is an example of input descriptor, a file that describes the input pattern fed into the network (see below).

### Configuration files
A standard config file looks like the following:
```shell
# Output files
# ------------
extinputs_file = "inputs"   # Default filename for of input descriptor
output_rate_file = "rates"  # Prefix for the firing rate filename
output_adaptation_file = "adaptation"      # and so on
output_facilitation_file = "facilitation"
normalize = 0

# Simulation parameters
# ---------------------
dt = 1e-2  # Time step for Euler method
T = 10.0   # Total simulated time (time units)

# Parameters sigmoid
# ------------------
N = 512         # Number of neurons
N_layers = 1    # Number of layers
max_rate = 1.0  # Maximum rate
beta = 3.0      # Temperature parameter
sigm_thr = 1.0  # Bias / soft threshold of sigmoid

# Negative feedback parameters
# ----------------------------
J_A = 0.0        # Adaptation amplitude (0: turned off)
tau_A = 10       # Adaptation timescale
facilitation = 0 # Facilitation amplitude (0: turned off)
depression = 0   # Depression amplitude (0: turned off)


# Connectivity parameters
# -----------------------
# Recurrent connections
J_rec_E = 3.5    # Coefficient excitatory profile
J_rec_I = 3.5    # Coefficient inhibitory profile
m_rec_E = 100    # Concentration parameter exc. profile
m_rec_I = 1      # Concentration parameter inh. profile

# Feedforward connections (irrelevant if N_layers = 1)
J_ff_E = 5.0     # Coefficient excitatory profile
J_ff_I = 0.2     # Coefficient inhibitory profile
m_ff_E = 10      # Concentration parameter exc. profile 
m_ff_I = 1       # Concentration parameter inh. profile
```
The example should be self-explanatory.

### Input descriptors
This type of file describes how the network is stimulated. A self-descriptive example of an input descriptor is the following:
```
#
# Configuration of external inputs 
#
# Each line specifies the characteristics of a particular pulse of stimulation,
# and is structured in 5 different fields separated by spaces. The order of the
# fields is: 
#     1. time onset of the pulse, 
#     2. duration of the pulse
#     3. angle (in degrees) where the pulse is applied, 
#     4. intensity of the pulse,
#     5. concentration of the pulse (in a von Mises distribution)
#
# time onset    duration      theta     Imax   concentration
# -----------------------------------------------------------
1.0              1.0        -45.0      1.0       10.0
4.0              0.5         10.5      0.2       10.0
```

## License
This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.
