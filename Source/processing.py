# python modules

# my modules
import TransitProject as tp

'''
Goals: Create a simple extensible model that can take system parameters and output a TTV sinusoid.
       Compare model to simulated TTV to verify its accuracy (Minimise simulation time where possible).
       Provide model in a form that can be used for MCMC fitting from real data.

       Take observed midtransit times and fit various curves to it.
       - ie: the model developed above, simple linear, polyfit, standard sine
       Determine whether the extensible model found is a better fit than others (ie: better than a line).

       If the model fit, find the additional parameters that best fit the data:
       - ie: fit for perturbing planetary period, for example.
       - can be found from fourier analysis, or from the MCMC model parameters probably.
'''

# model:
# required inputs:
# orbital elements of transiting planet
# physical parameters of transiting planet and parent star
# variable inputs:
# physical parameters and orbital elements of perturbing planet(s)

# output:
# sinusoidal curve with correct amplitude and frequency, as seen from simulation.
