''' Handles the creation and execution of different TTV Models '''
# Python modules
import pymc3 as pm
import numpy as np

# my modules
from . import *

def model1(t, target, *planets):
    ''' Model for a transiting planet exterior to 'n' coplanar, non-interacting, zero eccentricity perturbing planets. '''
    print(planets)
    # position of the barycentre
    rb = -ax * mux * np.sin(2*np.pi*(t -t0x) / px)
    # TTV due to the position
    TTV = (p1 / (2*np.pi*a1)) * rb

def dfToParams():

#TTV1 = -sma_1 * mu_1 * p_1 * np.sin(2*np.pi*(n*p_2 - t) / p_1) / (2*np.pi*sma_2)
model1(np.arange(9, 20), None)
# where n is the transit number, and t is the time of the transit, such that n*p2 gives the position of the outer planet

print(_planetVars_)

with pm.Model() as transitModel:
    # data
    data = pm.Data("data", y)
    # parameters
    mag = pm.Uniform("mag", lower=0, upper=1000)
    per = pm.Uniform("per", lower=0, upper=1000)
    phase = pm.Normal("phase", mu=np.pi, sigma=1.5)
    # expected outcome
    mu = model(x, mag, per, phase) #mag * pm.math.sin(x * (2 * np.pi / per) + phase)
    # observation
    obs = pm.Normal("y", mu=mu, sigma=1, observed=y)
