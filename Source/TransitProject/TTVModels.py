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

def dfToParams(df):
