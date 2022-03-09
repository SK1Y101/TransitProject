# import modules
from astropy import constants as const
from astropy import units as u
import matplotlib.pylab as plt
import matplotlib.animation as animation
from tqdm import trange, tqdm
import numpy as np

from random import random

# import my module
from ProjectModules import nbody

# initialise our system
def main():
    # fetch our system
    system = nbody.System()

    # add the sun
    sun = nbody.Body(mass = 1.98855E30*u.kg, colour = "yellow", name = "Sun")
    system.addBody(sun)

    # add the earth
    earth = nbody.Body(mass = 5.9722E24*u.kg, colour = "blue", name = "Earth")
    earth.keplerian(parent = sun, sma = 149598023*u.km, ecc = 0.0167086, inc = 0.00005)
    system.addBody(earth)

    # add venus
    venus = nbody.Body(mass = 4.8675E24 * u.kg, colour = "brown", name = "Venus")
    venus.keplerian(parent = sun, sma = 0.723332*const.au, ecc = 0.006772, inc = 3.39458)
    system.addBody(venus)

    # set the system timestep
    system.tstep = 1*u.day

    system.gravity_Func = system.do_gravity_vectorized

    # plot the positions of the system over time
    system.plot(interval=1, keepBarycentreConstant=True)

# execute if we are the top level code
if __name__ == "__main__":
    main()
