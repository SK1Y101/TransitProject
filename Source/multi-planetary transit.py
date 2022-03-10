# import modules
from astropy import constants as const
from astropy import units as u
import matplotlib.pylab as plt
import matplotlib.animation as animation
from tqdm import trange, tqdm
import numpy as np

from random import random
from decimal import *
D = Decimal

getcontext().prec = 50

# import my module
from ProjectModules import nbody, trig_funcs

# initialise our system
def main():
    # fetch our system
    system = nbody.System()

    # add the star
    star = nbody.Body(mass = 1.98855E30*u.kg, radius=695700*u.km, colour = "yellow", name = "Star")
    system.addBody(star)

    # add an earth-like planet at 0.02 AU
    primary_sma = 0.02*const.au
    p1 = nbody.Body(mass = 5.9722E24*u.kg, radius=6371*u.km, colour = "blue", name = "Planet 1")
    p1.keplerian(parent = star, sma = primary_sma, ecc = 0, inc = 0)
    system.addBody(p1)

    # add a super earth-like planet in 3:1 resonance
    resonance = 3/2
    p2 = nbody.Body(mass = 10 * 5.9722E24*u.kg, radius=6371*u.km * 10**(1/3), colour = "green", name = "Planet 2")
    p2.keplerian(parent = star, sma = primary_sma * resonance**(3/2), ecc = 0, inc = 0, arg=0)
    system.addBody(p2)

    # set the system timestep
    system.tstep = int(min(p1.period, p2.period) / 400)

    # set the position of the distant observer
    observer_pos = np.array([D(0), D(((20 * u.lyr).to(u.m)).value), D(0)])

    # define our list of transits,
    p1.phi, p1.transit = D(1000), 0
    p1.transit=0
    transit=[]

    # iterate for 60 orbits (roughly)
    for _ in trange(400 * 20):
        # propogate one timestep
        system.propogate(system.tstep, keepBarycentreConstant=True)

        # loop over all planets
        for body in system.bodies:
            # compute distance to body
            body.dist = (observer_pos - body.pos)
            body.dist_mag = np.linalg.norm(body.dist)

            # compute the angular radius
            body.theta = trig_funcs.arctan2(body.radius, body.dist_mag)

        # deal with the transit of the target planet

        # get the positions of the planet and star
        s_vec, s_mag = star.dist, star.dist_mag
        p_vec, p_mag = p1.dist, p1.dist_mag

        # dot and magnitude product
        dot_prod = sum(p_vec * s_vec)
        mag_prod = p_mag * s_mag

        # angle between the planet and star, as seen from earth
        cosphi = dot_prod / mag_prod
        phi = D(trig_funcs.arccos(cosphi))

        # don't do anything if we're not in front of the star
        if p1.dist_mag > star.dist_mag:
            p1.phi = phi
            continue

        # if phi bounces near zero, we had the midtransit occur
        if phi < p1.phi:
            p1.transit = 1
        elif (phi > p1.phi) and (p1.transit):
            p1.transit=0
            # linearly interpolate the time to get the midtransit point
            start, end = p1.phi, -phi
            slope = (end-start) / (system.tstep)
            time = end / slope + system.time
            # append to our list
            transit.append(time)

        # update our value for phi
        p1.phi = phi

    # get the transits
    transit = np.array([float(t) for t in np.array(transit)])
    # offset the begining
    transit -= transit[0]
    # compute the number of orbits for each transit
    orbit_num = np.round_(transit / float(p1.period))
    # fit a line to the transit timing
    p, residuals, _, _, _ = np.polyfit(orbit_num, transit, 1, full=True)
    print(p, residuals, residuals / (len(orbit_num) - 3))
    # and plot the residuals of that line
    plt.plot(transit - np.poly1d(p)(orbit_num))
    plt.show()

# execute if we are the top level code
if __name__ == "__main__":
    main()
