# import modules
from astropy import constants as const
from astropy import units as u
import matplotlib.pylab as plt
import matplotlib.animation as animation
from scipy.fft import rfft, rfftfreq, irfft
from scipy.optimize import curve_fit
from tqdm import trange, tqdm
import numpy as np
import scipy as sp

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

    # add a super earth-like planet in 1:2 resonance
    resonance = 1/2
    p2 = nbody.Body(mass = 4 * 5.9722E24*u.kg, radius=6371*u.km * 4**(1/3), colour = "green", name = "Planet 2")
    p2.keplerian(parent = star, sma = primary_sma * resonance**(3/2), ecc = 0, inc = 0, arg=0)

    # set the system timestep
    system.tstep = int(min(p1.period, p2.period) / 200)

    # set the position of the distant observer
    observer_pos = np.array([D(0), D(((20 * u.lyr).to(u.m)).value), D(0)])

    # propogate the system to find transits
    def find_transits(iterations = 200*200):
        # define our list of transits,
        p1.phi, p1.transit = D(1000), 0
        p1.transit=0
        transit=[]

        # iterate for 100 orbits (roughly)
        for _ in trange(iterations):
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
        # return the found transits
        return np.array([float(t) for t in np.array(transit)])

    # number of steps to integrate
    integration_steps = 100 * 200
    # ensure it's a power of two for fourier fitting
    integration_steps = int(2 ** np.ceil(np.log(integration_steps) / np.log(2)))
    # compute transits for the system with just the star and target planet
    transit_single = find_transits(integration_steps)
    # reset system positions
    system.resetPropogation()
    # compute transits for the system with a perturbing body
    system.addBody(p2)
    transit_full = find_transits(integration_steps)

    # fetch the number of transits in our data
    idxs = min(len(transit_single), len(transit_full))
    transit_full, transit_single = transit_full[:idxs], transit_single[:idxs]



    A = np.vstack([np.ones(idxs), range(idxs)]).T
    c, m = np.linalg.lstsq(A, transit_full, rcond=-1)[0]

    plt.plot(range(idxs), (transit_full-m*np.array(range(idxs))-c))
    plt.show()


    # create an array of orbits to plot against
    orbit_num = np.arange(len(transit_full))

    # find the difference in single vs perturbed transit
    perturbed = transit_single - transit_full

    # create a model to fit the data against
    def linefit(x, a, b):
        return a * x + b
    # fit a straight line to the perturbed value
    popt, pcov = curve_fit(linefit, orbit_num, perturbed)
    # and compute the residuals
    residuals = perturbed - linefit(orbit_num, *popt)
    plt.plot(orbit_num, residuals)
    plt.show()

    # compute the average sampling rate
    duration = max(transit_full[-1], transit_single[-1]) - min(transit_full[0], transit_single[0])
    sampling_rate = len(transit_full) / duration

    # compute the fourier transform of the residuals
    N = int(sampling_rate * duration)
    yf = rfft(residuals)
    xf = rfftfreq(N, 1 / sampling_rate)
    plt.plot(xf, np.abs(yf))
    plt.show()

    # compute the inverse fourier transform
    inv = irfft(yf)
    plt.plot(inv)
    plt.show()



    def sinfunc(x, a, b, c, d, e, f, g):
        return a * np.sin(b*(x - c)) + d * np.sin(e*(x - f)) + g
    params = (1, 0.5, 0.5, 1, 0.5, 0.5, 1)
    popt, pcov = curve_fit(sinfunc, orbit_num, residuals, p0=params)
    plt.scatter(orbit_num,residuals)
    plt.plot(orbit_num, sinfunc(orbit_num, *popt), 'r-',label='Fitted function')
    plt.legend()
    plt.show()

    # plot the two lines
    #plt.plot(transit)
    '''plt.plot(residuals)
    plt.xlabel("Orbit number")
    plt.ylabel("Transit variation (s)")
    plt.title("Transit variation from a simple N-Body simulation")
    plt.tight_layout()
    #plt.savefig("TTVSim2.png")
    plt.show()'''

# execute if we are the top level code
if __name__ == "__main__":
    main()
