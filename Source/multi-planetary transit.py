# import modules
from astropy import constants as const
from astropy import units as u
import matplotlib.pylab as plt
import numpy as np

# import my module
from ProjectModules import nbody

# initialise our system
def main():
    sun = nbody.Body(1.98855E30*u.kg)

    earth = nbody.Body(5.9722E24 * u.kg)
    earth.orbit(sun, 149598023*u.km, 0.0167086, 0.00005)

    venus = nbody.Body(4.8675E24 * u.kg)
    venus.orbit(sun, 0.723332 * const.au, 0.006772, 3.39458)

    moon = nbody.Body(7.342E22 * u.kg)
    moon.orbit(earth, 384299*u.km, 0.0549, 5.145)

    # add to the system
    system = nbody.System()
    system.addBody(sun, venus, earth, moon)

    for x in range(365):
        system.propogate(1*u.day)

    fig=plt.figure()
    ax = fig.add_subplot(projection="3d")
    ax.scatter(*sun.pos)
    for body in system.bodies:
        if body == sun:
            continue
        ax.plot(*np.array(system.positions[body]).T)
    ax.set_xlim(-const.au.value, const.au.value)
    ax.set_ylim(-const.au.value, const.au.value)
    ax.set_zlim(-const.au.value, const.au.value)
    plt.show()

# execute if we are the top level code
if __name__ == "__main__":
    main()
