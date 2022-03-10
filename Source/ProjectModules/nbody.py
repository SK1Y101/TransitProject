# solve the equations of motion numerically
# propogate the position of n-bodies in a system
# expose the ability to change the numerical method

# each body is derived from a single class
# the body class needs to store: mass, position, velocity
# it would be nice if it could also store: radius, luminosity, parent

# import python modules
from astropy import constants as const
from astropy import units as u
from tqdm import trange, tqdm
import numpy as np

from . import trig_funcs

from decimal import *
D = Decimal

# graphing modules
import matplotlib.animation as animation
import matplotlib.pylab as plt
from matplotlib import gridspec

# standard helper class for functions
class __Helper__:
    # ensure that a quantity has the correct unit type
    def asUnit(self, value, unit):
        ''' Take a value and unit, and ensure they are compatible.
            value: the quantity to check, can be python int/float, astropy quantity, or a list of the above.
            unit: the astropy unit to convert to.
            nparray: Whether the output should be converted to a numpy array with units.'''
        # ensure the unit given is a unit
        if not isinstance(unit, (u.core.Unit, u.core.IrreducibleUnit, u.core.CompositeUnit)):
            raise Exception("Provided unit is not an astropy unit quantity, given {} instead.".format(unit))
        # check if the value is a list of values
        if isinstance(value, list):
            # convert each value individually
            return [self.asUnit(val, unit) for val in value]
        # otherwise, check if this value has a unit already
        if isinstance(value, u.quantity.Quantity):
            # convert
            return value.to(unit)
        # if value is just a number
        return value * unit

    # convert a set of quantities with the same units to a numpy array
    def toNpArray(self, value):
        ''' convert the given value to a numpy array with the same units '''
        # fetch the units from the first value in the list
        if isinstance(value, list):
            unit = value[0].unit
        else:
            unit = value.unit
        # strip all units from the values
        value = [val.value for val in value]
        # and reconstruct as a np array
        return np.array(value) * unit

    def indefinite(self):
        while True:
            yield

# system class
class System(__Helper__):
    # initialisation
    def __init__(self):
        self.bodies = []
        self.propogator = self.euler
        self.grav_Func = self.do_gravity_vectorized
        self.positions = {"barycentre":[]}
        self.totalmass = D(0)
        self.largest_apo = D(0)
        self.barycentre = []
        self.time = D(0)

    # reset the positions of all bodies in the system
    def resetPropogation(self):
        for body in self.bodies:
            body.pos = np.copy(body._pos)
            body.vel = np.copy(body._vel)

            self.positions[body] = [body.pos]
        self.time = D(0)

    # start an updating plot of the system
    def plot(self, updatefunc=None, frames=None, interval=10, keepBarycentreConstant=True):
        # if we weren't given a function, use the inbuild
        updatefunc = updatefunc if updatefunc else self.plotpos
        # start a nice bar so we can see progress
        with tqdm() as pbar:
            # create the figure
            fig=plt.figure(figsize=(20,10))
            fig.tight_layout()
            spec = gridspec.GridSpec(ncols=1, nrows=2, height_ratios=[4, 1])
            # and axes
            ax = fig.add_subplot(spec[0], projection="3d")
            ax2 = fig.add_subplot(spec[1])
            # animation function
            def animate(i):
                # clear the graph
                ax.clear()
                ax2.clear()
                # run our function to update the graph
                updatefunc([ax, ax2], i, keepBarycentreConstant)
                # and update our bar too
                pbar.update(1)
            # execute the animation
            ani = animation.FuncAnimation(fig, animate, frames=frames, interval=interval)
            plt.show()

    # plot the position of each body over time
    def plotpos(self, ax, i, keepBarycentreConstant=True):
        # propogate one timestep
        self.propogate(self.tstep, keepBarycentreConstant)
        # for each body in the system
        for body in self.bodies:
            # plot the previous position of the body
            ax.plot(*self.positions[body].T, c=body.colour, lw=0.5)
            # plot the current position of the body
            ax.scatter(*body.pos, c=body.colour)

    # add bodies to the system
    def addBody(self, *bodies):
        # iterate on each body
        for body in bodies:
            # ensure we have a body object
            try:
                if body.isABody:
                    pass
                # append to our body list
                self.bodies.append(body)
                self.positions[body] = np.copy(body.pos)
                self.totalmass += body.mass
                # compute the apoapsis of this body, and add to our largest function
                apo = (1 + body.ecc) * body.sma
                self.largest_apo = max(self.largest_apo, apo)
            except:
                raise Exception("Object given is not a body: {}.".format(body))

    # compute gravitational forces between bodies
    def do_gravity(self, offset=1*u.m):
        # ensure our offset is valid
        offset = D(self.asUnit(offset, u.m).value)
        a_ = []
        # iterate on all the bodies
        for body1 in self.bodies:
            # reset acceleration from previous step
            a = np.array([D(0), D(0), D(0)])
            # itterate on all the bodies again
            for body2 in self.bodies:
                # distance between
                dist = body1.pos - body2.pos
                dist_mag = np.linalg.norm(dist)
                # compute total force, using the offset to avoid infinities
                F = -D(const.G.value) * body2.mass * body1.mass * dist / (dist_mag + offset)**3
                #compute force in each direction.
                a += F / body1.mass
            # append this bodies acceleration
            a_.append(a)
        # return the acceleration
        return np.array(a_)

    # compute the gravitational forces between bodies without using loops
    def do_gravity_vectorized(self, offset=1*u.m):
        # ensure our offset is valid
        offset = D(self.asUnit(offset, u.m).value)

        # mass and position of all bodies
        mass = np.vstack([body.mass for body in self.bodies])
        poss = np.vstack([body.pos for body in self.bodies])

        # compute m1 * m2
        mass_mat = (mass.reshape((1, -1, 1)) * mass.reshape((-1, 1, 1)).T).T

        # compute the distances between bodies in each axis
        dist = poss.reshape((1, -1, 3)) - poss.reshape((-1, 1, 3))

        # compute the magnitude of the distance between bodies
        dist_mag = np.linalg.norm(dist, axis=2)

        # also remove any zeros
        dist_mag[dist_mag == D(0)] = offset

        # compute the gravitational forces for each body in each direction
        F = D(const.G.value) * mass_mat * dist / np.expand_dims(dist_mag, 2)**3

        # convert to 3d-acceleration
        a = F.sum(axis=1) / mass.reshape(-1, 1)

        return a

    # propogator
    def propogate(self, step=1, keepBarycentreConstant=True):
        # ensure the step is a number
        step = D(self.asUnit(step, u.s).value)

        # compute the gravitational acceleration between bodies
        acc = self.grav_Func(offset = 1*u.m)

        # set the barycentre position to nothing
        self.barycentre = np.array([D(0), D(0), D(0)])

        # iterate on all bodies
        i=0
        for body in self.bodies:
            # apply the chosen propogator
            self.propogator(body, acc[i], step)
            i+=1
            # update historical positions
            self.positions[body] = np.vstack((self.positions[body], np.copy(body.pos)))
            # and add our barycentre to the mix
            self.barycentre += body.pos * (body.mass / self.totalmass)

        # if we are keeping the location of the barycentre fixed
        if keepBarycentreConstant:
            for body in self.bodies:
                body.pos -= self.barycentre
            self.barycentre -= self.barycentre

        # increase the system time
        self.time += step

    # propogator types

    # euler method, applies simple suvat to compute the bodies new position and velocity
    def euler(self, body, a, step):
        # y_{n+1} = y_n + h f(t_n, y_n)
        # ensure the step is a number
        step = D(self.asUnit(step, u.s).value)
        # compute velocity
        body.vel += a * step
        # compute position
        body.pos += body.vel * step

# body class
class Body(__Helper__):
    # initialisation
    def __init__(self, mass, radius=1000, name="Body", colour="black"):
        ''' Initialise the body class.
            mass: the mass of this body. This value will be converted to kilograms
            name: the name of this body for graphing.
            colour: the colour of body for graphing.'''
        self.isABody = True
        # graphing parameters
        self.name = name
        self.colour = colour
        # physical elements
        self.mass = D(self.asUnit(mass, u.kg).value)
        self.radius = D(self.asUnit(radius, u.m).value)
        # cartesian elements
        self.pos = np.array((D(0), D(0), D(0)))
        self.vel = np.array((D(0), D(0), D(0)))
        self._pos = np.array((D(0), D(0), D(0)))
        self._vel = np.array((D(0), D(0), D(0)))
        # keplerian elements
        self.sma = D(0)
        self.ecc = D(0)
        self.inc = D(0)
        self.lan = D(0)
        self.arg = D(0)
        self.tru = D(0)

    # set a body to be orbiting another
    def keplerian(self, parent, sma, ecc, inc, lan=0, arg=0, tru=0):
        ''' Take the orbital elements around the body, and store the position and velocity in the body
            parent: the body object that this orbits around.

            The following define the orbit around the parent:
            sma: Semi-major axis, ecc: Eccentricity, inc: Inclination,
            lan: Longitude of ascending node, arg: Argument of periapse
            tru: True anomaly'''
        # set the body elements
        sma = D(self.asUnit(sma, u.m).value)
        ecc = D(self.asUnit(ecc, u.m / u.m).value)
        inc, lan, arg, tru = [D(self.asUnit(val, u.radian).value) for val in [inc, lan, arg, tru]]
        # fetch the cartesian
        pos, vel = self.cartesian_from_kepler(parent, sma, ecc, inc, lan, arg, tru)
        # compute some orbital properties that will be helpful
        self.period = D(2) * trig_funcs.pi() * (sma**3 / (D(const.G.value) * (self.mass + parent.mass))).sqrt()
        # add to the parent position and velocity
        pos += parent.pos
        vel += parent.vel
        # and store
        self.pos = pos
        self.vel = vel
        self._pos = np.copy(pos)
        self._vel = np.copy(vel)
        # return for checking purposes
        return pos, vel

    # convert from keplerian elements to cartesian
    def cartesian_from_kepler(self, parent, sma, ecc, inc, lan=0, arg=0, tru=0):
        ''' Take the orbital elements, and convert them to cartesian coordinates.
            parent: the body object that this orbits around.

            The following define the orbit around the parent:
            sma: Semi-major axis, ecc: Eccentricity, inc: Inclination,
            lan: Longitude of ascending node, arg: Argument of periapse
            tru: True anomaly'''
        # ensure the parent is not itself
        if parent == self:
            raise Exception("Can't assign orbit of body around itself!")
        # ensure units
        sma = D(self.asUnit(sma, u.m).value)
        ecc = D(self.asUnit(ecc, u.m / u.m).value)
        inc, lan, arg, tru = [D(self.asUnit(val, u.radian).value) for val in [inc, lan, arg, tru]]
        # set of sine and cosine for the angular elements
        sininc, cosinc = trig_funcs.sin(inc), trig_funcs.cos(inc)
        sintru, costru = trig_funcs.sin(tru), trig_funcs.cos(tru)
        sinlan, coslan = trig_funcs.sin(lan), trig_funcs.cos(lan)
        sinarg, cosarg = trig_funcs.sin(arg), trig_funcs.cos(arg)

        # fetch the gravitational parameter of the system
        mu = D(const.G.value) * (parent.mass + self.mass)
        # sqrt(1- e^2)
        oneminusecc = (1 - ecc**2).sqrt()
        # eccentric anomaly
        eanom = trig_funcs.arctan( (oneminusecc * sintru) / (ecc + costru) )
        # distance to parent
        rad = (sma * oneminusecc) / (1 + ecc * costru)

        # position vector in orbital plane
        pos = rad * np.array([costru, sintru, 0])
        # velocity vector in orbital plane
        vel = ((mu * sma).sqrt() / rad) * np.array([-trig_funcs.sin(eanom), oneminusecc * trig_funcs.cos(eanom), 0])

        # rotation matrix values
        rx1 = cosarg*coslan - sinarg*cosinc*sinlan
        rx2 = sinarg*coslan + cosarg*cosinc*sinlan
        ry1 = cosarg*sinlan + sinarg*cosinc*coslan
        ry2 = cosarg*cosinc*coslan - sinarg*sinlan
        rz1 = sinarg*sininc
        rz2 = cosarg*sininc

        # transform position to bodycentric coordinate
        pos = np.array([pos[0] * rx1 - pos[1] * rx2,
                        pos[0] * ry1 + pos[1] * ry2,
                        pos[0] * rz1 + pos[1] * rz2])
        # transform velocity to bodycentric too
        vel = np.array([vel[0] * rx1 - vel[1] * rx2,
                        vel[0] * ry1 + vel[1] * ry2,
                        vel[0] * rz1 + vel[1] * rz2])

        # return the values
        return pos, vel
