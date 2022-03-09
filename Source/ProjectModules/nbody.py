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

# graphing modules
import matplotlib.animation as animation
import matplotlib.pylab as plt

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

# system class
class System(__Helper__):
    # initialisation
    def __init__(self):
        self.bodies = []
        self.propogator = self.euler
        self.grav_Func = self.do_gravity
        self.positions = {"barycentre":[]}
        self.totalmass = 0
        self.barycentre = []
        self.largest_apo = 0

    # start an updating plot of the system
    def plot(self, updatefunc=None, frames=None, interval=10, keepBarycentreConstant=True):
        # if we weren't given a function, use the inbuild
        updatefunc = updatefunc if updatefunc else self.plotpos
        # start a nice bar so we can see progress
        with tqdm() as pbar:
            # create the figure
            fig=plt.figure()
            # and axes
            ax = fig.add_subplot(projection="3d")
            # animation function
            def animate(i):
                # clear the graph
                ax.clear()
                # run our function to update the graph
                updatefunc(ax, i, keepBarycentreConstant)
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
        # limit the axes
        '''ax.set_xlim(-self.largest_apo, self.largest_apo)
        ax.set_ylim(-self.largest_apo, self.largest_apo)
        ax.set_zlim(-self.largest_apo, self.largest_apo)'''

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
                try:
                    apo = (1 + body.ecc) * body.sma
                    self.largest_apo = max(self.largest_apo, apo.to(u.m).value)
                except:
                    pass
            except:
                raise Exception("Object given is not a body: {}.".format(body))

    # compute gravitational forces between bodies
    def do_gravity(self, offset=1*u.m):
        # ensure our offset is valid
        offset = self.asUnit(offset, u.m)
        a_ = []
        # iterate on all the bodies
        for body1 in self.bodies:
            # reset acceleration from previous step
            a = np.zeros(3) * u.m / u.s**2
            # itterate on all the bodies again
            for body2 in self.bodies:
                # distance between
                dist = body1.pos - body2.pos
                dist_mag = np.linalg.norm(dist)
                # compute total force, using the offset to avoid infinities
                F = -const.G * body2.mass * body1.mass * dist / (dist_mag + offset)**3
                #compute force in each direction.
                a += F / body1.mass
            # append this bodies acceleration
            a_.append(a)
        # return the acceleration
        return self.toNpArray(a_)

    # compute the gravitational forces between bodies without using loops
    def do_gravity_vectorized(self, offset=1*u.m):
        # ensure our offset is valid
        offset = self.asUnit(offset, u.m)

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
        dist_mag[dist_mag == 0*u.m] = 1*u.m

        # compute the gravitational forces for each body in each direction
        F = const.G * mass_mat * dist / np.expand_dims(dist_mag, 2)**3

        # convert to 3d-acceleration
        a = F.sum(axis=1) / mass.reshape(-1, 1)

        return a

    # propogator
    def propogate(self, step=1, keepBarycentreConstant=True):
        # ensure the step is a number
        step = self.asUnit(step, u.s)

        # compute the gravitational acceleration between bodies
        acc = self.grav_Func(offset = 1*u.m)

        # set the barycentre position to nothing
        self.barycentre = self.toNpArray(self.asUnit([0,0,0], u.m))

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

    # propogator types

    # euler method, applies simple suvat to compute the bodies new position and velocity
    def euler(self, body, a, step):
        # y_{n+1} = y_n + h f(t_n, y_n)
        # ensure the step is a number
        step = self.asUnit(step, u.s)
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
        self.mass = self.asUnit(mass, u.kg)
        self.radius = self.asUnit(radius, u.m)
        # cartesian elements
        self.pos = np.zeros(3) * u.m
        self.vel = np.zeros(3) * (u.m / u.s)
        # keplerian elements
        self.sma = 0 * u.m
        self.ecc = 0 * (u.m / u.m)
        self.inc = 0 * u.degree
        self.lan = 0 * u.degree
        self.arg = 0 * u.degree
        self.tru = 0 * u.degree

    # set a body to be orbiting another
    def keplerian(self, parent, sma, ecc, inc, lan=0, arg=0, tru=0):
        ''' Take the orbital elements around the body, and store the position and velocity in the body
            parent: the body object that this orbits around.

            The following define the orbit around the parent:
            sma: Semi-major axis, ecc: Eccentricity, inc: Inclination,
            lan: Longitude of ascending node, arg: Argument of periapse
            tru: True anomaly'''
        # set the body elements
        self.sma, self.ecc, self.inc, self.lan, self.arg, self.tru = sma, ecc, inc, lan, arg, tru
        # fetch the cartesian
        pos, vel = self.cartesian_from_kepler(parent, sma, ecc, inc, lan, arg, tru)
        # compute some orbital properties that will be helpful
        self.period = 2*np.pi*np.sqrt(sma**3 / (const.G * (self.mass + parent.mass))).to(u.s)
        # add to the parent position and velocity
        pos += parent.pos
        vel += parent.vel
        # and store
        self.pos = pos
        self.vel = vel
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
        sma = self.asUnit(sma, u.m)
        ecc = self.asUnit(ecc, u.dimensionless_unscaled)
        inc, lan, arg, tru = self.asUnit([inc, lan, arg, tru], u.degree)
        # set of sine and cosine for the angular elements
        sininc, cosinc = np.sin(inc), np.cos(inc)
        sintru, costru = np.sin(tru), np.cos(tru)
        sinlan, coslan = np.sin(lan), np.cos(lan)
        sinarg, cosarg = np.sin(arg), np.cos(arg)

        # fetch the gravitational parameter of the system
        mu = const.G * (parent.mass + self.mass)
        # sqrt(1- e^2)
        oneminusecc = np.sqrt(1 - ecc**2)
        # eccentric anomaly
        eanom = np.arctan( (oneminusecc * sintru) / (ecc + costru) )
        # distance to parent
        rad = (sma * oneminusecc) / (1 + ecc * costru)

        # position vector in orbital plane
        pos = rad * np.array([costru, sintru, 0])
        # velocity vector in orbital plane
        vel = (np.sqrt(mu * sma) / rad) * np.array([-np.sin(eanom), oneminusecc * np.cos(eanom), 0])

        # rotation matrix values
        rx1 = cosarg*coslan - sinarg*cosinc*sinlan
        rx2 = sinarg*coslan + cosarg*cosinc*sinlan
        ry1 = cosarg*sinlan + sinarg*cosinc*coslan
        ry2 = cosarg*cosinc*coslan - sinarg*sinlan
        rz1 = sinarg*sininc
        rz2 = cosarg*sininc

        # transform position to bodycentric coordinate
        pos = self.toNpArray([pos[0] * rx1 - pos[1] * rx2,
                              pos[0] * ry1 + pos[1] * ry2,
                              pos[0] * rz1 + pos[1] * rz2])
        # transform velocity to bodycentric too
        vel = self.toNpArray([vel[0] * rx1 - vel[1] * rx2,
                              vel[0] * ry1 + vel[1] * ry2,
                              vel[0] * rz1 + vel[1] * rz2])

        # return the values
        return pos, vel
