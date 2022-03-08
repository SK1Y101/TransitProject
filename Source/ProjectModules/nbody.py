# solve the equations of motion numerically
# propogate the position of n-bodies in a system
# expose the ability to change the numerical method

# each body is derived from a single class
# the body class needs to store: mass, position, velocity
# it would be nice if it could also store: radius, luminosity, parent

# import python modules
from astropy import constants as const
from astropy import units as u
import numpy as np

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
        self.positions = {}

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
                self.positions[body] = []
            except:
                raise Exception("Object given is not a body: {}.".format(body))

    # propogator
    def propogate(self, step=1):
        # iterate on all bodies in the system
        # compute their acceleration due to every other body
        # move their position forward by one time step

        # ensure the step is a number
        step = self.asUnit(step, u.s)

        # iterate on all the bodies
        for body1 in self.bodies:
            # reset acceleration from previous step
            body1.a = np.zeros(3) * u.m / u.s**2
            # itterate on all the bodies again
            for body2 in self.bodies:
                # skip self interactions
                if body1 == body2:
                    continue
                # distance between
                dist = body1.pos - body2.pos
                dist_mag = np.linalg.norm(dist)
                # compute total force
                a = -const.G * body2.mass / dist_mag**2
                #compute force in each direction.
                body1.a += a * dist / dist_mag
        # update position using euler method
        for body in self.bodies:
            # velocity
            body.vel += body.a * step
            # position
            body.pos += body.vel * step
            # update historical positions
            self.positions[body] += [np.copy(body.pos)]

    # propogator types
    def euler(self):
        pass

# body class
class Body(__Helper__):
    # initialisation
    def __init__(self, mass, position=[0,0,0], velocity=[0,0,0]):
        ''' Initialise the body class.
            mass: the mass of this body. This value will be converted to kilograms
            position: the starting 3-position for this body.
            velocity: the starting 3-velocity for this body.'''
        self.isABody = True
        self.mass = self.asUnit(mass, u.kg)
        self.pos = self.toNpArray(self.asUnit(position, u.m))
        self.vel = self.toNpArray(self.asUnit(velocity, u.m/u.s))

    # set a body to be orbiting another
    def orbit(self, parent, sma, ecc, inc, lan=0, arg=0, tru=0):
        ''' Take the orbital elements around the body, and store the position and velocity in the body
            parent: the body object that this orbits around.

            The following define the orbit around the parent:
            sma: Semi-major axis, ecc: Eccentricity, inc: Inclination,
            lan: Longitude of ascending node, arg: Argument of periapse
            tru: True anomaly'''
        # fetch the cartesian
        pos, vel = self.cartesian_from_kepler(parent, sma, ecc, inc, lan, arg, tru)
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
