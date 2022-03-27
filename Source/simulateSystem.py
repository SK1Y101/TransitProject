import rebound
import numpy as np
from tqdm import tqdm
import matplotlib.pylab as plt

#from matplotlib.lines import Line2D

def simulateSystem(*objects):
    req_sims = 2**sum([sum(["_e" in x for x in list(object.keys())]) for object in objects])
    print("This requires {} simulations for all error values".format(req_sims))

    # create binary tree of all possible choices
    out=[]
    for object in objects:
        for y in object:
            if "_e" in y:
                thisval = object[y.replace("_e", "")], object[y]
                out.append([thisval[0]+thisval[1][0], thisval[0]+thisval[1][1]])

    def tree(lis):
        if lis:
            return [[lis[0][0]] + tree(lis[1:]), [lis[0][1]] + tree(lis[1:])]
        return []

    a = tree(out)

    for x in a:
        print(x)

'''simulateSystem({"name":"HAT-P-13",
                "mass":1.261, "mass_e":[0.029, 0.023]},
               {"name":"HAT-P-13b",
                "mass":0.883, "mass_e":[0.027, 0.047],
                "period":2.9162431, "period_e":0.0000006,
                "eccentricity":0.013, "eccentricity_e":[0, 0.013],
                "omega":197,
                "inclination":83.4, "inclination_e":0.26},
               {"name":"HAT-P-13c",
                "mass":14.61, "mass_e":[0.46, 0.48],
                "period":445.82, "period_e":0.11,
                "eccentricity":0.6554, "eccentricity_e":[0.0021, 0.0020],
                "omega":175.4})'''

simulateSystem({"name":"HAT-P-13",
                "mass":1.261, "mass_e":[0.029, 0.023]},
               {"name":"HAT-P-13b",
                "mass":0.883, "mass_e":[0.027, 0.047]})

sim = rebound.Simulation()
sim.add(m=1, P=None)

'''

sim = rebound.Simulation()
# HAT-P-13
sim.add(m=1.261)
# HAT-P-13 b
sim.add(m=1e-5 * 0.851, P=2.9162431+0.0000006, e=0.021, omega=181, inc=np.radians(90-83.4))
# HAT-P-13 b
sim.add(m=1e-5 * 14.28, P=445.82+0.11, e=0.6616, omega=176.7)
# move to the centre of mass
sim.move_to_com()

# step forward by one thousandth of an orbit
tstep = 2.91626 / 360

# ensure we collect 1000 transits
N=2000
transittimes = np.zeros(N)
p = sim.particles
i = 0

def gen():
    while i<N:
        yield

for _ in tqdm(gen()):
    y_old = p[1].y - p[0].y  # (Thanks to David Martin for pointing out a bug in this line!)
    t_old = sim.t
    sim.integrate(sim.t+tstep) # check for transits every 0.5 time units. Note that 0.5 is shorter than one orbit
    t_new = sim.t
    if y_old*(p[1].y-p[0].y)<0. and p[1].x-p[0].x>0.:   # sign changed (y_old*y<0), planet in front of star (x>0)
        while t_new-t_old>1e-7:   # bisect until prec of 1e-5 reached
            if y_old*(p[1].y-p[0].y)<0.:
                t_new = sim.t
            else:
                t_old = sim.t
            sim.integrate( (t_new+t_old)/2.)
        transittimes[i] = sim.t
        i += 1
        sim.integrate(sim.t+tstep)       # integrate 0.05 to be past the transit

#print(TTV)
print(transittimes)
_transits = np.copy(transittimes)


sim = rebound.Simulation()
# HAT-P-13
sim.add(m=1.261)
# HAT-P-13 b
sim.add(m=1e-5 * 0.851, P=2.9162431-0.0000006, e=0.021, omega=181, inc=np.radians(90-83.4))
# HAT-P-13 b
sim.add(m=1e-5 * 14.28, P=445.82-0.11, e=0.6616, omega=176.7)
# move to the centre of mass
sim.move_to_com()

# step forward by one thousandth of an orbit
tstep = 2.91626 / 360

# ensure we collect 1000 transits
transittimes = np.zeros(N)
p = sim.particles
i = 0

def gen():
    while i<N:
        yield

for _ in tqdm(gen()):
    y_old = p[1].y - p[0].y  # (Thanks to David Martin for pointing out a bug in this line!)
    t_old = sim.t
    sim.integrate(sim.t+tstep) # check for transits every 0.5 time units. Note that 0.5 is shorter than one orbit
    t_new = sim.t
    if y_old*(p[1].y-p[0].y)<0. and p[1].x-p[0].x>0.:   # sign changed (y_old*y<0), planet in front of star (x>0)
        while t_new-t_old>1e-7:   # bisect until prec of 1e-5 reached
            if y_old*(p[1].y-p[0].y)<0.:
                t_new = sim.t
            else:
                t_old = sim.t
            sim.integrate( (t_new+t_old)/2.)
        transittimes[i] = sim.t
        i += 1
        sim.integrate(sim.t+tstep)       # integrate 0.05 to be past the transit

print(transittimes)
__transits = np.copy(transittimes)

sim = rebound.Simulation()
# HAT-P-13
sim.add(m=1.261)
# HAT-P-13 b
sim.add(m=1e-5 * 0.851, P=2.9162431, e=0.021, omega=181, inc=np.radians(90-83.4))
# HAT-P-13 b
sim.add(m=1e-5 * 14.28, P=445.82, e=0.6616, omega=176.7)
# move to the centre of mass
sim.move_to_com()

# step forward by one thousandth of an orbit
tstep = 2.91626 / 360

# ensure we collect 1000 transits
transittimes = np.zeros(N)
p = sim.particles
i = 0

def gen():
    while i<N:
        yield

for _ in tqdm(gen()):
    y_old = p[1].y - p[0].y  # (Thanks to David Martin for pointing out a bug in this line!)
    t_old = sim.t
    sim.integrate(sim.t+tstep) # check for transits every 0.5 time units. Note that 0.5 is shorter than one orbit
    t_new = sim.t
    if y_old*(p[1].y-p[0].y)<0. and p[1].x-p[0].x>0.:   # sign changed (y_old*y<0), planet in front of star (x>0)
        while t_new-t_old>1e-7:   # bisect until prec of 1e-5 reached
            if y_old*(p[1].y-p[0].y)<0.:
                t_new = sim.t
            else:
                t_old = sim.t
            sim.integrate( (t_new+t_old)/2.)
        transittimes[i] = sim.t
        i += 1
        sim.integrate(sim.t+tstep)       # integrate 0.05 to be past the transit

slices=0.5

# recover transit timing variation from linear regression
A = np.vstack([np.ones(N), range(N)]).T
c, m = np.linalg.lstsq(A, transittimes, rcond=-1)[0]

TTV = (transittimes-m*np.array(range(N))-c)*(24.*365./2./np.pi)

print(TTV)
upper = _transits-transittimes
lower = __transits-transittimes
#plt.scatter(__transits)
plt.fill_between(np.arange(N), TTV+upper, TTV+lower,facecolor='gray',alpha=0.5)
plt.show()
'''
