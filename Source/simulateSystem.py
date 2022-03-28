import rebound
import numpy as np
from tqdm import tqdm
import matplotlib.pylab as plt
from time import sleep as wait

#from matplotlib.lines import Line2D

def simulateSystem(*objects):
    # add params to object
    params, out, noneVal = ["mass", "period", "eccentricity", "omega", "inclination"], [[], []], 0.123456789
    # construct the 2d array of values
    for obj in objects:
        for param in params:
            vs = [noneVal]*2
            if param in obj:
                if not isinstance(obj[param+"_e"], list):
                    obj[param+"_e"] = [obj[param+"_e"]]*2
                vs = [obj[param]+obj[param+"_e"][0], obj[param]-obj[param+"_e"][1]]
            out[0].append(vs[0])
            out[1].append(vs[1])

    # construct the entire possibility space
    out, maxperms, _out = np.array(out), len(out[0]), []
    for x in range(2**maxperms):
        # fetch a binary number
        thisbin = "{:0{}b}".format(x, maxperms)
        # convert to array of truthy/falsey
        _out.append([out[int(thisbin[x])][x] for x in range(maxperms)])
    _out = np.unique(_out, axis=0)

    print(len(_out))

    def val(x):
        if x!=noneVal:
            return x

    TTV, N = [], 1000

    # construct a simulation for each possibility
    for x in tqdm(_out):
        sim = rebound.Simulation()
        i = 0
        for obj in objects:
            try:
                inc=np.radians(90-val(x[i+4]))
            except:
                inc=None
            sim.add(m=val(x[i+0]), P=val(x[i+1]), e=val(x[i+2]), omega=val(x[i+3]), inc=inc)
            i+=5
        sim.move_to_com()
        # simulate till TTV
        TTV100 = simulateTTV(sim, N)
        plt.scatter(np.arange(N), TTV100)
    plt.show()

def simulateTTV(sim, N):
    tstep = 1 / 360
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

    A = np.vstack([np.ones(N), range(N)]).T
    c, m = np.linalg.lstsq(A, transittimes, rcond=-1)[0]

    TTV = (transittimes-m*np.array(range(N))-c)*(24.*365./2./np.pi)

    return TTV


simulateSystem({"name":"HAT-P-13",
                "mass":1.220, "mass_e":[0.059, 0.1]
                },
               {"name":"HAT-P-13b",
                "mass":0.883E-5, "mass_e":0,#[0.027E-5, 0.047E-5],
                "period":2.9162431, "period_e":0.0000006,
                "eccentricity":0.013, "eccentricity_e":0,#[0, 0.013],
                "omega":197,"omega_e":0,#[32, 37],
                "inclination":83.4, "inclination_e":0.26,
                },
               {"name":"HAT-P-13c",
                "mass":14.61E-5, "mass_e":0,#[0.46E-5, 0.48E-5],
                "period":445.82, "period_e":0.11,
                "eccentricity":0.6554, "eccentricity_e":0,#[0.0021, 0.0020],
                "omega":175.4, "omega_e":0,#0.22,
                })

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
sim.integrate(1)

fig, ax_main, ax_top, ax_right = rebound.OrbitPlot(sim, slices=0.5, xlim=[-2.,2], ylim=[-2.,2], color=True, unitlabel="[AU]")
fig.suptitle("HAT-P-13")
fig.tight_layout()
fig.savefig("test.png")
plt.show()
fig.show()'''
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
