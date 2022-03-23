import rebound
import numpy as np
from tqdm import tqdm
import matplotlib.pylab as plt
from scipy.fft import rfft, rfftfreq, irfft
from scipy.signal import find_peaks
import scipy as sp

from matplotlib.lines import Line2D

sim = rebound.Simulation()
# HAT-P-13
sim.add(m=1.261)
# HAT-P-13 b
sim.add(m=1e-5 * 0.851, P=2.91626, e=0.021, omega=181, inc=np.radians(90-83.4))
# HAT-P-13 b
sim.add(m=1e-5 * 14.28, P=445.81, e=0.6616, omega=176.7)
# move to the centre of mass
sim.move_to_com()

# step forward by one thousandth of an orbit
tstep = 2.91626 / 360

# ensure we collect 1000 transits
N=520
print(N)
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

# plot the position of bodies
fig, ax_main, ax_top, ax_right = rebound.OrbitPlot(sim, slices=slices, xlim=[-2.,2], ylim=[-2.,2], color=True, unitlabel="[AU]")
fig.suptitle("HAT-P-13")

legend_elements = [Line2D([0], [0], color=(1.,0.,0.), lw=1.5, label='HAT-P-13b Orbit'),
                   Line2D([0], [0], color=(0.,0.75,0.75), lw=1.5, label='HAT-P-13c Orbit'),]

ax_main.legend(handles = legend_elements, loc=2)

fig.tight_layout()
fig.savefig("HAT-P-13 layout.png", transparent=True, bbox_inches='tight')
fig.show()

print(sim.cite())

# recover transit timing variation from linear regression
A = np.vstack([np.ones(N), range(N)]).T
c, m = np.linalg.lstsq(A, transittimes, rcond=-1)[0]

TTV = (transittimes-m*np.array(range(N))-c)*(24.*365./2./np.pi)

#fetch the limit from the maximum TTV
lim = max(abs(TTV))

# plot the variation
fig = plt.figure(figsize=(10,5))
ax = plt.subplot()
plt.scatter(range(N), TTV, marker="x", c="red", s=15, label="Observed variation from Least Squares Regression")

ax.set_xlim([-10,N-10])
ax.set_ylim([-lim*1.1, lim*1.1])
ax.set_xlabel("Transit number")
ax.set_ylabel("TTV [hours]")
ax.tick_params(axis="x",direction="in", pad=-15)
ax.xaxis.labelpad = -20

ax2 = ax.twiny()
ax2.set_xlabel("Time [days]")
ax2.xaxis.set_label_position("bottom")
ax2.xaxis.tick_bottom()
ax2.xaxis.labelpad = 0
ax2.set_xlim([-10 * 2.91626, sim.t - 10*2.91626])

ax.legend(loc="best")

fig.suptitle("Simulated Transit Timing Variation in the HAT-P-13 System", fontsize="x-large")
fig.tight_layout()
plt.legend()
plt.savefig("Simulated TTV", transparent=True, bbox_inches='tight')
plt.show()

# fourier analysis
duration = transittimes[-1] - transittimes[0]
sampling_rate = N / duration

from scipy.fft import rfft, rfftfreq, irfft
# compute the fourier transform of the residuals
yf = rfft(TTV)
xf = rfftfreq(N, 1 / sampling_rate)
peaks = find_peaks(np.abs(yf), distance=10)[0]

print(peaks)
fig = plt.figure(figsize=(10,5))
ax = plt.subplot()
ax.plot(xf, np.abs(yf))
plt.show()


fig = plt.figure(figsize=(10,5))
ax = plt.subplot()
# remove the largest frequency
plt.scatter(range(N), TTV, marker="x")
plt.plot(irfft(yf))
ax.set_xlim([0,N])
ax.set_xlabel("Transit number")
ax.set_ylabel("TTV [hours]")
ax.set_ylim([-lim, lim])

plt.show()
