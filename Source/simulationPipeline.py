# python modules
from multiprocessing import current_process as curProc
import matplotlib.pylab as plt
from tqdm import tqdm, trange
import numpy as np
import rebound

# my modules
import TransitProject.Simulation as ts
import TransitProject as tp

'''
Goals: Collect system information from a local database.
       Setup simulation from known system parameters.
       Compute timing transit variations from simulation.
'''

def computeTTV(transitTimes):
    ''' Compute the transit timing variation from an array of transit times using linear regression.
        transitTimes: the array of transit times. '''
    N = len(transitTimes)
    # construct an empty 2d array for regression
    A = np.vstack([np.ones(N), range(N)]).T
    # perform regression to get coefficients of linear fit
    c, m = np.linalg.lstsq(A, transitTimes, rcond=-1)[0]
    # subtract from the transit times to get the variation
    TTV = transitTimes - m*np.arange(N) - c
    # and return
    return TTV

# define the simulation parameters
''' choice of mass, period or semimajor axis, eccentiricty, inclination, argument of periapsis '''
''' if both a period and semimajor axis is given, the script will default to SMA '''
params = ["mass", "sma", "ecc", "inc", "arg"]

# fetch the parameters for the system
df = ts.fetchParams("HAT-P-13 b")

# compute the number of transits needed
times = 2008, 2022
N=int(np.ceil((times[1]-times[0])*365.25/df.iloc[1]["per"]))

# manually removing errors
a, b, N = "ecc", .3, 1000
for x in params:
    for y in [1,2]:
        df[x+"_e{}".format(y)] = 0
df.loc[2, a] = df.loc[2, a+"_e1"] = df.loc[2, a+"_e2"] = b

# construct the aray of simulation parameters
simArray = ts.constructSimArray(df, params)

# fetch the transit times from simulation
TT = ts.fetchTT(simArray, params, N, returnAll=True)

# compute the TTV for all values
TTV = [computeTTV(tt) for tt in TT]
# and the error bounds (upper and lower)
TTVa, TTVl, TTVu = tp.avMinMax(TTV, 0)

# error area
#plt.fill_between(range(N), TTVl, TTVu, color="gray", label="TTV error")
# main value
for TTVb in TTV:
    plt.plot(TTVb[:300], color="black", label="Predicted TTV")
#labels
plt.xlabel("Epoch")
plt.ylabel("Time [Years]")
plt.legend()
plt.show()
