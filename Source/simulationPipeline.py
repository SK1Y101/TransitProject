# python modules
from multiprocessing import current_process as curProc
import matplotlib.pylab as plt
from tqdm import tqdm, trange
import numpy as np
import rebound

# my modules
import TransitProject.Simulation as ts

# define the simulation parameters
''' choice of mass, period or semimajor axis, eccentiricty, inclination, argument of periapsis '''
''' if both a period and semimajor axis is given, the script will default to SMA '''
params = ["mass", "sma", "ecc", "inc", "arg"]

# fetch the parameters for the system
df = ts.fetchParams("HAT-P-13 b")

# compute the number of transits needed
times = 2008, 2022
N=int(np.ceil((times[1]-times[0])*365.25/df.iloc[1]["per"]))

# construct the aray of simulation parameters
simArray = ts.constructSimArray(df, params)

# fetch the transit times from simulation
TT = ts.fetchTT(simArray, params, N)

# linear regression to fetch transit times
A = np.vstack([np.ones(N), range(N)]).T
c, m = np.linalg.lstsq(A, TT[0], rcond=-1)[0]

# compute the TTV
TTV = TT[0] - m*np.array(range(N)) - c

# error area
#plt.fill_between(range(N), TT[0]-TT[1], TT[0]+TT[2], color="gray")
#plt.fill_between(range(N), TTV-TT[1]*1E-8, TTV+TT[2]*1E-8, color="gray", label="TTV error")
# main value
#plt.plot(range(N), TT[0], label="From archive", color="black")
plt.plot(TTV, color="black", label="Predicted TTV")
#labels
plt.xlabel("Epoch")
plt.ylabel("Time [Years]")
plt.legend()
plt.show()


'''
# Benchmarking suite
# Timing how long things take
params = ["mass", "sma"]#, "ecc", "inc", "arg"]
df = ts.fetchParams("HAT-P-13 b")

import time
times = []
Ns = [5, 10, 25]#, 50, 100]
workers = np.arange(4, 16)
for N in Ns:
    subtimes = []
    for worker in workers:
        start = time.time()
        simArray = ts.constructSimArray(df, params, tqdmLeave=False)
        TT = ts.fetchTT(simArray, params, N, workers=worker, tqdmLeave=False)
        subtimes.append(time.time()-start)
    times.append(subtimes)


#print(times)
for x in range(len(Ns)):
    plt.scatter(workers, times[x], label="{} transits".format(Ns[x]))
    if len(workers)>1:
        plt.plot(workers, times[x][0] / (workers/workers[0]))#, label="Idealised runtime")
plt.yscale("log")
plt.legend()
plt.show()
'''
