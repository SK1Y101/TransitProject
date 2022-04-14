# python modules
from multiprocessing import current_process as curProc
import rebound, argparse, scipy.optimize
import matplotlib.pylab as plt
from tqdm import tqdm, trange
import numpy as np

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
    # number of transits
    N = len(transitTimes)
    # linear function
    f = lambda x,a,b: a*x+b
    # decompose the inputs to x,y
    x, y = np.arange(N), transitTimes
    # fit to transit times
    pars, corr = scipy.optimize.curve_fit(f, x, y, maxfev=5000)
    # subtract from the transit times to get the variation
    TTV = transitTimes - f(x, *pars)
    # and return
    return TTV

# fetch the planet from command line input
parser = argparse.ArgumentParser(description="Simulate transit timing variations.")
parser.add_argument("--planet", help="The name of the planet to simulate TTV for.")
parser.add_argument("--years", help="The number of years to simulate TTV for.\nDefaults to 4.")
parser.add_argument("--workers", help="The number of workers to use for multiprocessing.\nDefaults to half the available CPU cores")
parser.add_argument("--precision", help="The precision (in seconds) of the simulation.\nDefaults to 1 second.")
parser.add_argument("--useerror", help="Include error bounds in the calculation.", action="store_true")
parser.add_argument("--forceUse", help="Force the simulation to be executed if no midtransit data is available", action="store_true")
args = parser.parse_args()
years, workers = float(args.years) if args.years else 4, int(args.workers) if args.workers else None
precision = float(args.precision if args.precision else 1) / 31557600
if not args.planet:
    raise Exception("No planet provided")

# define the simulation parameters
''' choice of mass, period or semimajor axis, eccentiricty, inclination, argument of periapsis '''
''' if both a period and semimajor axis is given, the script will default to SMA '''
params = ["mass", "sma", "ecc", "inc", "arg"]

# fetch the parameters for the system
df = ts.fetchParams(args.planet, forceUse=args.forceUse)

# compute the number of transits needed
N, a = int(np.ceil(years/df.iloc[1]["per"])), 1
# ensure the outermost planet transits at least ten times (to bring out the TTV signal)
N = max(N, int(np.ceil(N * 10 * max(df.per[df.per.notnull()]) / years)))

# construct the aray of simulation parameters
simArray = ts.constructSimArray(df, params, useerror=args.useerror)

# fetch the transit times from simulation
TT = ts.fetchTT(simArray, params, N, prec=precision, returnAll=True, workers=workers)

# compute the TTV for all values, convert to minutes
TTV = [computeTTV(tt) * 365.25 * 1440 for tt in TT]
# and the error bounds (upper and lower)
TTVa, TTVl, TTVu = tp.avMinMax(TTV, 0)

# error area
if args.useerror:
    plt.plot(TTVl, color="black", label="TTV Range (Lower estimate)")
    plt.plot(TTVu, color="black", label="TTV Range (Upper estimate)")
    plt.fill_between(range(N), TTVl, TTVu, color="gray", label="TTV error")
# main value
plt.plot(TTVa, color="black", label="Predicted TTV")
#labels
plt.xlim([0, int(np.ceil(years/df.iloc[1]["per"]))])
plt.xlabel("Epoch")
plt.ylabel("Time [Minutes]")
plt.legend()
plt.show()
