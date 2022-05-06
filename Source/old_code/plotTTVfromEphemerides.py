import TransitProject as tp
import TransitProject.Simulation as ts

import matplotlib.pylab as plt
import argparse
import numpy as np
import pandas as pd

def fetchArgs():
    ''' fetch the arguments from command line input. '''
    parser = argparse.ArgumentParser(description="Process transit timing variations.")
    # add any arguments
    parser.add_argument("--planet", help="The name of the planet to process TTV for.", required=True)
    args = parser.parse_args()
    # return the argument input
    return args

def fetchMidTransitTimes(exoplanet):
    ''' fetch the midtransit times for a given exoplanet.
        exoplanet: the name of the exoplanet to search for the data of.'''
    # convert to capitalised with no space for exopclock
    exoclockTarget = exoplanet.capitalize().replace(" ", "")
    # if the target does not contain a distinction between name and number, attempt to do so
    #if "-" not in exoclockTarget:
    #    numidx = exoclockTarget.index(tp.findFloats(exoclockTarget)[0])
    #    exoclockTarget = exoclockTarget[:numidx]+"-"+exoclockTarget[numidx:]
    # fetch the datafil name for if this planet is available
    dataFile = "/raw_data/midTransitTimes/{}.csv".format(exoclockTarget)
    # load if available
    df = tp.loadIfAvailable(exoclockTarget, "/raw_data/exoplanetList.csv", dataFile)
    # ensure we have enough datapoints
    if len(df) < 2:
        raise Exception("Cannot perform analysis with {} datapoint{}!".format(len(df), "s"*(len(df)!=1)))
    # warn if we have few datapoints
    if len(df) < 20:
        print("Warning! planet only has {} mid transit points. Acuracy may be diminished.".format(len(df)))
    # assign the exoplanet to the dataframe index
    df.index.name = exoclockTarget
    # and return
    return df

args = fetchArgs()

# fetch midtransit times
df = fetchMidTransitTimes(args.planet)

# fetch system details
sdf = ts.fetchParams(args.planet, params=["per", "t0"])[0]
target = sdf.iloc[1]

# shorthands to date information
dates = pd.to_datetime(df["date"])
oced = [pd.to_timedelta(df["ocel"]), pd.to_timedelta(df["oceu"])]
# shorthands to OC information
oce = [df["ocel"], df["oceu"]]
oc = df["oc"]

totalsec = lambda x: np.array([t.total_seconds() for t in x])
rounded = lambda x, a=0: np.array([round(t, a) for t in x])

# assume the zero epoch is the first observation, and fetch period from the system details
t0, P = target["t0"], target["per"] * 365.25
t0 = pd.to_datetime(t0 - 2400000, unit="d", origin=pd.Timestamp("1858-11-16 12:00"))
# compute the times difference between the zero epoch
dt = dates - t0
# compute the transit number
transit = rounded(totalsec(dt/P)/86400)

#xlimit of the graphs
limoffset = np.ceil((max(transit)-min(transit))*0.05)
xlim = [min(transit)-limoffset, max(transit)+limoffset]

# plot the midtransit times
fig = plt.figure(figsize=(12,8))
plt.errorbar(transit, dates, yerr=oced, fmt="x", zorder=10, label="observed transits")
plt.xlim(xlim)
plt.xlabel("Transit Number")
plt.ylabel("Date")
plt.title("Transits of {}".format(args.planet.capitalize()))
plt.legend()
plt.savefig("{}_Transits.png".format(args.planet.capitalize()), transparent=True)
plt.show()

# plot the midtransit times with linear ephemerides
fig = plt.figure(figsize=(12,8))
plt.errorbar(transit, dates, yerr=oced, fmt="x", zorder=10, label="observed transits")
plt.plot(transit, t0+pd.to_timedelta(transit*P, "d"), "--", zorder=2, label="linear fit", mfc="white")
middle = np.average(transit)
plt.text(middle, t0+pd.to_timedelta(np.median(middle)*P, "d"), \
         "\nT(n) = t0 + nP\n\nt0 = {}\nP = {:0.5f} d\n".format(t0, P), va="top")
plt.xlim(xlim)
plt.xlabel("Transit Number")
plt.ylabel("Date")
plt.title("Transits of {}".format(args.planet.capitalize()))
plt.legend()
plt.savefig("{}_Ephemerides.png".format(args.planet.capitalize()), transparent=True)
plt.show()

# plot the residuals to the ephemeride fit
fig = plt.figure(figsize=(12,8))
plt.errorbar(transit, oc, yerr=oce, fmt="o", zorder=10, mfc="white")
plt.plot(xlim, [0, 0], "--", zorder=2)
plt.xlim(xlim)
plt.xlabel("Transit Number")
plt.ylabel("O-C (minutes)")
plt.title("Transit Timing Variation of {}".format(args.planet.capitalize()))
plt.savefig("{}_TTV.png".format(args.planet.capitalize()), transparent=True)
plt.show()
