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
sdf = ts.fetchParams(args.planet)[0]

# assume the zero epoch is the first observation, and fetch period from the system details
t0, P = min(df["date"]), sdf.iloc[1]["per"] * 365.25
# compute the times difference between the zero epoch
dt = pd.to_datetime(df["date"]) - pd.to_datetime(t0)
# compute the transit number
transit = np.array([round(t.total_seconds()/86400) for t in (dt / P)])

# shorthands to information
dates = pd.to_datetime(df["date"])
oc  =
oce = [pd.to_timedelta(df["ocel"]), pd.to_timedelta(df["oceu"])]

'''
# plot the midtransit times
fig = plt.figure(figsize=(12,4))
plt.errorbar(transit, pd.to_datetime(df["date"]), fmt="x", zorder=10, label="observed transits")
plt.xlabel("Transit Number")
plt.ylabel("Date")
plt.title("Transits of {}".format(args.planet.capitalize()))
plt.legend()
plt.show()'''

# plot the midtransit times with linear ephemerides
fig = plt.figure(figsize=(12,4))
plt.errorbar(transit, dates, yerr=oce, fmt="o", zorder=10, label="observed transits")
plt.plot(transit, pd.to_datetime(t0)+pd.to_timedelta(transit*P, "d"), "--", zorder=2, label="linear fit")
plt.xlabel("Transit Number")
plt.ylabel("Date")
plt.title("Transits of {}".format(args.planet.capitalize()))
plt.legend()
plt.show()

# plot the residuals to the ephemeride fit
fig = plt.figure(figsize=(12,4))
