# python modules
import matplotlib.pylab as plt
import pandas as pd
import numpy as np

# my modules
import TransitProject as tp

'''
Goals: Create a simple extensible model that can take system parameters and output a TTV sinusoid.
       Compare model to simulated TTV to verify its accuracy (Minimise simulation time where possible).
       Provide model in a form that can be used for MCMC fitting from real data.

       Take observed midtransit times and fit various curves to the residuals.
       - ie: the model developed above, simple linear, polyfit, standard sine
       Determine whether the extensible model found is a better fit than others (ie: better than a line).

       If the model fit, find the additional parameters that best fit the data:
       - ie: fit for perturbing planetary period, for example.
       - can be found from fourier analysis, or from the MCMC model parameters probably.
'''

# model:
# required inputs:
# orbital elements of transiting planet
# physical parameters of transiting planet and parent star
# variable inputs:
# physical parameters and orbital elements of perturbing planet(s)

# output:
# sinusoidal curve with correct amplitude and frequency, as seen from simulation.

#-> eccentricity  == 0  -> sinusoid
#-> eccentricity --> 1  -> smooth sawtooth (or skewed sine)

#-> period --> 1/ frequency of the periodic function

#-> inclination -> sineness or tanness of the function

def fetchMidTransitTimes(exoplanet):
    # convert to capitalised with no space for exopclock
    exoclockTarget = exoplanet.capitalize().replace(" ", "")
    # fetch the datafil name for if this planet is available
    dataFile = "/raw_data/midTransitTimes/{}.csv".format(exoclockTarget)
    # load if available
    df = tp.loadIfAvailable(exoclockTarget, "/raw_data/exoplanetList.csv", dataFile)
    # assign the exoplanet to the dataframe index
    df.index.name = exoclockTarget
    # and return
    return df

def plotMidtransits(df):
    # fetch the sources list
    sources = sorted(df["source"].unique())
    # if there are any
    if sources:
        # fetch the data
        for source in sources:
            data = df.loc[df["source"]==source]
            data = data.sort_values(by="date", ascending=True)
            plt.errorbar(pd.to_datetime(data["date"]), data["oc"], yerr=data["oce"], label=source, fmt="o")

        #dates = pd.to_datetime(df["date"])
        #plt.plot([min(dates), max(dates)], [0,0])
        plt.title(df.index.name)
        plt.gcf().autofmt_xdate()
        plt.ylabel("O-C (minutes)")
        plt.xlabel("Date of observation")
        plt.legend()
        plt.show()

df = fetchMidTransitTimes("HAT-P-13b")
plotMidtransits(df)
