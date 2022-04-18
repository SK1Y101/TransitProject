# python modules
import matplotlib.pylab as plt
import scipy.optimize
import pandas as pd
import numpy as np
import argparse

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

def fetchArgs():
    ''' fetch the planet from command line input. '''
    parser = argparse.ArgumentParser(description="Process transit timing variations.")
    # add any arguments
    parser.add_argument("--planet", help="The name of the planet to process TTV for.")
    args = parser.parse_args()
    # ensure we are not missing any required
    if not args.planet:
        raise Exception("No planet provided")
    # return the argument input
    return args

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

def dateToYearFloat(date):
    nanoSecToDay = 86400 * 1E9
    if not isinstance(date, str):
        return pd.to_datetime(date).astype(int) / nanoSecToDay
    return pd.to_datetime(date).value / nanoSecToDay

def fetchFit(x, y, err=None, f=lambda x,a,b:a*x+b):
    pars, corr = scipy.optimize.curve_fit(f, x, y, sigma=err, maxfev=5000)
    return lambda x:f(x, *pars)

def plotMidtransits(df):
    # fetch the sources list
    sources = sorted(df["source"].unique())
    # if there are any
    if sources:
        df = df.loc[df["source"]!="ETD"]
        df = df.sort_values(by="date", ascending=True)

        def f(x, a, b, c, d):
            return a*np.sin(x*(2*np.pi/b)+c) + d

        fig = plt.figure()
        ax = fig.add_gridspec(2, hspace=0).subplots(sharex=True)
        [ax1, ax2] = ax

        dateAsFloat = dateToYearFloat(df["date"])
        dateAsDate = pd.to_datetime(df["date"])
        dateminmax = [min(dateAsFloat), max(dateAsFloat)]
        daterange = np.linspace(*dateminmax, 10000)

        linefit = fetchFit(dateAsFloat, df["oc"], df["oce"], f=lambda x,a,b:a*x+b)
        df["oc"] -= linefit(dateAsFloat)
        fitfunc = fetchFit(dateAsFloat, df["oc"], df["oce"], f)
        ax1.plot(daterange, fitfunc(daterange), "b", label="Model fit")
        ax2.plot(dateminmax, [0, 0], "b", label="Model fit")
        # fetch the data
        for source in sources:
            data = df.loc[df["source"]==source]
            dateAsFloat = dateToYearFloat(data["date"])
            ax1.errorbar(pd.to_datetime(data["date"]), data["oc"], yerr=data["oce"], label=source, fmt="o")
            ax2.errorbar(pd.to_datetime(data["date"]), data["oc"]-fitfunc(dateAsFloat), yerr=data["oce"], label=source+" residuals", fmt="o")

        #dates = pd.to_datetime(df["date"])
        #plt.plot([min(dates), max(dates)], [0,0])
        for axs in ax:
            axs.set_title(df.index.name)
            #axs.gcf().autofmt_xdate()
            axs.set_ylabel("O-C (minutes)")
            axs.set_xlabel("Date of observation")
            axs.legend()
        plt.show()

# execte the program if called
if __name__ == "__main__":
    # fetch the program arguments
    args = fetchArgs()
    # fetch the planetary data
    df = fetchMidTransitTimes(args.planet)
    # plot the data
    plotMidtransits(df)
