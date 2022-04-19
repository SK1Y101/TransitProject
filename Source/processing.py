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
    ''' fetch the arguments from command line input. '''
    parser = argparse.ArgumentParser(description="Process transit timing variations.")
    # add any arguments
    parser.add_argument("--planet", help="The name of the planet to process TTV for.", required=True)
    parser.add_argument("--dataSource", help="The list of Sources to use data from.", nargs="+", default="exoclock", required=True)
    args = parser.parse_args()
    # return the argument input
    return args

def trimToSources(df, sources):
    ''' Remove any sources for data not wanted.
        df: The dataframe of midtransit points.
        sources: The list of allowed sources. '''
    # fetch the sources list
    dfsources = sorted(df["source"].unique())
    # allowed sources
    allowed = set()
    # iterate on the list of desired
    for source in sources:
        # itterate on the list of posisble
        for dfsource in dfsources:
            # if we have any kind of match
            if source.capitalize() in dfsource.capitalize():
                # add to our allowed sourced
                allowed.add(dfsource)
    # locate all the source matches
    df = df.loc[df["source"].isin(list(allowed))]
    # if we have no data remaining, we should stop
    if len(df) == 0:
        raise Exception("Not enough data for {} with given sources.".format(df.index.name))
    # fetch the index name
    name = df.index.name
    # reset the indices
    df = df.reset_index(drop=True)
    # and reset the index name
    df.index.name = name
    # return the dataframe
    return df

def fetchMidTransitTimes(exoplanet):
    ''' fetch the midtransit times for a given exoplanet.
        exoplanet: the name of the exoplanet to search for the data of.'''
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
    ''' Convert a date string to a unix float (ie: 2008-06-31 ~2008.5)
        date: The date string to convert. '''
    # number of nanoseconds in a day
    nanoSecToDay = 86400 * 1E9
    # if the date isn't a string, this works
    if not isinstance(date, str):
        return pd.to_datetime(date).astype(int) / nanoSecToDay
    # otherwise do this
    return pd.to_datetime(date).value / nanoSecToDay

def fetchFit(x, y, err=None, f=lambda x,a,b:a*x+b, bounds=(-np.inf, np.inf)):
    ''' fit a set of datapoints (with or without errors) to a function.
        x, y: The individual datapoints.
        err: The error on the y data. (We assume x data is controlled)
        f: The function to fit. if not givem , will default to a line.'''
    # fetch the parameters and the correlation of the optimisation
    pars, cov = scipy.optimize.curve_fit(f, x, y, sigma=err, maxfev=5000, bounds=bounds)
    # return the function that fits to this curve
    return lambda x:f(x, *pars), pars, cov

def plotMidtransits(df):
    ''' plot the midtransit times to a nice chart.
        df: The dataframe of midtransit times.'''
    # fetch the sources list
    sources = sorted(df["source"].unique())
    # if there are any
    if sources:
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

        linefit, lineparam, linecov = fetchFit(dateAsFloat, df["oc"], df["oce"], f=lambda x,a,b:a*x+b)
        print(lineparam)
        df["oc"] -= linefit(dateAsFloat)
        fitfunc, _, _ = fetchFit(dateAsFloat, df["oc"], df["oce"], f, bounds=((0, 0, -10, -100), (1000, 500, 10, 100)))
        ax1.plot(dateminmax, [0, 0], "lightgray")
        ax2.plot(dateminmax, [0, 0], "lightgray")
        ax1.plot(daterange, fitfunc(daterange), "b", label="Model fit")
        # fetch the data
        for source in sources:
            data = df.loc[df["source"]==source]
            dateAsFloat = dateToYearFloat(data["date"])
            ax1.errorbar(pd.to_datetime(data["date"]), data["oc"], yerr=data["oce"], label=source, fmt="o")
            ax2.errorbar(pd.to_datetime(data["date"]), data["oc"]-fitfunc(dateAsFloat), yerr=data["oce"], label=source, fmt="o")

        #dates = pd.to_datetime(df["date"])
        #plt.plot([min(dates), max(dates)], [0,0])
        for axs in ax:
            axs.set_title(df.index.name)
            #axs.gcf().autofmt_xdate()
            axs.set_ylabel("O-C (minutes)")
            axs.set_xlabel("Date of observation")
            axs.legend()
        ax2.set_title("Residuals from fit")
        plt.show()

# execte the program if called
if __name__ == "__main__":
    # fetch the program arguments
    args = fetchArgs()
    # fetch the planetary data
    df = fetchMidTransitTimes(args.planet)
    # include the datasources we are using
    df = trimToSources(df, args.dataSource)
    # plot the data
    plotMidtransits(df)
