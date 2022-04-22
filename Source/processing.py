# python modules
import matplotlib.pylab as plt
import scipy.optimize
import pandas as pd
import numpy as np
import argparse

# my modules
import TransitProject as tp
import TransitProject.Simulation as ts

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
    parser.add_argument("--dataSource", help="The list of Sources to use data from.", nargs="+", default=["exoclock"])
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
    # if the target does not contain a distinction between name and number, attempt to do so
    if "-" not in exoclockTarget:
        numidx = exoclockTarget.index(tp.findFloats(exoclockTarget)[0])
        exoclockTarget = exoclockTarget[:numidx]+"-"+exoclockTarget[numidx:]
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
    pars, cov = scipy.optimize.curve_fit(f, x, y, sigma=err, maxfev=1000000, bounds=bounds)
    # return the function that fits to this curve
    return lambda x:f(x, *pars), pars, cov

def checkEphemerides(transitdf, systemdf):
    ''' fit a line to the midtransit times to determine if the ephemerides agree with literature.
        transitdf: the dataframe of midtransit mid-transit times.
        systemdf: the dataframe of the system information.'''
    # fetch the date as a float (days since 1970-01-01)
    dateAsFloat = dateToYearFloat(transitdf["date"])
    # planet orbit in days
    period = systemdf.iloc[1]["per"] * 365.25
    print("Planet orbital period is:         {}.".format(tp.totimestring(period)))
    # fit a line to the o-c data
    linefit, lineparam, linecov = fetchFit(dateAsFloat, transitdf["oc"], transitdf["oce"])
    # offset per orbit in seconds
    offset = lineparam[0] * period
    print("Linear fit suggests it should be: {}.".format(tp.totimestring(period+offset)))
    # update the midtransit times with the linefit
    transitdf["oc"] -= linefit(dateAsFloat)

def plotMidtransits(transitdf, fitfunc):
    ''' plot the midtransit times to a nice chart.
        df: The dataframe of midtransit times.'''
    # fetch the sources list
    sources = sorted(transitdf["source"].unique())
    # ensure the transits are sorted by date
    transitdf = transitdf.sort_values(by="date", ascending=True)

    # fetch the x axis details
    dateAsFloat = dateToYearFloat(transitdf["date"])
    dateAsDate = pd.to_datetime(transitdf["date"])
    dateminmax = [min(dateAsFloat), max(dateAsFloat)]
    daterange = np.linspace(*dateminmax, 10000)

    # compute the fit
    fitfunc, _, _ = fetchFit(dateAsFloat, transitdf["oc"], transitdf["oce"], fitfunc,
                             # ensure the orbital period is reasonable
                             bounds=((-np.inf, 1, -np.inf, -np.inf), (np.inf, np.inf, np.inf, np.inf)))

    # plot laout
    fig = plt.figure()
    ax = fig.add_gridspec(2, hspace=0).subplots(sharex=True)
    [ax1, ax2] = ax
    # plot the grey central area
    ax1.plot(dateminmax, [0, 0], "lightgray")
    ax2.plot(dateminmax, [0, 0], "lightgray")
    # plot the model
    ax1.plot(daterange, fitfunc(daterange), "b", label="Model fit")
    secondfit = lambda x:6.7*np.sin(x*(2*np.pi/300) + np.pi)
    ax1.plot(daterange, secondfit(daterange), "r")
    # print the chi squared value
    print("Chi squared value of data:", sum([x**2 for x in transitdf["oc"]]))
    print("Chi squared value of computer fit:", sum([x**2 for x in transitdf["oc"]-fitfunc(dateAsFloat)]))
    print("Chi squared value of manual fit:", sum([x**2 for x in transitdf["oc"]-secondfit(dateAsFloat)]))
    # plot each datapoint
    for source in sources:
        data = transitdf.loc[transitdf["source"]==source]
        dateAsFloat = dateToYearFloat(data["date"])
        ax1.errorbar(pd.to_datetime(data["date"]), data["oc"], yerr=data["oce"], label=source, fmt="o")
        ax2.errorbar(pd.to_datetime(data["date"]), data["oc"]-fitfunc(dateAsFloat), yerr=data["oce"], label=source, fmt="o")
        ax2.errorbar(pd.to_datetime(data["date"]), data["oc"]-secondfit(dateAsFloat), yerr=data["oce"], label="myfit", fmt="o")

    #dates = pd.to_datetime(df["date"])
    #plt.plot([min(dates), max(dates)], [0,0])
    for axs in ax:
        axs.set_title(transitdf.index.name)
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
    transitdf = fetchMidTransitTimes(args.planet)
    # include the datasources we are using
    transitdf = trimToSources(transitdf, args.dataSource)
    # fetch the system data
    systemdf = ts.fetchParams(args.planet)[0]
    print(systemdf)
    # check the ephemerides
    checkEphemerides(transitdf, systemdf)
    # fitting function
    def fitfunc(x, a, b, c, d):
        return a*np.sin(x*(2*np.pi/b)+c) + d
    # plot the data with the fitting residuals
    plotMidtransits(transitdf, fitfunc)
