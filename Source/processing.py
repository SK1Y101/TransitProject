# python modules
import matplotlib.pylab as plt
import astropy.units as u
import scipy.optimize
import pandas as pd
import numpy as np
import argparse

# my modules
import TransitProject as tp
import TransitProject.Simulation as ts
import TransitProject.modelling as tm

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
    parser.add_argument("--dataSource", help="The list of Sources to use data from.", nargs="+")
    parser.add_argument("--addSim", help="Compute the TTV for this system and overlay it on the O-C graph.", action="store_true")
    args = parser.parse_args()
    # return the argument input
    return args

def trimToSources(df, sources):
    ''' Remove any sources for data not wanted.
        df: The dataframe of midtransit points.
        sources: The list of allowed sources. '''
    # fetch the sources list
    dfsources = sorted(df["source"].unique())
    # if the provided sources were empty, use all of them
    if not sources:
        sources = dfsources
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

def fetchMidTransitTimes(exoplanet, dataSource):
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
    # trim to the desired sources
    df = trimToSources(df, dataSource)
    # and return
    return df

def simulateMidTransitTimes(exoplanet):
    import simulationPipeline as simP
    # define an arbitary start time
    start_time, end_time = 0, 60000
    # simulate and fetch the TTV
    TTV, TTVl, TTVu, TT = simP.runMinimalSim(args.planet, start_time, end_time/365.25, False)
    # convert parameters and return
    return TT*365.25*u.d, TTV, np.average((TTVl, TTVu))

def dfToParams(df):
    # fetch the star paramsfetchPara
    star = ts.bodyParams(df)[0]
    # and the params for each planet > (a, P, mu, ecc, inc, arg, t0)
    planets = ts.planetParams(df)
    # vstack the bodies
    bodies = np.vstack([star, planets])
    # remove orbital period as we're calculating it in models
    bodies = np.delete(bodies, 1, axis=1)
    # convert reduced mass
    # used to be m / mTot, change to raw mass in Msun
    bodies[:, 1] = df["mass"].to_numpy()
    # convert from m to AU, and days to years (This ensure the parameter space search is smaller)
    bodies[:, 0] = bodies[:, 0] * (1*u.m).to(u.AU).value
    bodies[:, 5] = bodies[:, 5] * (1*u.d).to(u.yr).value
    # return > (a, mu, ecc, inc, arg, t0)
    return bodies

# execte the program if called
if __name__ == "__main__":
    # fetch the program arguments
    args = fetchArgs()

    # fetch the system data
    systemdf = ts.fetchParams(args.planet)[0]

    # fetch the parameters for each body in the standard format
    bodies = dfToParams(systemdf)

    # if we passed a real planet:
    if ".csv" not in args.planet:
        # fetch the planetary data
        transitdf = fetchMidTransitTimes(args.planet, args.dataSource)
        # convert to useful units
        x = tp.DatetoHJD(transitdf["date"]).to_numpy() * u.d
        y = (transitdf["oc"].to_numpy() * u.min).to(u.s).value
        yerr = (0.5*(transitdf["ocel"]+transitdf["oceu"]).to_numpy() * u.min).to(u.s).value
    else:
        # compute TTV from simulation instead
        x, y, yerr = simulateMidTransitTimes(args.planet)
    # fetch the orbital periods of the bodies
    P = tm._extraModelParam_(tm.pSpaceToReal(bodies))[1].to(u.d)

    # define the models used. From the paper, they are:
    '''
        Model 1: Eq.13, 'n' interior, coplanar, non-interacting, zero-eccentricity, perturbing planets
        Model 2: Eq.16, 'n' interior, coplanar, weakly-interacting, non-zero-eccentricity, perturbing planets
        Model 3: Eq.19, 'n' interior, coplanar, weakly-interacting, non-zero-eccentricity, perturbing planets, including the transiting planets eccentricity
        Model 4: Eq.47, 'n' exterior, coplanar, non-interacting, non-zero-eccentricity, perturbing planets, with a zero eccentricity transiting planet
    '''
    models = [tm.model1, tm.model2, tm.model3]

    # plot the initial states of the models used
    xlim, xlimz = [0.5, 15.5], [5.5, 10.5]
    #tm.plotModels(x, y, yerr, P, models, tm.pSpaceToReal(bodies), xlim=xlim, xlimz=xlimz, fname="TTVTestModelComparison")

    # determine the optimal soltion of the fit of 'n' models with 'p' parameters and combined optimisation methods.
    solutions = tm.optimiser(x, y, yerr, models, bodies)
    # and graph them too
    tm.plotModels(x, y, yerr, P, solutions["models"], solutions["solutions"], xlim=xlim, xlimz=xlimz, fname="TTVTestModelFitting")

    # fetch the MCMC values for the models, as well as evaluating their information criterion
    MCMCVal = tm.determineUncertainties(x, y, yerr, solutions)

    # create a corner plot and generate a system chart for the best model
    tm.plotModelSystem(MCMCVal, bodies[:2], systemdf.iloc[:2]["name"])
