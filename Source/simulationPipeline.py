# python modules
from multiprocessing import current_process as curProc
from matplotlib.ticker import AutoMinorLocator
import rebound, argparse, scipy.optimize
import matplotlib.pylab as plt
from tqdm import tqdm, trange
import pandas as pd
import numpy as np

# my modules
import TransitProject.Simulation as ts
import TransitProject as tp

'''
Goals: Collect system information from a local database.
       Setup simulation from known system parameters.
       Compute timing transit variations from simulation.
'''

def fetchArgs():
    ''' fetch the program arguments from command line input. '''
    parser = argparse.ArgumentParser(description="Simulate transit timing variations.")
    # select planet
    parser.add_argument("--planet", help="The name of the planet to simulate TTV for.", required=True)
    # select timeframe
    parser.add_argument("--years", help="The number of years to simulate TTV for.\nDefaults to 4.", default=4)
    parser.add_argument("--simYears", help="The minimum number of years to simulate, so that TTV signals are more stable.")
    parser.add_argument("--precision", help="The precision (in seconds) of the simulation.\nDefaults to 0.1 second.", default=1E-6 / 31557600)
    parser.add_argument("--startTime", help="The date string that dictates the starting time.", default="2000-01-01")
    parser.add_argument("--endTime", help="The date string that dictates the ending time. Requires startTime to be set")
    # select plotting and simulation setups
    parser.add_argument("--useError", help="Include error bounds in the calculation.", action="store_true")
    parser.add_argument("--plotErrorArea", help="Requires error bounds in the calculation, will plot the error output.", action="store_true")
    parser.add_argument("--workers", help="The number of workers to use for multiprocessing.\nDefaults to half the available CPU cores", default=None)
    parser.add_argument("--forceUse", help="Force the simulation to be executed if no midtransit data is available", action="store_true")
    # limit the simulation
    parser.add_argument("--unperturbed", help="Run the simulation with only the target planet", action="store_true")
    parser.add_argument("--limitError", help="Run the simulation with only errors on the mass and period/semimajor axis", action="store_true")
    # save the simulation output
    parser.add_argument("--saveAs", help="Filename to save the simulation results as.")
    args = parser.parse_args()
    # if we have a start and end time, then overwrite the years
    if args.endTime:
        # ensure we end after we start
        if pd.to_datetime(args.startTime) >= pd.to_datetime(args.endTime):
            raise Exception("Simulation cannot end before it has begun! (No time travelling shenanigans allowed!)")
        # compute the new years value
        args.years = (pd.to_datetime(args.endTime)-pd.to_datetime(args.startTime)).total_seconds() / 31557600
    # the year quantity is finicky if i don't do this
    args.years = float(args.years) if args.years else 4
    args.simYears = float(max(args.years, args.simYears if args.simYears else 0))
    # set use error to true if limit error is enabled
    args.useError = args.limitError or args.useError
    # convert to float
    args.precision = float(args.precision)
    # return the argument input
    return args

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

def runTTVPipeline(df, params, args):
    ''' run the pipeline that takes a set of requirements and produces a s imulates TTV.
        df: The planetary data dataframe.
        params: The parameters to simulate.
        args: The program arguments. '''
    # number of desired transits
    Nd = int(np.ceil(args.years/df.iloc[1]["per"]))
    # number of minimum simulated transits
    Ns = int(np.ceil(args.simYears/df.iloc[1]["per"]))
    # number of transits to ensure the outermost planet orbits 10 times
    No = int(np.ceil(10 * max(df.per[df.per.notnull()]) / df.per.iloc[1]))
    # set number of transits to simulate to the maximum of the three
    N = max(Nd, Ns, No)

    # construct the aray of simulation parameters
    simArray = ts.constructSimArray(df, params, useerror=args.useError, limitError=args.limitError)
    #ts.plotSystem(ts._arrayToSim_(simArray[0], params), df)

    # fetch the transit times from simulationsimYears
    TT = ts.fetchTT(simArray, params, N, prec=args.precision, returnAll=True, workers=args.workers)

    # compute the TTV for all values, convert to seconds
    TTV = [computeTTV(tt) * 31557600 for tt in TT]

    # return the TTV and raw transit times
    return TTV, TT

def plotSimulation(TTVMinMaxAv, df, args):
    ''' plot the results of a TTV simulation. '''
    # decompose TTV min, max, and average
    TTVa, TTVl, TTVu = TTVMinMaxAv
    # compute the y axis scaling.
    rescale = ts.computeScale(TTVu)
    N = len(TTVa)
    # plot an error area
    if args.plotErrorArea:
        plt.plot(TTVl*rescale[1], color="black", label="TTV Range (Lower estimate)")
        plt.plot(TTVu*rescale[1], color="black", label="TTV Range (Upper estimate)")
        plt.fill_between(range(N), TTVl*rescale[1], TTVu*rescale[1], color="gray", label="TTV error")
    # main value
    if args.useError:
        err = np.vstack([(TTVa-TTVl), (TTVu - TTVa)]) * rescale[1]
        plt.errorbar(range(N), TTVa*rescale[1], yerr=err, color="black", label="Predicted TTV", fmt="o")
    else:
        plt.plot(range(N), TTVa*rescale[1], color="black", label="Predicted TTV")#, s=15)
    # exoplanet name
    eName = args.planet.split(",")[1] if ".csv" in args.planet else args.planet
    # title
    plt.title("Simulated TTV for {} over a period of {}.".format(eName, tp.totimestring(args.years*365.25)))

    # ---< Axis shenanigans >---
    # fetch axis objects
    ax = plt.gca()
    # y axis
    ax.set_ylabel("TTV [{}]".format(rescale[0].capitalize()))
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    # plot epoch number
    ax.set_xlabel("Epoch")
    ax.tick_params(axis="x",direction="in", pad=-15)
    ax.xaxis.labelpad = -20
    ax.set_xlim([-1, int(np.ceil(args.years/df.iloc[1]["per"]))])
    plt.legend()

    # pot as years if the numbers start getting unweildy
    asYears = args.years >= 20
    orbitscale = df.iloc[1]["per"]*(1 if asYears else 365.25)
    # plot day/year number
    ax2 = ax.twiny()
    ax2.set_xlabel("Time [{}]".format("Years" if asYears else "Days"))
    ax2.xaxis.set_label_position("bottom")
    ax2.xaxis.tick_bottom()
    ax2.xaxis.set_minor_locator(AutoMinorLocator())
    ax2.xaxis.labelpad = 0
    ax2.set_xlim([axval*orbitscale for axval in ax.get_xlim()])

    # show the final plot
    plt.show()

def saveSimulation(system, TTV, args):
    ''' save the system layout and computed ttv for a given simulation.
        system: the system data dataframe
        ttv: the list of computed ttvs.
        args: the program commandline arguments.'''
    # if we are svaing
    if args.saveAs:
        # create the simulation save location
        tp.makeFolder("/sim_data/")
        # save the system information
        tp.saveDataFrame(system.replace(np.nan, 0), "/sim_data/"+args.saveAs+"_system")
        # convert ttv times to dataframe
        TTVt = pd.DataFrame(np.array(TTV).T)
        # ensure the columns are strings, not numbers
        TTVt.columns = [str(x) for x in list(TTVt.columns)]
        # save the TTV times
        tp.saveDataFrame(TTVt, "/sim_data/"+args.saveAs+"_TTV")

def runMinimalSim(target, startTime="2000-01-01", years=4, unperturbed=False):
    ''' Run a simulation with minimal input parameters, so that other modules can fetch.
        target: The transiting planet to simulate.
        startTime: The date to start.
        years: The number of years to simulate for.
        unperturbed: Whether we should only simulate the target planet.
        extraPlanet: Include an extra planet in the simulation. '''
    # construct fake commanline inputs
    class args:
        pass
    args = args()
    args.years = args.simYears = years
    args.precision = float(1E-6 / 31557600)
    args.limitError = args.useError = False
    args.workers = None

    # fetch the parameters for the system
    df, params = ts.fetchParams(target, forceUse=False, startTime=startTime)

    # if we are only simulating the target planet
    if unperturbed:
        df = df.iloc[:2]

    # run the TTV pipeline
    TTV, TT = runTTVPipeline(df, params, args)

    # as we only used one sim setup, fetch the timing
    TTV, TTVl, TTVu = tp.avMinMax(TTV, 0)
    TT = TT[0]

    # fetch only the TTV within the observation window
    endWindow = np.where(TT > years)[0][0] if max(TT) > years else len(TT)
    TTV, TTVl, TTVu, TT = TTV[:endWindow], TTVl[:endWindow], TTVu[:endWindow], TT[:endWindow]

    # and return the TTV and times
    return TTV, TTVl, TTVu, TT

# define the simulation parameters
# choice of mass, period or semimajor axis, eccentiricty, inclination, argument of periapsis
# if both a period and semimajor axis is given, the script will default to SMA
params = ["mass", "sma", "ecc", "inc", "arg", "mean"]

# start the script
if __name__ == "__main__":
    # fetch the program arguments
    args = fetchArgs()

    # fetch the parameters for the system
    df, params = ts.fetchParams(args.planet, forceUse=args.forceUse, startTime=args.startTime)

    # predict TTV if there are additonal planets
    if len(df) > 2:
        # compute the predicted TTV
        predTTV = ts.predictTTVMagnitude(df)
        # convert to useful value
        print("The predicted TTV for {} is on the order of {}".format(df.iloc[1]["name"], tp.totimestring(predTTV[0])))
        print("The predicted TTV is cyclic on the order of {}".format(tp.totimestring(predTTV[1])))

        # compute the effect due to GR
        GRTTV = ts.predictGREffect(df, args.years)
        print("\nThe effect due to GR is on the order of {:0.2f} {}".format(*tp.tosi(GRTTV[0], "s")))
        print("The effect due to GR is cyclic on the order of {}".format(tp.totimestring(GRTTV[1])))

        # compute the effect due to barycentre shift
        BCTTV = ts.predictBCEffect(df)
        print("\nThe effect due to Barycentre shift is on the order of {:0.2f} {}".format(*tp.tosi(BCTTV[0], "s")))
        print("The effect due to Barycentre shift is cyclic on the order of {}".format(tp.totimestring(BCTTV[1])))

    # if the planet is to be simulated alone, do that here
    if args.unperturbed:
        df = df.iloc[:2]

    # run the pipeline
    TTV = runTTVPipeline(df, params, args)[0]

    # compute the minimum, maximum, and average
    TTVMinMaxAv = tp.avMinMax(TTV, 0)

    # if we want to save
    saveSimulation(df, TTV, args)

    # plot the output
    plotSimulation(TTVMinMaxAv, df, args)
