# python modules
from multiprocessing import current_process as curProc
from matplotlib.ticker import AutoMinorLocator
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

def fetchArgs():
    ''' fetch the program arguments from command line input. '''
    parser = argparse.ArgumentParser(description="Simulate transit timing variations.")
    # select planet
    parser.add_argument("--planet", help="The name of the planet to simulate TTV for.")
    # select timeframe
    parser.add_argument("--years", help="The number of years to simulate TTV for.\nDefaults to 4.")
    parser.add_argument("--precision", help="The precision (in seconds) of the simulation.\nDefaults to 0.1 second.")
    # select plotting and simulation setups
    parser.add_argument("--useError", help="Include error bounds in the calculation.", action="store_true")
    parser.add_argument("--plotErrorArea", help="Requires error bounds in the calculation, will plot the error output.", action="store_true")
    parser.add_argument("--workers", help="The number of workers to use for multiprocessing.\nDefaults to half the available CPU cores")
    parser.add_argument("--forceUse", help="Force the simulation to be executed if no midtransit data is available", action="store_true")
    # limit the simulation
    parser.add_argument("--unperturbed", help="Run the simulation with only the target planet", action="store_true")
    parser.add_argument("--limitError", help="Run the simulation with only errors on the mass and period/semimajor axis", action="store_true")
    args = parser.parse_args()
    # ensure we are not missing any required
    if not args.planet:
        raise Exception("No planet provided")
    # default arguments if not provided
    args.years = float(args.years) if args.years else 4
    args.workers = int(args.workers) if args.workers else None
    args.precision = float(args.precision if args.precision else 0.1) / 31557600
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
    # compute the number of transits needed
    N, a = int(np.ceil(args.years/df.iloc[1]["per"])), 1
    # ensure the outermost planet transits at least ten times (to bring out the TTV signal)
    N = max(N, int(np.ceil(N * 10 * max(df.per[df.per.notnull()]) / args.years)))

    # construct the aray of simulation parameters
    simArray = ts.constructSimArray(df, params, useerror=args.useError, limitError=args.limitError)

    # plot the current system layout
    #ts.plotSystem(ts._arrayToSim_(simArray[0], params), df, False)

    # fetch the transit times from simulation
    TT = ts.fetchTT(simArray, params, N, prec=args.precision, returnAll=True, workers=args.workers)

    # plot the system layout after the simulated time has passed
    #ts.plotSystem(ts._arrayToSim_(simArray[0], params), df, False)

    # compute the TTV for all values, convert to seconds
    TTV = [computeTTV(tt) * 31557600 for tt in TT]
    # and the error bounds (upper and lower)
    TTVMinMaxAv = tp.avMinMax(TTV, 0)

    # return the TTV
    return TTVMinMaxAv

def computeScale(val):
    ''' Compute the scale factor for a value. '''
    # maximum value
    mval = val if isinstance(val, (int, float)) else max(abs(val))
    # scale factors
    scales = {"milliseconds":0.001, "seconds":1, "minutes":60, "hours":3600, "days":86400}
    # standard scale output
    rescale = ("milliseconds", 0.001)
    # for each scale
    for scale in scales:
        # if the scale is smaller than the maximum value
        if scales[scale] < mval:
            # update the rescale factor
            rescale = (scale, 1/scales[scale])
    # return
    return rescale

def plotSimulation(TTVMinMaxAv, args):
    ''' plot the results of a TTV simulation. '''
    # decompose TTV min, max, and average
    TTVa, TTVl, TTVu = TTVMinMaxAv
    # compute the y axis scaling.
    rescale = computeScale(TTVa)
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
    plt.title("Simulated TTV for {} over a period of {} years.".format(eName, args.years))

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

# start the script
if __name__ == "__main__":
    # fetch the program arguments
    args = fetchArgs()

    # define the simulation parameters
    # choice of mass, period or semimajor axis, eccentiricty, inclination, argument of periapsis
    # if both a period and semimajor axis is given, the script will default to SMA
    params = ["mass", "sma", "ecc", "inc", "arg"]

    # fetch the parameters for the system
    df, params = ts.fetchParams(args.planet, forceUse=args.forceUse)
    if args.unperturbed:
        df = df.iloc[:2]

    # compute the predicted TTV
    predTTV = ts.predictTTVMagnitude(df)
    predScale = computeScale(predTTV)
    # convert to useful value
    print("The predicted TTV for {} is on the order of {:0.2f} {}".format(df.iloc[1]["name"], predTTV*predScale[1], predScale[0]))

    # compute the effect due to GR
    GRTTV = ts.predictGREffect(df, args.years)
    print("The predicted TTV due to GR is on the order of {:0.2f} {}".format(*tp.tosi(GRTTV[0], "s")))

    # run the pipeline
    TTVMinMaxAv = runTTVPipeline(df, params, args)

    # plot the output
    plotSimulation(TTVMinMaxAv, args)
