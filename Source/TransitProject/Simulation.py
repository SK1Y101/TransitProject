# python modules
from multiprocessing import current_process as curProc
from tqdm import tqdm, trange
import numpy as np
import rebound

# my modules
from . import *

# -----< Methods >-----

def mapToArchive(params):
    ''' Take parameter names as used in my code, and convert them to
        exoplanet archive names.
        params: The parameter, or list of parameters, to convert. '''
    # function to replace names
    def replace_error_name(archive, name):
        if "_e" in name:
            name, err = name.split("_e")
            return archive[name]+"err{}".format(err)
        return archive[name]
    # dictionary of mappings
    archive =  {"mass":"pl_bmassj",
                "per":"pl_orbper",
                "sma":"pl_orbsmax",
                "inc":"pl_orbincl",
                "ecc":"pl_orbeccen",
                "arg":"pl_orblper",}
    # map and return
    return [replace_error_name(archive, param) for param in params]

def _fetchSystem_(target, params=["mass", "sma", "per", "ecc", "inc", "arg"]):
    ''' Fetch the system information for a given target planet.
        target: The name of the planet to fetch the information for.
        params: List of parameters to fetch from the stored databases. '''
    # fetch the central body name
    star = target[:-1]
    # add errors to the paramaters, and map to archive names
    params = ",".join([",".join([param,param+"_e1",param+"_e2"]) for param in params]).split(",")
    aparams = mapToArchive(params)
    # load the dataframe
    archive = loadDataFrame("/raw_data/pscomppars.csv")
    # fetch all planets with the required hostname
    data = archive.loc[archive["hostname"].str.replace(" ", "") == star].reset_index()
    # fetch the target planet index
    targetIdx = np.where(data["pl_name"].str.replace(" ","") == target)[0][0]
    # ensure it is the first value in the dataframe
    if targetIdx != 0:
        # move the specified row to the top of the dataFrame
        data = moveRowToTop(data, targetIdx)
    # construct a new dataframe of only required parameters
    systemData = data.loc[:, ["pl_name"]+aparams]
    # convert mass to solar mass
    if "pl_bmassj" in systemData.columns:
        systemData.loc[:, ["pl_bmassj", "pl_bmassjerr1", "pl_bmassjerr2"]] *= 9.55E-4
    # convert mass to solar mass
    #convert(systemData, "pl_orbincl", lambda x: np.radians(90-x))
    # populate a new list with the stars data
    d1 = data.iloc[0]
    starIdx, starData = len(systemData.index), [d1["hostname"], d1["st_mass"], d1["st_masserr1"], d1["st_masserr2"]]
    # append it to the dataframe, filling any unspecified values with NaN
    systemData.loc[starIdx] = starData + [np.nan]*(len(systemData.iloc[0])-len(starData))
    # and shift it to the top
    systemData = moveRowToTop(systemData, starIdx)
    # rename all the columns to be more usefull
    systemData.columns = ["name"]+params
    # return the data
    return systemData

def _possibilitySpace_(poss):
    ''' Construct a 2D array of all possible values, if each value can have one of two states.
        poss: A (2,N) array of possible input values, where each column is the two values for a single state.'''
    # compute the maximum number of parameters
    totalParams = out.shape[1]
    # initial possibility space array
    pos = np.full((2**totalParams, totalParams), np.nan)
    # itterate over all possibility space
    for idx in trange(2**totalParams, desc="Filling possibility space"):
        # convert the current itteration number to a binary with as many positions as needed
        idx_b = "{:0{}b}".format(idx, totalParams)
        # use the binary representation of this iteration to select the zeroth or first row of the output
        pos[idx,:] = np.array([poss[int(idx_b[x])][x] for x in range(totalParams)])
    # return all possibilities
    return pos

# construct all simulation possibilities
def constructSimArray(df, params=["mass", "sma", "ecc", "inc", "arg"]):
    ''' Construct an array of all possible simulation parameters.
        Each row corresponds to a new simulation setup, where the columns are the parameters for each object.
        df: The dataframe of system information to fetch data from.
        params: The list of parameters to include in this simulation. '''
    # compute the total number of required parameters
    totalParams = len(params)*len(df.index)
    # create an output array for each objects error values
    out = np.full((2, totalParams), np.nan)
    # and an output array for the default parameters
    out_ = np.full((1, totalParams), np.nan)
    # index counter for the loop
    idx = 0

    # construct the array of error values
    for obj in range(len(df.index)):
        object = df.loc[obj]
        # for each parameter
        for param in params:
            # fetch the value and errors
            val = object[param]
            er1 = val + object[param+"_e1"]
            er2 = val + object[param+"_e2"]
            # add to the output
            out[:,idx] = np.hstack([er1, er2])
            out_[:,idx] = val
            # inrement the index counter
            idx+=1

    # fetch the entire possibility space of our combinations
    pos = _possibilitySpace_(out)
    # add the default simulation settings to the begining
    pos = np.vstack([out_, pos])
    # replace nan values with infinitiy so we know they are wrong
    pos = np.nan_to_num(pos, nan=np.inf)
    # remove nonunique values
    pos = np.unique(pos, axis=0)
    # remove completely empty rows
    pos = pos[~np.isinf(pos).all(axis=1)]
    # return the possibility space
    return pos

# fetch relevant data given exoplanet
def fetchParams(exoplanet, params=["mass", "sma", "per", "ecc", "inc", "arg"]):
    ''' Fetch all relevant data for an exoplanet. Will raise an exception if no transit data exists for the target.
        exoplanet: The name of the exoplanet to search for.
        params: The list of parameters to be used in this simulation. '''
    # convert to capitalised with no space for exopclock
    exoclockTarget = exoplanet.capitalize().replace(" ", "")
    # convert name to uppercase (leaving planet delimiter lowercase) and remove spaces for the archive
    if exoplanet[-1].isalpha():
        archiveTarget = exoplanet[:-1].upper().replace(" ", "")+exoplanet[-1].lower()
    else:
        # if there wasn't a planet delimiter, add one
        archiveTarget = exoplanet.upper().replace(" ", "")+"b"

    # check the planet has midtransit data
    if not inDataFrame("/raw_data/exoplanetList.csv", exoclockTarget):
        raise Exception("Mid-Transit data does not exist for target planet.")
    # fetch the archive data for the system
    archiveData = _fetchSystem_(archiveTarget, params)

    # return the archive dataFrame
    return archiveData

def valFromParam(value, array, param):
    ''' Fetch the value of a given parameter given an array to search over.
        If a parameter doesn't exist, or is infinity, this will return None.
        value: A parameter or list of parameters to search for.
        array: The array containing the information to search within.
        param: The list of parameters for the simulation. (This is used to index the value in the array).'''
    # if the value given is a list
    if isinstance(value, list):
        # reevaluate for each list value and return
        return [valFromParam(val, array, param) for val in value]
    # if the value is not in the parameters
    if value not in param:
        # return none
        return None
    # otherwise, return the part of the array that corresponds
    return noneIfInf(array[param.index(value)])

def _simulateTT_(sim, timestep, transits=1000, i=1, prec=1E-7):
    ''' Simulate a system and find the transit times.
        sim: A rebound Simulation object to integrate. (particle 1 is assumed to be the target)
        timestep: The timestep to integrate over.
        transits: The number of transits to find the times for.
        i: The vertical position of the TQDM bar.
        prec: The temporal precision of the transit. Defaults to 1 second.'''
    # transit time array
    tarray, n = np.zeros(transits), 0
    # reference to the particles in the simulation
    p = sim.particles
    # create a progress bar for this simulation
    pbar = tqdm(total=transits, desc="Simulating transits", leave=False, position=i)
    # keep iterating until we have the required number of transits
    while n < transits:
        # fetch the position of our target relative to the central body and the current time
        y_0, t_0 = p[1]-p[0], sim.t
        # step forward in time
        sim.integrate(sim.t+timestep)
        # fetch the new position and time
        y_1, t_1 = p[1]-p[0], sim.t
        # if the target has passed infront of the central body (ie: the sign has changed from positive to negative)
        # and the target is infron of the central body
        if y_0.y*y_1.y < 0 and y_1.x > 0:
            # iterate until we have the desired precision
            while t_1-t_0 > prec:
                # if the sign of the relative y position changed this timestep
                if y_0.y*(p[1]-p[0]).y < 0:
                    # set this time as our upper bound
                    t_1 = sim.t
                else:
                    # set this time as our lower bound
                    t_0 = sim.t
                # step forward by the difference in our times
                sim.integrate(0.5 * (t_0+t_1))
            # add the transit time to the array
            tarray[n] = sim.t
            # increment our transit counter
            n += 1
            # increment the progress bar
            pbar.update(1)
            # and step forward past the transit we just found
            sim.integrate(sim.t+timestep)
    # return the found transit times
    # because we are working with G=1, a=1 [AU], then t [year]=2pi
    return tarray / (2*np.pi)

def _arrayToSim_(thisSim, params):
    ''' Generates a list of simulations from an array of simulation parameters.
        thisSim: The array of simulation parameters.
        params: The list of parameter names. used to map the array to usable values. '''
    # fetch the number of objects
    objs = len(thisSim) / len(params)
    # create the simulation
    sim = rebound.Simulation()
    # iterate on each object
    for obj in np.split(thisSim, objs):
        # fetch this objects paramaters
        mass, sma, per, inc, ecc, arg = valFromParam(["mass", "sma", "per", "inc", "ecc", "arg"], obj, params)
        # if we have a semimajor axis
        if sma:
            sim.add(m=mass, a=sma, e=ecc, omega=arg, inc=inc)
        else:
            sim.add(m=mass, P=per, e=ecc, omega=arg, inc=inc)
    # move to the centre of mass
    sim.move_to_com()
    # and return the simulation
    return sim

def _simulateTransitTimes_(simArray, params, transits=1000):
    ''' Locate all transit times for a given simulation setup.
        simArray: The array of simulation parameters.
        params: The list of parameter names. used to map the array to usable values.
        transits: The number of transits to locate. '''
    # fetch this simulation
    sim = _arrayToSim_(simArray, params)
    # fetch the smallest orbital period in the system
    p_orb = min([x.P for x in sim.particles[1:]])
    # ensure our timestep is either hourly, or one tenth the smallest orbital period
    timestep = min(p_orb*0.1, 1 / (24*365))
    # fetch the position from the worker ID
    pos = curProc()._identity[0] if curProc().name != "MainProcess" else None
    # simulate for the chosen number of transits
    TT = _simulateTT_(sim, timestep, transits, pos)
    # and return the transit times
    return TT

def fetchTT(simArray, params, transits=1000):
    ''' Create a multiprocessing job to simulate all possible configurations of a planetary system.
        returns the transit times, the upper, and the lower error bounds.
        simArray: The 2D array of all possible simulation parameters.
        params: The list of parameter names. used to map the array to usable values.
        transits: The number of transits to simulate for. '''
    inputs = toItterableInput(simArray, params, transits, keep=(1,))
    # run the multiprocessing job
    TT = parallelJob(_simulateTransitTimes_, inputs, outType=np.array)
    # compute the mininimum time and maximum time for each transit
    av = TT[0] if TT.shape[0] > 1 else TT
    mn = TT.min(axis=0)
    mx = TT.max(axis=0)
    # return the transit time, and error values
    return av, av-mn, mx-av
