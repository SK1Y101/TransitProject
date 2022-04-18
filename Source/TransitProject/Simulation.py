# python modules
from multiprocessing import current_process as curProc
from tqdm import tqdm, trange
import rebound, reboundx
import numpy as np

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

def convert(df, col, convFunc):
    ''' Convert a dataframe column using a specified function.
        df: The dataframe to convert.
        col: The column(s) to convert.
        convFunc: The function that would convert from old to new (ie: np.radians)'''
    # if we had a list of columns
    if isinstance(col, list):
        for co in col:
            df = convert(df, co, convFunc)
    # oterwise
    else:
        # fetch the column values
        val = convFunc(df.loc[:,col])
        # and the insert position
        idx = list(df.columns).index(col)
        # Remove the old
        df = df.drop([col], axis=1)
        # and add the new
        df.insert(loc=idx, column=col, value=val)
    # return the new dataframe
    return df

def removeWhitespace(df):
    ''' Removes all whitespace from a dataframe. usefull for converting human readable into machine readable.
        df: The dataframe to convert. '''
    # remove the whitespace in the columns row
    df.columns = df.columns.str.strip()
    # for each column in the dataframe
    for col in df.columns:
        # remove whitespace
        df[col] = df[col].str.strip()
    # return the dataframe
    return df

def nestToLinear(lis):
    ''' convert a nested list to a linear one.
        lis: The list to convert. '''
    return list(np.array(lis).flatten())

def _fetchCustomSystem_(file, exoplanet):
    ''' Fetch the system information for a given planet in a custom system.
        file: The csv file containing the system information.
        exoplanet: The name of the target exoplanet.'''
    # define the allowed parameters
    params=["mass", "sma", "per", "ecc", "inc", "arg"]
    # get the dataset
    data = loadDataFrame(file)
    # remove any whitespace
    data = removeWhitespace(data)
    # define some conversion factors to change things into the simulation units
    conv = {"mass": {"":1, "S":1, "J": 1/1047.348644, "E":1/332946.0487, "M":1/27068702.8061, "Kg":1/1.98855E30},
            "sma":  {"":1, "AU":1, "Km":1/149597870.7, "m":1/149597870700},
            "per":  {"":1, "Y":1, "M": 1/12.008, "D":1/365.25, "h":1/8766, "m":1/52600, "s":1/31557600}}
    # fetch the columns that could contain a unit shorthand
    cols = nestToLinear([[cols, cols+"_e1", cols+"_e2"] for cols in list(conv.keys())])
    # for each column
    for col in cols:
        # for each value in the column
        for v, val in enumerate(data[col]):
            # fetch all the numbers in that row
            floats = findFloats(val)
            # if there are none, move on
            if not floats:
                continue
            # the starting value, and conversion key
            thisval, convKey = 0, col.split("_")[0]
            # work backwards through the floats
            for f in floats[::-1]:
                # fetch the position of the unit shorthand
                unitPos = val.index(f)+len(f)
                # fetch the unit and value
                value, unit = f, val[unitPos:]
                # add to thisvalue
                thisval += float(value)*conv[convKey][unit]
                # remove from the val
                val = val[:val.index(f)]
            # update the dataframe
            data.loc[v, [col]] = thisval
    # fetch all the columns that are meant to be numeric
    numCol = data.columns[data.columns!="name"]
    # convert to numeric
    data.loc[:, numCol] = data.loc[:, numCol].apply(lambda x : pd.to_numeric(x, "coerce"))
    # return the dataframe
    return data

def _fetchSystem_(target):
    ''' Fetch the system information for a given target planet.
        target: The name of the planet to fetch the information for.
        params: List of parameters to fetch from the stored databases. '''
    # define the allowed parameters
    params=["mass", "sma", "per", "ecc", "inc", "arg"]
    # if the target does not contain a distinction between name and number, attempt to do so
    if "-" not in target:
        numidx = target.index(findFloats(target)[0])
        target = target[:numidx]+"-"+target[numidx:]
    # fetch the central body name
    star = target[:-1]
    # add errors to the paramaters, and map to archive names
    def allParams(params):
        params = nestToLinear([[param,param+"_e1",param+"_e2"] for param in params])
        aparams = mapToArchive(params)
        return params, aparams
    param, aparam = allParams(params)
    # load the dataframe
    archive = loadDataFrame("/raw_data/pscomppars.csv")
    # fetch all planets with the required hostname
    data = archive[archive.hostname.str.replace(" ", "").isin([star, star.lower(), star.upper(), star.capitalize()])].reset_index()
    # fetch the target planet index
    targetIdx = np.where(data.pl_name.str.replace(" ", "").isin([target, target.lower(), target.upper(), target.capitalize()]))[0][0]
    # ensure it is the first value in the dataframe
    if targetIdx != 0:
        # move the specified row to the top of the dataFrame
        data = moveRowToTop(data, targetIdx)
    # construct a new dataframe of only required parameters
    systemData = data.loc[:, ["pl_name"]+aparam]
    # convert mass to solar mass
    if "pl_bmassj" in systemData.columns:
        systemData.loc[:, ["pl_bmassj", "pl_bmassjerr1", "pl_bmassjerr2"]] *= 9.547919E-4
    # convert radians to degrees
    if ("pl_orbincl" in systemData.columns) or ("pl_orblper" in systemData.columns):
        degcols = mapToArchive(["inc", "inc_e1", "inc_e2", "arg", "arg_e1", "arg_e2"])
        systemData = convert(systemData, degcols[1:], np.radians)
        # because inclination is perpendicular to earth, we need to deal with this slightly differently
        systemData = convert(systemData, degcols[0], lambda x: np.radians(90-x))
    # convert period to years
    if "pl_orbper" in systemData.columns:
        systemData = convert(systemData, mapToArchive(["per", "per_e1", "per_e2"]), lambda x:x/365.25)
    # populate a new list with the stars data
    d1 = data.iloc[0]
    starData = [d1["hostname"], d1["st_mass"], d1["st_masserr1"], d1["st_masserr2"]]
    # append it to the dataframe, filling any unspecified values with NaN
    systemData = insertRow(systemData, starData + [np.nan]*(len(systemData.iloc[0])-len(starData)), 0)
    # remove the index column if given
    if "index" in systemData.columns:
        systemData = systemData.drop(["index"], axis=1)
    # rename all the columns to be more usefull
    systemData.columns = ["name"]+param
    # return the data
    return systemData

def _possibilitySpace_(poss, tqdmLeave=True):
    ''' Construct a 2D array of all possible values, if each value can have one of two states.
        poss: A (2,N) array of possible input values, where each column is the two values for a single state.
        tqdmLeave: Whether to leave the overarchive tqdm bar when completed.'''
    # compute the maximum number of parameters
    totalParams = poss.shape[1]
    # initial possibility space array
    pos = np.full((2**totalParams, totalParams), np.nan)
    # itterate over all possibility space
    for idx in trange(2**totalParams, desc="Filling possibility space", leave=tqdmLeave):
        # convert the current itteration number to a binary with as many positions as needed
        idx_b = "{:0{}b}".format(idx, totalParams)
        # use the binary representation of this iteration to select the zeroth or first row of the output
        pos[idx,:] = np.array([poss[int(idx_b[x])][x] for x in range(totalParams)])
    # return all possibilities
    return pos

# construct all simulation possibilities
def constructSimArray(df, params=["mass", "sma", "ecc", "inc", "arg"], tqdmLeave=True, useerror=True):
    ''' Construct an array of all possible simulation parameters.
        Each row corresponds to a new simulation setup, where the columns are the parameters for each object.
        df: The dataframe of system information to fetch data from.
        params: The list of parameters to include in this simulation.
        tqdmLeave: Whether to leave the overarchive tqdm bar when completed.'''
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
        thisobj = df.loc[obj]
        # for each parameter
        for param in params:
            # fetch the value and errors
            val = thisobj[param]
            er1 = val + abs(thisobj[param+"_e1"])
            er2 = val - abs(thisobj[param+"_e2"])
            # add to the output
            out[:,idx] = np.hstack([er1, er2])
            out_[:,idx] = val
            # inrement the index counter
            idx+=1
    # if we are using error values
    if useerror:
        # fetch the entire possibility space of our combinations
        pos = _possibilitySpace_(out, tqdmLeave=tqdmLeave)
        # add the default simulation settings to the begining
        pos = np.vstack([out_, pos])
    else:
        pos = out_
    # replace nan values with infinitiy so we know they are wrong
    pos = np.nan_to_num(pos, nan=np.inf)
    # remove nonunique values
    pos = np.unique(pos, axis=0)
    # remove completely empty rows
    pos = pos[~np.isinf(pos).all(axis=1)]
    # return the possibility space
    return pos

# fetch relevant data given exoplanet
def fetchParams(exoplanet, params=["mass", "sma", "per", "ecc", "inc", "arg"], forceUse=False):
    ''' Fetch all relevant data for an exoplanet. Will raise an exception if no transit data exists for the target.
        exoplanet: The name of the exoplanet to search for.
        params: The list of parameters to be used in this simulation. '''
    # if we were given a file to simulate from
    if ".csv" in exoplanet:
        # fetch the file and planet
        file, exoplanet = exoplanet.split(",")
        # fetch the data using the manual system different handler
        archiveData = _fetchCustomSystem_(file, exoplanet)
    else:
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
            if forceUse:
                print("Mid-Transit data does not exist for target planet.")
            else:
                raise Exception("Mid-Transit data does not exist for target planet.")
        # fetch the archive data for the system
        archiveData = _fetchSystem_(archiveTarget)

    # check if sma or per are nan
    nonNan = ~archiveData.loc[1:, ["sma", "per"]].isna().any(axis=0)
    # if all semimajor axis values are non-Nan, remove the period section
    if nonNan[0]:
        params.remove("per")
    # else, if all period values are non-Nan, remove semimajor axis.
    elif nonNan[1]:
        params.remove("sma")
    # otherwise, this will just use a mixture of the two

    # return the archive dataFrame
    return archiveData, params

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

def _simulateTT_(sim, timestep, transits=1000, i=1, prec=1/31557600.0):
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
                # find the zero point by the arithmetic mean
                new_t = 0.5 * (t_0+t_1)
                # find the zero point by the geometric mean
                #new_t = np.sqrt(t_0 * t_1)
                # find the zero point by the harmonic mean
                #new_t = 1/(1/t_0 + 1/t_1)
                # step forward by the difference in our times
                sim.integrate(new_t)
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
    return tarray

def _arrayToSim_(thisSim, params):
    ''' Generates a list of simulations from an array of simulation parameters.
        thisSim: The array of simulation parameters.
        params: The list of parameter names. used to map the array to usable values. '''
    # fetch the number of objects
    objs = len(thisSim) / len(params)
    # create the simulation
    sim = rebound.Simulation()
    # set the units for the simulation
    sim.units = ('yr', 'AU', 'Msun')
    # because our units are in years, the default 1E-9 may not be small enough
    sim.ri_ias15.min_dt = 1E-9**3
    # iterate on each object
    for obj in np.split(thisSim, objs):
        # fetch this objects paramaters
        mass, sma, per, inc, ecc, arg = valFromParam(["mass", "sma", "per", "inc", "ecc", "arg"], obj, params)
        # if we have a semimajor axis
        if sma:
            sim.add(m=mass, a=sma, e=ecc, omega=arg, inc=inc)
        elif per:
            # because of our units of G=1, a=1 [AU], then t [year]=2pi
            sim.add(m=mass, P=per, e=ecc, omega=arg, inc=inc)
        else:
            # the star has no interesting properties beyond mass
            sim.add(m=mass)
    # move to the centre of mass
    sim.move_to_com()
    # add the reboundx extras
    rebx = reboundx.Extras(sim)
    # introduce general relatiity
    gr = rebx.load_force("gr")
    rebx.add_force(gr)
    gr.params["c"] = (299792458) * ((365.25*86400) / 149597870700) # speed of light in Au/Year
    # and return the simulation
    return sim

def _simulateTransitTimes_(simArray, params, transits=1000, prec=1/31557600.0):
    ''' Locate all transit times for a given simulation setup.
        simArray: The array of simulation parameters.
        params: The list of parameter names. used to map the array to usable values.
        transits: The number of transits to locate.
        prec: The precision, defaults to 1 second. '''
    # fetch this simulation
    sim = _arrayToSim_(simArray, params)
    # fetch the smallest orbital period in the system
    p_orb = min([x.P for x in sim.particles[1:]])
    # ensure our timestep is resonable: Take the smallest orbital period:
    # if it is very small, step in units of 0.1 periods.
    # if very large, step in units of 0.001 periods.
    # if niehter, step hourly
    timestep = max(p_orb*0.001, min(p_orb*0.1, 1/8760))
    # fetch the position from the worker ID
    pos = curProc()._identity[0] if curProc().name != "MainProcess" else None
    # simulate for the chosen number of transits
    TT = _simulateTT_(sim, timestep, transits, pos, prec)
    # and return the transit times
    return np.array(TT)

def fetchTT(simArray, params, transits=1000, prec=1/31557600.0, workers=None, tqdmLeave=True, returnAll=False):
    ''' Create a multiprocessing job to simulate all possible configurations of a planetary system.
        returns the transit times, the upper, and the lower error bounds.
        simArray: The 2D array of all possible simulation parameters.
        params: The list of parameter names. used to map the array to usable values.
        transits: The number of transits to simulate for.
        prec: The precision, defaults to 1 second.
        workers: The number of worker processes to use for the simulation (defaults to cpucount / 2).
        tqdmLeave: Whether to leave the overarchive tqdm bar when completed.
        returnAll: Whether to return the entire output, or just the average/minimum/maximum.'''
    inputs = toItterableInput(simArray, params, transits, prec, keep=(1,))
    # run the multiprocessing job
    TT = parallelJob(_simulateTransitTimes_, inputs, outType=np.array, workers=workers, tqdmLeave=tqdmLeave)
    # if we want to return everything, do that
    if returnAll:
        return TT
    return avMinMax(TT, 0)

def _planetVars_(star, planet, mTot):
    # gravitational parameter
    Gmu = (planet["mass"]+star["mass"])*6.67408E-11*1.98855E30
    # we require a period or semimajor axis
    if ("sma" not in planet) and ("per" not in planet):
        raise Exception("{} requires a defined semimajor axis, or period. Has neither!".format(planet["name"]))
    hassma, hasper = ~np.isnan(planet["sma"]), ~np.isnan(planet["per"])
    # if we have a semimajor axis
    if hassma:
        # convert to metres
        a = planet["sma"]*149597870700
        # if we don't have an orbital period, compute
        if ~hasper:
            p = 2*np.pi*np.sqrt(a**3 / Gmu)
    # if we have a period
    if hasper:
        # convert to seconds
        p = planet["per"]*31557600
        # if we didn't have a semimajor axis
        if ~hassma:
            a = (Gmu * p**2 / (4*np.pi**2))**(1/3)
    # fetch everything else
    arg, inc, ecc, mMu = np.radians(planet["arg"]), np.radians(90-planet["inc"]), planet["ecc"], planet["mass"] / mTot
    # variables in an np array
    var = np.array([a, p, mMu, ecc, inc, arg])
    # remove NaN
    var[np.isnan(var)] = 0
    # and return the variables
    return var

def _innerTTV_(target, perturber):
    ''' Compute the predicted magnitude of a TTV from an inner perturbing planet.
        target: The transiting planet to observe for TTV of.
        perturber: The perturbing planet that will cause a TTV.'''
    # decompose variables
    a1, p1, m1, e1, i1, arg1 = perturber
    a2, p2, m2, e2, i2, arg2 = target
    # compute the predicted TTV
    return (p2 * m1 * a1 * np.cos(arg1) * np.sqrt(1 - e2**2)) / (2*np.pi * a2 * (1 + e2*np.sin(arg2)))

def _outerTTV_(target, perturber):
    ''' Compute the predicted magnitude of a TTV from an outer perturbing planet.
        target: The transiting planet to observe for TTV of.
        perturber: The perturbing planet that will cause a TTV.'''
    # decompose variables
    a1, p1, m1, e1, i1, arg1 = target
    a2, p2, m2, e2, i2, arg2 = perturber
    # compute the predicted TTV
    return m2 * e2 * p2 * (a1 / a2)**3

def predictTTVMagnitude(df):
    ''' Calculate the predicted TTV Magnitude due to planets in the system.
        df: The dataframe holding the system information. '''
    # compute the reduced mass for each object
    mTotal = sum(df["mass"])
    # and fetch the rquired variables
    planets = [_planetVars_(df.iloc[0], df.iloc[idx], mTotal) for idx in df.index[1:]]
    # reference to the target planet
    target, TTV = planets[0], []
    # compute the TTV due to each planet
    for otherPlanet in planets[1:]:
        # if the semimajor axis of the target is larger
        if target[0] > otherPlanet[0]:
            TTV.append(_innerTTV_(target, otherPlanet))
        # else
        else:
            TTV.append(_outerTTV_(target, otherPlanet))
    # show an inacuracy warning if there are multiple planets in the system
    if len(planets) > 2:
        print("Warning! TTV Prediction may be innacurate due to system containing more than two planets!")
    # compute the TTV by summation
    #print(np.array(TTV))
    #print(np.average(1 / np.array(TTV)))
    #print(sum(TTV))
    TTV = sum(TTV)
    # return the TTV prediction
    return TTV
