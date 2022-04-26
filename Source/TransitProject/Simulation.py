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
                "arg":"pl_orblper",
                "mean":"pl_orbtper"}
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
    params=["mass", "sma", "per", "ecc", "inc", "arg", "mean"]
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

def _fetchSystem_(target, startTime="2000-01-01"):
    ''' Fetch the system information for a given target planet.
        target: The name of the planet to fetch the information for.
        params: List of parameters to fetch from the stored databases. '''
    # define the allowed parameters
    params=["mass", "sma", "per", "ecc", "inc", "arg", "mean"]
    # if the target does not contain a distinction between name and number, attempt to do so
    #if "-" not in target:
    #    numidx = target.index(findFloats(target)[0])
    #    target = target[:numidx]+"-"+target[numidx:]
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
    data = archive[archive.hostname.str.replace(" ", "").str.lower().isin([star.lower()])].reset_index()
    # fetch the target planet index
    targetIdx = np.where(data.pl_name.str.replace(" ", "").str.lower().isin([target.lower()]))[0][0]
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
    # fetch the epoch of periapse passage and starting time
    EoPP = HJDtoDate(systemData["mean"])
    EoPP_e1 = HJDtoDate(systemData["mean"]+systemData["mean_e1"])
    EoPP_e2 = HJDtoDate(systemData["mean"]+systemData["mean_e2"])
    startTime = pd.to_datetime(startTime)
    # compute the difference from starting time
    tdif, tdif_e1, tdif_e2 = EoPP - startTime, EoPP_e1 - startTime, EoPP_e2 - startTime
    # compute the mean anomaly due to the start time
    m = ((tdif.dt.total_seconds() / (systemData["per"] * 31557600)) % 1) * 2 * np.pi
    m_e1 = ((tdif_e1.dt.total_seconds() / ((systemData["per"]+systemData["per_e2"]) * 31557600)) % 1) * 2 * np.pi
    m_e2 = ((tdif_e2.dt.total_seconds() / ((systemData["per"]+systemData["per_e1"]) * 31557600)) % 1) * 2 * np.pi
    # and change to mean anomaly
    systemData["mean"] = m
    systemData["mean_e1"] = m_e1
    systemData["mean_e2"] = m_e2
    # return the data
    return systemData

def _possibilitySpace_(poss, tqdmLeave=True):
    ''' Construct a 2D array of all possible values, if each value can have one of two states.
        poss: A (2,N) array of possible input values, where each column is the two values for a single state.
        tqdmLeave: Whether to leave the overarchive tqdm bar when completed.'''
    # compute the maximum number of parameters
    totalParams = poss.shape[1]
    # compute the size (In bytes) of a (2**totalParams, totalParams) array
    asize = 2**totalParams * totalParams * 8
    # if the array would be too large for memory, do the slower but more memory efficient search
    if asize > 4*2**30:
        # compute the array size in powers of two
        print("Fast approach requires {:0.3f} {} of memory, using slow approach instead".format(*tosi(asize, "B", 1024)))
        # initialise as empty
        pos = []#np.zeros((2, totalParams))
        # itterte
        for idx in trange(2**totalParams, desc="Filling possibility space", leave=tqdmLeave):
            # convert the current itteration number to a binary with as many positions as needed
            idx_b = "{:0{}b}".format(idx, totalParams)
            # use the binary representation of this iteration to select the zeroth or first row of the output
            val = np.array([poss[int(idx_b[x])][x] for x in range(totalParams)])
            # remove nan
            val = np.nan_to_num(val, nan=np.inf)
            # if this doesnt match any previous rows
            if list(val) not in pos:#~np.isin(pos, val).all(axis=1).any():
                # add to the possibilities area
                #pos = np.vstack([pos, val])
                pos.append(list(val))
        # remove the initial zero
        pos = np.array(pos)
    # otherwise
    else:
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
def constructSimArray(df, params=["mass", "sma", "ecc", "inc", "arg", "mean"], tqdmLeave=True, useerror=True, limitError=False):
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
            # if we are limiting errors, and this is an error to limit
            if limitError and param not in ["mass", "sma", "per"]:
                er1, er2 = np.nan, np.nan
            # add to the output
            out[:,idx] = np.hstack([er1, er2])
            out_[:,idx] = val
            # inrement the index counter
            idx+=1
    # remove NaN for default
    out_ = np.nan_to_num(out_, nan=np.inf)
    # if we are using error values
    if useerror:
        # fetch the entire possibility space of our combinations
        pos = _possibilitySpace_(out, tqdmLeave=tqdmLeave)
        # replace nan values with infinitiy so we know they are wrong
        pos = np.nan_to_num(pos, nan=np.inf)
        # remove completely empty rows
        pos = pos[~np.isinf(pos).all(axis=1)]
        # remove nonunique values
        pos = np.unique(pos, axis=0)
        # add the default simulation settings to the begining
        pos = np.vstack([out_, pos])
    else:
        # set the possibility space to only the default values
        pos = out_
    # return the possibility space
    return pos

# fetch relevant data given exoplanet
def fetchParams(exoplanet, params=["mass", "sma", "per", "ecc", "inc", "arg", "mean"], forceUse=False, startTime="2000-01-01"):
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
        # if the target does not contain a distinction between name and number, attempt to do so
        #if "-" not in exoplanet:
        #    numidx = exoplanet.index(findFloats(exoplanet)[0])
        #    exoplanet = exoplanet[:numidx]+"-"+exoplanet[numidx:]
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
            err = "Mid-Transit data does not exist for {}".format(exoclockTarget)
            if forceUse:
                print(err)
            else:
                raise Exception(err)
        # fetch the archive data for the system
        archiveData = _fetchSystem_(archiveTarget, startTime=startTime)

    # check if sma or per are nan
    nonNan = ~archiveData.loc[1:, ["sma", "per"]].isna().any(axis=0)
    # if all semimajor axis values are non-Nan, remove the period section
    if nonNan[0]:
        if "per" in params:
            params.remove("per")
    # else, if all period values are non-Nan, remove semimajor axis.
    elif nonNan[1]:
        if "sma" in params:
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
        mass, sma, per, inc, ecc, arg, mean = valFromParam(["mass", "sma", "per", "inc", "ecc", "arg", "mean"], obj, params)
        # if we have a semimajor axis
        if sma:
            sim.add(m=mass, a=sma, e=ecc, omega=arg, inc=inc, M=mean)
        elif per:
            # because of our units of G=1, a=1 [AU], then t [year]=2pi
            sim.add(m=mass, P=per, e=ecc, omega=arg, inc=inc, M=mean)
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
    #sim.cite()
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
    # if neither, step hourly
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
    # fetch the inputs for the multiprocessing
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

def _resonanceOrder_(target, perturber):
    ''' Compute the resonance order (j) of two planets.
        target: The transiting planet to observe for TTV of.
        perturber: The perturbing planet that will cause a TTV.'''
    # decompose variables
    p1 = target[1]
    p2 = perturber[1]
    # greatest common divisor
    gcd = np.gcd(int(p1), int(p2))
    # decompose to (any) common multiples
    j = int(p1 / gcd), int(p2 / gcd)
    # an nth order resonance has value m:m-n, and has j=(m-n) * n ie: n=1, first order resonance
    order = max(j) - min(j)
    # return
    return order, min(j) * order

def _resonantTTV_(target, perturber):
    ''' Compute the predicted magnitude of a TTV from a perturbing planet in a resonant orbit.
        target: The transiting planet to observe for TTV of.
        perturber: The perturbing planet that will cause a TTV.'''
    # decompose variables
    a1, p1, m1, e1, i1, arg1 = target
    a2, p2, m2, e2, i2, arg2 = perturber
    # resonance order.
    j = _resonanceOrder_(target, perturber)[1]
    # compute predicted TTV
    return (p2 / (4.5*j)) * (m1 / (m1+m2))

def predictTTVMagnitude(df):
    ''' Calculate the predicted TTV Magnitude due to planets in the system.
        df: The dataframe holding the system information. '''
    # compute the reduced mass for each object
    mTotal = sum(df["mass"])
    # and fetch the rquired variables
    planets = np.array([_planetVars_(df.iloc[0], df.iloc[idx], mTotal) for idx in df.index[1:]])
    # reference to the target planet
    target, TTV, TTVt = planets[0], [], []
    # compute the TTV due to each planet
    for otherPlanet in planets[1:]:
        # compute any resonances
        TTV.append(_resonantTTV_(target, otherPlanet))
        # and the resonance timescale
        TTVt.append(int(target[1]) * int(otherPlanet[1]) / np.gcd(int(target[1]), int(otherPlanet[1])))
        # if the semimajor axis of the target is larger
        if target[0] > otherPlanet[0]:
            TTV.append(_innerTTV_(target, otherPlanet))
        # else
        else:
            TTV.append(_outerTTV_(target, otherPlanet))
        # the orbital period
        TTVt.append(otherPlanet[1])
    # show an inacuracy warning if there are multiple planets in the system
    if len(planets) > 2:
        print("Warning! TTV Prediction may be innacurate due to system containing more than two planets!")
    # compute which TTV effect was largest
    largesteffect = np.array(TTV).argmax()
    # return the TTV prediction and the timescale of the largest effect
    return float(sum(TTV))/86400, TTVt[largesteffect]/86400

def _timeFromTrueAnomaly_(f, t, e=0):
    ''' compute the time of flight from periapse to a given true anomaly.
        f: true anomaly (radians)
        t: orbital period.
        e: eccentricity. '''
    # compute eccentric anomaly
    E = np.arccos((e+np.cos(f)) / (1+e*np.cos(f)))
    # compute the mean anomaly
    M = E - e*np.sin(E)
    # compute the mean motion
    n = t / (2*np.pi)
    # and compute this time
    return M / n

def predictGREffect(df, observationYears=4):
    ''' Compute the TTV magnitude due to general relativity.
        df: The system data dataframe.
        observationYears: Self explanatory. '''
    # fetch the variables needed (use zero for total mass as we don't acually care about mass for this)
    a, t, _, e, _, _ = _planetVars_(df.iloc[0], df.iloc[1], 1)
    # compute the number of transits
    Ntransits = int(np.ceil(31557600 * observationYears/t))
    # speed of light
    c = 299792458
    # compute epsilon (this is in radians of precession per orbit)
    epsilon = 24 * np.pi**3 * a**2 * t**-2 * c**-2 * (1-e**2)**-1
    # the timescale this effect is cyclic is a full circle / epsilon
    cycle_time = (np.pi * 2 / epsilon) * t
    # convert from seconds to days
    cycle_time /= 86400
    # this twists the apsides by epsilon each orbit, which is equivalent to saying the true anomaly of the transit is
    # offset by epsilon each orbit. compute the time from periapse to epsilon times the number of observed transits
    # , also return the cycle time
    return _timeFromTrueAnomaly_(epsilon*Ntransits, t, e), cycle_time

def predictBCEffect(df):
    ''' compute the TTV magnitude due to other planets in the system shifting the barycentre position.
        df: the system data dataframe.'''
    # fetch the total mass and reduced mass of the star
    mTotal = sum(df["mass"])
    mStar = df.iloc[0]["mass"] / mTotal
    # fetch the planetary system variables
    planets = np.array([_planetVars_(df.iloc[0], df.iloc[idx], mTotal) for idx in df.index[1:]])
    # fetch the centre of mass offset due to all other planets in the system.
    com = sum(planets[1:, 0] * (1+planets[1:, 3]) * planets[1:, 2]) / (sum(planets[1:, 2]) + mStar)
    # compute the true anomaly the transiting planet must travel to correct for the com shift
    com_ang = np.arctan2(com, planets[0][0] * (1+planets[0][3]))
    # compute the time this is cyclic over (synodic period of all planets in the system)
    cycleTime = 1/max(1 / planets[1:, 1])-min(1 / planets[1:, 1])
    # if there was only one planet
    cycleTime = planets[1, 1] if cycleTime==0 else cycleTime
    # convert from seconds to days
    cycleTime /= 86400
    # and compute the time of flight for this angle
    return _timeFromTrueAnomaly_(com_ang, planets[0][1], planets[0][3]), cycleTime

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

def plotSystem(sim, df, savefig=False):
    ''' Plot the layout of a system.
        sim: The simulation to plot.
        df: The dataframe of the system. '''
    from matplotlib.lines import Line2D
    import matplotlib.pylab as plt
    # plot the orbits
    fig, ax_main, ax_top, ax_right = rebound.OrbitPlot(sim, slices=0.5, xlim=[-2.,2], ylim=[-2.,2], color=True, unitlabel="[AU]")
    # use the stars name as tge system name
    fig.suptitle(df.iloc[0]["name"])
    # colours that rebound uses
    colours = [(1.,0.,0.),(0.,0.75,0.75),(0.75,0.,0.75),(0.75, 0.75, 0,),(0., 0., 0.),(0., 0., 1.),(0., 0.5, 0.)]
    # create a new line element for each planet
    legend_elements = [Line2D([0], [0], color=colours[idx-1], lw=1.5, label="{} Orbit".format(df.iloc[idx]["name"])) for idx in df.index[1:]]
    # add to the legend
    ax_main.legend(handles = legend_elements, loc=2)
    # create the tight layout
    fig.tight_layout()
    # save the figure if desired
    if savefig:
        fig.savefig("{} layout.png".format(df.iloc[0]["name"]), transparent=True, bbox_inches='tight')
    # show the figure
    plt.show()
