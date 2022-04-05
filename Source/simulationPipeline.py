# python modules
from multiprocessing.pool import Pool
from tqdm import tqdm, trange
import numpy as np
import rebound

# my modules
import TransitProject as tp
import TransitProject.webScraping as ws

def moveRowToTop(df, targetIdx):
    # append an empty row
    df.loc[len(df.index)] = np.nan
    # shift down
    df = df.shift(1)
    # copy target to first row
    df.iloc[0] = df.iloc[targetIdx+1]
    # remove the old value and return
    return df[df.index != targetIdx+1]


# map my parameter names to archive names
def mapToArchive(params):
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

def fetchSystem(target, params):
    # fetch the central body name
    star = target[:-1]
    # add errors to the paramaters, and map to archive names
    params = ",".join([",".join([param,param+"_e1",param+"_e2"]) for param in params]).split(",")
    aparams = mapToArchive(params)
    # load the dataframe
    archive = tp.loadDataFrame("/raw_data/pscomppars.csv")
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
    return systemData

# construct a 2d array of all possible values
# takes a (2,N) array, where the two different possible values for a column are (1, N) and (2, N)
def possibilitySpace(out, totalParams):
    # initial possibility space array
    pos = np.full((2**totalParams, totalParams), np.nan)
    # itterate over all possibility space
    for idx in trange(2**totalParams, desc="Filling possibility space"):
        # convert the current itteration number to a binary with as many positions as needed
        idx_b = "{:0{}b}".format(idx, totalParams)
        # use the binary representation of this iteration to select the zeroth or first row of the output
        pos[idx,:] = np.array([out[int(idx_b[x])][x] for x in range(totalParams)])
    # return all possibilities
    return pos

# construct all simulation possibilities
def constructSimArray(df, params=["mass", "sma", "ecc", "inc", "arg"]):
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
    pos = possibilitySpace(out, totalParams)
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
def fetchParams(exoplanet, params):
    # convert to capitalised with no space for exopclock
    exoclockTarget = exoplanet.capitalize().replace(" ", "")
    # convert name to uppercase (leaving planet delimiter lowercase) and remove spaces for the archive
    if exoplanet[-1].isalpha():
        archiveTarget = exoplanet[:-1].upper().replace(" ", "")+exoplanet[-1].lower()
    else:
        # if there wasn't a planet delimiter, add one
        archiveTarget = exoplanet.upper().replace(" ", "")+"b"

    # check the planet has midtransit data
    if not tp.inDataFrame("/raw_data/exoplanetList.csv", exoclockTarget):
        raise Exception("Mid-Transit data does not exist for target planet.")
    # fetch the archive data for the system
    archiveData = fetchSystem(archiveTarget, params)

    # return the archive dataFrame
    return archiveData

# if a value is infinity, return None, else return the value
def noneIfInf(x):
    if np.isinf(x):
        return None
    return x

def valFromParam(value, array, param):
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

def fetchSims(simArray, params):
    # fetch the number of objects
    objs = int(len(simArray[0]) / len(params))
    # simulation list
    sims = []
    # construct the required simulations
    for thisSim in tqdm(simArray, desc="Creating simulations"):
        # create the base simulation
        sim = rebound.Simulation()
        i=0
        # for each object
        for obj in np.split(thisSim, objs):
            # fetch the parameters for this object
            mass, sma, per, inc, ecc, arg = valFromParam(["mass", "sma", "per", "inc", "ecc", "arg"], obj, params)
            # add it to the simulation
            if sma:
                # default to using sma if possible
                sim.add(m=mass, a=sma, e=ecc, omega=arg, inc=inc)
            else:
                # else use period
                sim.add(m=mass, P=per, e=ecc, omega=arg, inc=inc)
            i+=1
        # move the simulation to the centre of mass
        sim.move_to_com()
        # and append to the simulation list
        sims.append(sim)
    # return the simulations
    return sims

def simulateTT(sim, timestep, transits=1000, i=0, prec=1E-7):
    # transit time array
    tarray, n = np.zeros(transits), 0
    # reference to the particles in the simulation
    p = sim.particles
    # create a progress bar for this simulation
    pbar = tqdm(total=transits, desc="Simulating transits", leave=False)
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
    return tarray

def fetchTT(sims, transits=1000):
    # fetch the orbital period of our transiting target
    p_orb = sims[0].particles[1].P
    # step either hourly, or per tenth of an orbit, depending on which is smaller
    tstep = min(p_orb*0.1, 1 / (24*365))
    # empty array of transit times
    TT = np.zeros((len(sims), transits))
    i=0
    # fetch all the simulated transit times
    for sim in tqdm(sims, desc="Simulating"):
        # fetch this transit time
        tt = simulateTT(sim, tstep, transits, sims.index(sim))
        # and add to the transit time array
        TT[i, :] = tt
        i+=1

    # compute the average time, mininimum time, and maximum time for each transit
    av = TT[0]
    mn = TT.min(axis=0)
    mx = TT.max(axis=0)

    # return the average transit time, and the error values
    return av, np.vstack([av-mn, mx-av])

# testing area

# define the simulation parameters
''' choice of mass, period or semimajor axis, eccentiricty, inclination, argument of periapsis '''
''' if both a period and semimajor axis is given, the script will default to SMA '''
params = ["mass", "sma"]#, "ecc"]#, "inc"]

N=1000


import matplotlib.pylab as plt

# fetch the parameters for the system
df = fetchParams("HAT-P-13 b", params)

simArray = constructSimArray(df, params)

sims = fetchSims(simArray, params)

TT = fetchTT(sims, N)
A = np.vstack([np.ones(N), range(N)]).T
c, m = np.linalg.lstsq(A, TT[0], rcond=-1)[0]

plt.errorbar(range(N), TT[0], yerr=TT[1], label="From archive", fmt="o")
plt.legend()
plt.show()
