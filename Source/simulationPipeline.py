# python modules
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

def fetchSystem(target):
    # fetch the central body name
    star = target[:-1]
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
    systemData = data.loc[:, ["pl_name",
                              # mass*sin(i) or mass in jupiter mass, depending on if sin(i) is known
                              "pl_bmassj", "pl_bmassjerr1", "pl_bmassjerr2",
                              # orbital period in days
                              "pl_orbper", "pl_orbpererr1", "pl_orbpererr2",
                              # eccentricity
                              "pl_orbeccen", "pl_orbeccenerr1", "pl_orbeccenerr2",
                              # inclination
                              "pl_orbincl", "pl_orbinclerr1", "pl_orbinclerr2",
                              # argument of periapsis
                              "pl_orblper", "pl_orblpererr1", "pl_orblpererr2"]]
    # convert mass to solar mass
    systemData.loc[:, ["pl_bmassj", "pl_bmassjerr1", "pl_bmassjerr2"]] *= 9.55E-4
    # convert inclination and arg. per. to radians
    systemData.loc[:, ["pl_orbincl", "pl_orbinclerr1", "pl_orbinclerr2",
                       "pl_orblper", "pl_orblpererr1", "pl_orblpererr2"]] *= (np.pi/180)
    # populate a new list with the stars data
    d1 = data.iloc[0]
    starIdx, starData = len(systemData.index), [d1["hostname"], d1["st_mass"], d1["st_masserr1"], d1["st_masserr2"]]
    # append it to the dataframe, filling any unspecified values with NaN
    systemData.loc[starIdx] = starData + [np.nan]*(len(systemData.iloc[0])-len(starData))
    # and shift it to the top
    systemData = moveRowToTop(systemData, starIdx)
    # rename all the columns to be more usefull
    systemData.columns = ["name", "mass", "mass_e1", "mass_e2", "per", "per_e1", "per_e2",
                          "ecc", "ecc_e1", "ecc_e2", "inc", "inc_e1", "inc_e2", "arg", "arg_e1", "arg_e2"]
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
    # replace nan values with infinitiy so we know they are wrong
    pos = np.nan_to_num(pos, nan=np.inf)
    # fetch and return the unique permutations
    return np.unique(pos, axis=0)

# construct all simulation possibilities
def constructSimArray(df):
    # parameters needed for the simulation
    params = ["mass", "per", "ecc", "inc", "arg"]
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

    # return the possibility space
    return pos

# fetch relevant data given exoplanet
def fetchParams(exoplanet):
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
    archiveData = fetchSystem(archiveTarget)

    # return the archive dataFrame
    return archiveData

# testing area
df = fetchParams("HAT-P-13 b")

simArray = constructSimArray(df)
print(simArray.shape)

simArray = constructSimArray(df.iloc[:2])
print(simArray.shape)
