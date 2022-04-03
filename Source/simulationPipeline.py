# python modules
from tqdm import tqdm
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
    systemData.columns = ["name", "mass", "masse1", "masse2", "per", "pere1", "pere2",
                          "ecc", "ecce1", "ecce2", "inc", "ince1", "ince2", "arg", "arge1", "arge2"]
    return systemData

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

    # fetch the archive data for the planet
    archiveData = fetchSystem(archiveTarget)
    print(archiveData)

fetchParams("HAT-P-13 b")
fetchParams("Epic-203771098 b")
