# python modules
import pandas as pd
import numpy as np

# my modules
import TransitProject as tp
import TransitProject.Simulation as ts

'''
Goals: Determine which, if any, of the known exolanets are suitable candidates for TTV control tests.
       They must:
         be in a planetary system with 1 other planet.
         have known midtransit times.
'''

def checkTransiting(systemname, systemdf):
    ''' check if any planets in a system have midtransit data. '''
    # fetch the system
    thissystem = systemdf[systemdf["hostname"]==systemname]
    # and fetch each planet
    theseplanets = list(thissystem["pl_name"].str.capitalize().str.replace(" ",""))
    # check if either planet has midtransit times
    transits = [tp.inDataFrame("/raw_data/exoplanetList", planet) for planet in theseplanets]
    # if they do, return the planet name(s)
    return list(np.array(theseplanets)[transits])

def searchIdeal(systemdf):
    ''' search through a dataframe to find how many planetary systems are ideal TTV control candidates. '''
    # look through the hostname section
    stars = sorted(systemdf["hostname"].unique())
    # for each star, compute the planetary systsem size
    sysSize = np.array([len(np.where(systemdf["hostname"]==star)[0]) for star in stars])
    # pull out ideal and nonideal candidates
    single = np.where(sysSize==1)[0]
    ideal = np.where(sysSize==2)[0]
    nonideal = np.where(sysSize>2)[0]
    # fetch the transiting systems in the ideal and nonideal list
    singleTransiting = [checkTransiting(stars[star], systemdf) for star in single]
    idealTransiting = [checkTransiting(stars[star], systemdf) for star in ideal]
    nonidealTransiting = [checkTransiting(stars[star], systemdf) for star in nonideal]
    # remove empty values and flatten
    singleTransiting = sum([transit for transit in singleTransiting if transit], [])
    idealTransiting = sum([transit for transit in idealTransiting if transit], [])
    nonidealTransiting = sum([transit for transit in nonidealTransiting if transit], [])
    # return
    return singleTransiting, idealTransiting, nonidealTransiting

def midTransitTimeNumber(planets):
    ''' compute the number of midtransit times for each planet given. '''
    # dictionary of planet times
    planetTimes = {}
    # for each planet
    for planet in planets:
        # load the midtransit times
        df = tp.loadDataFrame("/raw_data/midTransitTimes/{}".format(planet))
        # add the number of midtransit times to the planet dictionary
        planetTimes[planet] = len(df)
    # return the dictionary
    return planetTimes

if __name__ == "__main__":
    # load the exoplanet archive
    systemdf = tp.loadDataFrame("/raw_data/pscomppars")
    # search the exoplanet archive for ideal planetary systems
    single, ideal, nideal = searchIdeal(systemdf)
    print("Of the planets with mid-transit times:")
    print("{} are in single-planet systems".format(len(single)))
    print("{} are in two planet systems".format(len(ideal)))
    print("{} are in multi-planet systems".format(len(nideal)))
    # fetch the number of midtransit times per planet
    mtts = midTransitTimeNumber(single)
    mtti = midTransitTimeNumber(ideal)
    mttn = midTransitTimeNumber(nideal)
    # sort by value (second value in dict.items)
    mtts = dict(sorted(mtts.items(), key=lambda x:x[1], reverse=True))
    mtti = dict(sorted(mtti.items(), key=lambda x:x[1], reverse=True))
    mttn = dict(sorted(mttn.items(), key=lambda x:x[1], reverse=True))
    # show as a nice output
    print("\nThe TTV search candidates in singe-planet systems (in decreasing order of datapoints) are:")
    print(",\t".join(["{}: {}".format(k,v) for k,v in mtts.items()]))
    print("\nThe TTV test candidates in 2-planet systems (in decreasing order of datapoints) are:")
    print(",\t".join(["{}: {}".format(k,v) for k,v in mtti.items()]))
    print("\nThe TTV test candidates in 2+-planet systems (in decreasing order of datapoints) are:")
    print(",\t".join(["{}: {}".format(k,v) for k,v in mttn.items()]))
    # combine dictionaries into 2d array of key and value
    combined = np.array(list(mtts.items()) + list(mtti.items()) + list(mttn.items()))
    # create a dataframe
    controlCandidates = pd.DataFrame.from_dict({"Exoplanet":combined[:, 0], "Datapoints":combined[:, 1]})
    # and save it
    tp.saveDataFrame(controlCandidates, "/ControlCandidates")
