# import python modules
from tqdm import tqdm
import pandas as pd
import numpy as np

# import my modules
import TransitProject as tp
import TransitProject.webScraping as ws

# fetch the ephemerides for all transiting planets, and all planets in the same system
planets = tp.loadDataFrame("/raw_data/exoplanetList.csv")
exoclock = tp.loadDataFrame("/raw_data/exoClockEphemerides.csv")
archive = tp.loadDataFrame("/raw_data/pscomppars.csv")
for x in tqdm(planets["name"]):
    planet = x.lower()
    star = x[:-1].upper()
    print(x)
    print(exoclock.loc[exoclock["name"].str.lower()==planet])
    print(archive.loc[archive["hostname"].str.replace(" ","") == star])
