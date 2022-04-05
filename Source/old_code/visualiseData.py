# python modules
from datetime import datetime as dt
import matplotlib.pylab as plt
from tqdm import tqdm
import pandas as pd
import numpy as np

# my modules
import TransitProject as tp
import TransitProject.webScraping as ws

planetNames = tp.loadDataFrame("/raw_data/exoplanetList.csv")
planets = tp.loadDataFrame("/raw_data/pscomppars.csv")
for x in tqdm(planetNames["name"]):
    targetName = x
    target = planets.loc[planets["hostname"].str.replace(" ","") == targetName.upper()[:-1]]
    #print(target)

    #print(target)
    #print(target["pl_orbper"], target["pl_orbpererr1"], target["pl_orbpererr2"])
    #print(target["st_mass"])

    #print(target["pl_orbpererr1"] * 365 / target["pl_orbper"] * 10 * 1440)

    # plot a timeseries of data
    df = tp.loadDataFrame("/raw_data/midTransitTimes/"+targetName.replace(" ","").capitalize())

    sources = sorted(df["source"].unique())
    if "ETD" in sources:
        sources.remove("ETD")

    if sources:
        # fetch the data
        for source in sources:
            data = df.loc[df["source"]==source]
            data = data.sort_values(by="date", ascending=True)
            plt.errorbar(pd.to_datetime(data["date"]), data["oc"], yerr=data["oce"], label=source, fmt="o")

        #dates = pd.to_datetime(df["date"])
        #plt.plot([min(dates), max(dates)], [0,0])

        plt.title(targetName)
        plt.gcf().autofmt_xdate()
        plt.legend()
        plt.show()
