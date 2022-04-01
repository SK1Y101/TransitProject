# python modules
from datetime import datetime as dt
import matplotlib.pylab as plt
import pandas as pd
import numpy as np

# my modules
import TransitProject as tp
import TransitProject.webScraping as ws

planets = tp.loadDataFrame("/raw_data/pscomppars.csv")
targetName = "Hat-p-13 b"
target = planets.loc[planets["hostname"].str.lower() == targetName.lower().split()[0]]
print(target)

#print(target)
#print(target["pl_orbper"], target["pl_orbpererr1"], target["pl_orbpererr2"])

#print(target["pl_orbpererr1"] * 365 / target["pl_orbper"] * 10 * 1440)

# plot a timeseries of data
df = tp.loadDataFrame("/raw_data/midTransitTimes/"+targetName.replace(" ","").capitalize())

for source in sorted(df["source"].unique()):
    data = df.loc[df["source"]==source]
    data = data.sort_values(by="date", ascending=True)
    plt.errorbar(pd.to_datetime(data["date"]), data["oc"], yerr=data["oce"], label=source, fmt="o")

plt.gcf().autofmt_xdate()
plt.legend()
plt.show()
