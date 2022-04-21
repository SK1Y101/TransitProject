# python modules
import matplotlib.pylab as plt
import scipy.optimize
import pandas as pd
import numpy as np
import argparse, os

# my modules
import TransitProject as tp
import TransitProject.Simulation as ts

#open the sims
files = ["/sim_data/"+file for file in os.listdir(os.getcwd()+"/sim_data/") if "TTV" in file and "eccentricity" in file]

# for each file
for file in files:
    # load them
    TTV = tp.loadDataFrame(file)
    # fetch the average TTV
    TTV = TTV["0"]
    # plot
    plt.plot(TTV.iloc[:400], label=file)
plt.legend()
plt.ylim((-30, 30))
plt.show()
