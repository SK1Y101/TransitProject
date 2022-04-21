# python modules
import matplotlib.pylab as plt
import scipy.optimize
import pandas as pd
import numpy as np
import argparse, os

# my modules
import TransitProject as tp
import TransitProject.Simulation as ts
import processing

#open the sims
files = ["/sim_data/"+file for file in os.listdir(os.getcwd()+"/sim_data/") if "TTV" in file and "period" in file]

maxval=0
# for each file
for file in files:
    # load them
    TTV = tp.loadDataFrame(file)
    # fetch the average TTV
    TTV = TTV["0"].iloc[:400]
    x = np.arange(len(TTV))
    # fit a curve to it
    fit = processing.fetchFit(x, TTV, f=lambda x,a,b,:a*x+b)[0]
    TTV = TTV-fit(x)
    # plot
    plt.plot(TTV, label=file)
    maxval = max(max(abs(TTV)), maxval)
plt.legend()
plt.ylim((-np.ceil(maxval)-1, np.ceil(maxval)+1))
plt.show()
