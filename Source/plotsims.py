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

#open the TTV and sims
parser = argparse.ArgumentParser(description="Plot simulations and models")
parser.add_argument("--files", help="The names of files to use", nargs="+")
args = parser.parse_args()

allowedfiles = lambda x : True in [file in x for file in args.files]

files = ["/sim_data/"+file for file in os.listdir(os.getcwd()+"/sim_data/") if "TTV" in file and allowedfiles(file)]
systems = ["/sim_data/"+file for file in os.listdir(os.getcwd()+"/sim_data/") if "system" in file and allowedfiles(file)]

maxval=0
# for each file
for file in files:
    # load them
    TTV = tp.loadDataFrame(file)
    # fetch the average TTV
    TTV = TTV["0"]#.iloc[:1000]
    x = np.arange(len(TTV))
    # fit a curve to it
    fit = processing.fetchFit(x, TTV, f=lambda x,a,b,:a*x+b)[0]
    TTV = TTV-fit(x)
    # plot
    plt.plot(TTV, label=file)
    maxval = max(max(abs(TTV)), maxval)

for system in systems:
    sys = tp.loadDataFrame(system)
    # fetch the params
    target, perturber = ts.planetParams(sys)[:2]
    # fetch the outer planet index
    a1, p1, mu1, e1, i1, arg1 = target
    a2, p2, mu2, e2, i2, arg2 = perturber
    # and compute the TTV analytically
    def inner_dt(a1, a2, p1, p2, mu1, mu2, tc=0):
        return -a2 * mu2 * p1 * np.sin(2*np.pi * tc / p2) / (2*np.pi*a1)
    def outer_dt(a1, a2, p1, p2, mu1, mu2, tc=0):
        return -a1 * mu1 * p2 * np.sin(2*np.pi * tc / p1) / (2*np.pi*a2)
    # convert from transit number to time
    x = np.arange(len(TTV), step=.1)
    TTVpi = inner_dt(a1, a2, p1, p2, mu1, mu2, x * p1)
    TTVpo = outer_dt(a1, a2, p1, p2, mu1, mu2, x * p2)
    print("Inner Perturber TTV magnitude:")
    print(p2 * a1 * mu1 / (2*np.pi * a2))
    print(p2 * a1 * mu1 * np.cos(arg1) * np.sqrt(1-e2 ** 2) / (2*np.pi * a2 * (1 - e2 * np.sin(arg1))))
    print("Outer Perturber TTV magnitude:")
    print(p1 * a2 * mu2 / (2*np.pi * a1))
    print(p1 * a2 * mu2 * np.cos(arg2) * np.sqrt(1-e1 ** 2) / (2*np.pi * a1 * (1 - e1 * np.sin(arg2))))
    print(mu2 * e2 * (a1/a2)**3 * p2)
    print(target, perturber)
    print("Simlulation")
    print(max(abs(TTV)))
    #print(p2 * )
    # TTV is in seconds
    #plt.plot(x, TTVp, label="Model given zero eccentricity")
    plt.plot(x, -max(abs(TTV))*np.sin(x ), label="Model given zero eccentricity")

plt.legend()
plt.ylim((-np.ceil(maxval)-1, np.ceil(maxval)+1))
plt.show()
