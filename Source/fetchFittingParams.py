import pandas as pd
import numpy as np

df = pd.read_csv("TTVTestModelFittingParameters.csv")

p1 = pd.DataFrame(columns=["Model", "Method", "mu", "a", "e", "omega"])
p2 = pd.DataFrame(columns=["Model", "Method", "mu", "a", "e", "omega"])
p3 = pd.DataFrame(columns=["Model", "Method", "mu", "a", "e", "omega"])

for idx in df.index:
    thisfit = df.iloc[idx]

    name = thisfit["models"]
    params = [x.strip() for x in thisfit["solutions"].replace("[", "").replace("]", "").split(" ") if x][:-1]

    per = int(name.split(":")[1].split(",")[0].replace("perturbers", "").strip())

    model = name.split(":")[0].replace("model", "Model ")
    method = name.split(",")[1].replace("differential", "diff.").replace("evolution", "evo.").replace("annealing", "an.").replace("fitting", "")

    for i, x in enumerate(np.reshape(params, (-1, 6))):
        mu = round(float(x[1])*10**6, 3)
        a = round(float(x[0])*10**6, 3)
        e = round(float(x[2]), 3)
        omega = round(float(x[4]), 3)

        thisrow = {"Model":"" if i else model, "Method":"" if i else method, "mu":mu, "a":a, "e":e, "omega":omega}

        if per == 1:
            p1 = p1.append(thisrow, ignore_index=True)
        if per == 2:
            p2 = p2.append(thisrow, ignore_index=True)
        if per == 3:
            p3 = p3.append(thisrow, ignore_index=True)

print(p1.to_latex(index=False))
print(p2.to_latex(index=False))
print(p3.to_latex(index=False))

udf = pd.read_csv("TTVTestModelFittingUncertainties.csv")
u1 = pd.DataFrame(columns=["Model","Method","k","ln L","AIC","AICc","BIC"])
for idx in udf.index:
    thisfit = udf.iloc[idx]
    name = thisfit["models"]

    thisrow = {}

    thisrow["Model"] = name.split(":")[0].replace("model", "Model ")
    thisrow["Method"] = name.split(",")[1].replace("differential", "diff.").replace("evolution", "evo.").replace("annealing", "an.").replace("fitting", "")
    thisrow["k"] = thisfit["k"]
    thisrow["ln L"] = thisfit["ln L"]
    thisrow["AIC"] = thisfit["AIC"]
    thisrow["AICc"] = thisfit["AICc"]
    thisrow["BIC"] = thisfit["BIC"]

    u1 = u1.append(thisrow, ignore_index=True)
print(u1.to_latex(index=False))

# generate an arbitary corner plot
from TransitProject import modelling as tm
from TransitProject import Simulation as ts

import simulationPipeline as simP
from astropy import units as u
# define an arbitary start time
start_time, end_time = 0, 60000
# simulate and fetch the TTV
TTV, TTVl, TTVu, TT = simP.runMinimalSim("testSystem.csv,transiting_planet", start_time, end_time/365.25, False)
# convert parameters and return
x = TT*365.25*u.d
y = TTV
yerr = np.average((TTVl, TTVu))

best = 8
#best = np.argmin(udf["AIC"])
import ast
sol = df["solutions"].iloc[best].replace(" ", ",").replace("\n", ",").replace("[,", "[").replace(",]", "]").replace(",,", ",").replace(",,", ",")
sol = np.array(ast.literal_eval(sol))
print(sol)
print(sol.shape)
priors = df["priors"].iloc[best].replace(" ", ",").replace("\n", ",").replace("[,", "[").replace(",]", "]").replace(",,", ",").replace(",,", ",")
priors = np.array(ast.literal_eval(priors)).reshape((-1, 2))
print(priors)
print(priors.shape)
labels = np.array(ast.literal_eval(udf["labels"][best]))
print(labels)
print(labels.shape)

def dfToParams(df):
    from astropy import units as u
    # fetch the star paramsfetchPara
    star = ts.bodyParams(df)[0]
    # and the params for each planet > (a, P, mu, ecc, inc, arg, t0)
    planets = ts.planetParams(df)
    # vstack the bodies
    bodies = np.vstack([star, planets])
    # remove orbital period as we're calculating it in models
    bodies = np.delete(bodies, 1, axis=1)
    # convert reduced mass
    # used to be m / mTot, change to raw mass in Msun
    bodies[:, 1] = df["mass"].to_numpy()
    # convert from m to AU, and days to years (This ensure the parameter space search is smaller)
    bodies[:, 0] = bodies[:, 0] * (1*u.m).to(u.AU).value
    bodies[:, 5] = bodies[:, 5] * (1*u.d).to(u.yr).value
    # return > (a, mu, ecc, inc, arg, t0)
    return bodies

bodies = dfToParams(ts.fetchParams("testSystem.csv,transiting_planet")[0])
model = tm.model1

default = bodies[:2].flatten()
# use default parameters with the model
def modelHandler(x, *theta):
    # prepend the default values
    params = np.append(default, theta).reshape((-1, 6))
    # convert from parameter space
    params = tm.pSpaceToReal(params)
    # return the model with the default values included
    return model(x, *params.flatten())

modelHandler.__name__ = df["models"][best]

import emcee
walkers = max(4*len(sol), 200)
pos = sol + 1e-4 * np.random.randn(walkers, len(sol))
nwalk, ndim = pos.shape
sampler = emcee.EnsembleSampler(nwalk, ndim, tm.log_prob, args=(x, y, yerr, priors, modelHandler))
# and run to find parameters
sampler.run_mcmc(pos, 20000, progress=True)


flat_samples = sampler.get_chain(discard=100, flat=True)
best, besterror = [], []
# for each parameter
for idx, val in enumerate(sol):
    mcmc = np.percentile(flat_samples[:, idx], [16, 50, 84])
    q = np.diff(mcmc)
    best.append(mcmc[1])
    besterror.append([[q[0], q[1]]])

import corner
import matplotlib.pylab as plt
fig = corner.corner(flat_samples[:, :-1], quantiles=(0.16, 0.5, 0.84), show_titles=True, \
                    labels=labels, title_fmt=".2E")
plt.savefig("TTVTestModelFittingCorner.pdf", bbox_inches="tight")
plt.savefig("TTVTestModelFittingCorner_transparent.pdf", bbox_inches="tight", transparent=True)
plt.show()

# combine default system params and found solution
sysParams = np.append(bodies[:2], sol[:-1]).reshape((-1, 6))
# construct dataframe for system
names = ["Example System", "Transiting Planet"]
df = pd.DataFrame(columns=["name", "mass", "sma", "ecc", "inc", "arg"])
df["name"] = list(names) + [f"perturber {p+1}" for p in range(sysParams.shape[0]-2)]
df["mass"] = sysParams[:, 1]
df["sma"] = sysParams[:, 0]
df["ecc"] = sysParams[:, 2]
df["inc"] = sysParams[:, 3]
df["arg"] = sysParams[:, 4]
df["mean"] = np.zeros(sysParams.shape[0])
# to sim array
simArray = ts.constructSimArray(df, useerror=False)
# convert to simulation
sim = ts._arrayToSim_(simArray[0], list(df.columns[1:]))

# plot the simulation
import rebound
from matplotlib.lines import Line2D
fig, ax_main, ax_top, ax_right = rebound.OrbitPlot(sim, slices=0.5, xlim=[-2.,2], ylim=[-2.,2], color=True, unitlabel="[AU]")
fig.suptitle(df.iloc[0]["name"].capitalize())
# create the legend
cols = [(1.,0.,0.),(0.,0.75,0.75),(0.75,0.,0.75),(0.75, 0.75, 0,),(0., 0., 0.),(0., 0., 1.),(0., 0.5, 0.)]
leg_elem = [Line2D([0], [0], color=cols[i], lw=1.5, label="{} Orbit".format(df.iloc[i+1]["name"].capitalize())) \
            for i, item in enumerate(df.iloc[1:].index)]
# apply the legend
plt.legend(handles = leg_elem, framealpha=0, bbox_to_anchor=(0, 1.5), loc="upper left")
# and save
fig.tight_layout()
fig.savefig("{}_BestFitLayout.pdf".format(df.iloc[0]["name"].replace("_","")), transparent=False, bbox_inches='tight')
fig.savefig("{}_BestFitLayoutTransparent.pdf".format(df.iloc[0]["name"].replace("_","")), transparent=True, bbox_inches='tight')
plt.show()
