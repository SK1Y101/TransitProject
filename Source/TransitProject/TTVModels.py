''' Handles the creation and execution of different TTV Models '''
# Python modules
from multiprocessing import Pool
import matplotlib.pylab as plt
import emcee, corner
import numpy as np

# my modules
from . import *

# compute an orbital period given a semimajor axis
def pfroma(a, mu, smass=1):
    return 2*np.pi*np.sqrt(a**3 / (6.67408E-11 * 1.98855E30 * smass * (1 + mu)))

def model1(t, star, target, planets):
    ''' Model for a transiting planet exterior to 'n' coplanar, non-interacting, zero eccentricity perturbing planets. '''
    # position of the barycentre
    rb = [-a * mu * np.sin(2*np.pi*(t-t0) / pfroma(a, mu, star["mass"])) \
          for a, mu, e, i, arg, t0 in planets]
    # TTV due to the position
    TTV = (pfroma(target[0], target[1], star["mass"]) / (2*np.pi*target[0])) * sum(rb)
    # return the TTV
    return TTV

def dfToParams(df):
    # fetch the star params
    star = df.iloc[0]
    # and the params for each planet > (a, P, mu, ecc, inc, arg, t0)
    planets = Simulation.planetParams(df)
    # remove orbital period as we're calculating it in models
    planets = [np.delete(pl, 1) for pl in planets]
    # return
    return star, planets

def log_likelihood(theta, x, y, yerr, model=None):
    # evaluate model
    mod = model(x, theta)
    # sigma
    sigma2 = yerr ** 2 + mod ** 2
    # evaluate and return
    return -0.5 * np.sum((y-mod)**2 / sigma2 + np.log(sigma2))

def log_prior(theta, priors=None):
    # ensure we have enough priors
    assert len(priors) == len(theta) if isinstance(theta, (list, tuple)) else 1
    # check values are between priors
    if all(p[0]<=t<=p[1] for t,p in zip(theta, priors)):
        return 0
    return -np.inf

def log_prob(theta, x, y, yerr, priors=None, model=None):
    # evaluate priors
    lp = log_prior(theta, priors=priors)
    # if not infinite, return infinite
    if not np.isinf(lp):
        return -np.inf
    # return the log likelihood
    return lp + log_likelihood(theta, x, y, yerr, model=model)

def BIC(k, n, L):
    ''' compute the Bayesian information criterion of the fit.
        k: Number of free parameters.
        n: Number of data points.
        L: Maximum observed value of the log_likelihood function.'''
    return k*np.log(n) - 2*L

def priorsToInit(p, walker=32):
    # construct the entire parameter space choices
    inits = [np.linspace(pr[0], pr[1], walker*100) for pr in p]
    # select random elements for each walker
    return np.vstack([np.random.choice(init, walker) for init in inits]).T

def toSI(theta):
    conv = np.reshape(theta, (-1, 6))
    conv[:,0] *= 149597870700
    conv[:,5] = (conv[:,5]+2400000)*86400
    return conv

def modelToMCMC(x, y, yerr, model, star, target, perturbers=1):
    ''' Convert observed data, and TTV models to emcee samplers. '''
    # create parameter space
    priors, planets, labels = [], [], []
    # parameter space for perturbers
    for p in range(perturbers):
        # reduced mass
        mu, mup = target[2], [0, 0.25]
        # semimajor axis between 0 and 100 AU
        a, ap = target[0] / 149597870700, [0, 100]
        # eccentricity
        e, ep = target[3], [0, 1]
        # initial epoch
        t0, t0p = target[5]-2400000, [0, 60000]
        # inclination and arg are fixed for now
        i, ip, arg, argp = 0, [0, 0], 0, [0, 0]
        # add values to planet
        planets+=[a, mu, e, i, arg, t0]
        # add values to priors
        priors+=[ap, mup, ep, ip, argp, t0p]
        # add labels
        labels+="a_{0},mu_{0},e_{0},i_{0},arg_{0},t0_{0}".format(p).split(",")
    # randomly sample from the parameter space
    pos = planets + np.random.randn(32, len(planets))
    nwalk, ndim = pos.shape
    # construct the unchanging parts of the model, and ensure the planet params are Nx6 array
    mod = lambda x, theta : model(x, star, target, toSI(theta))
    # pass the model to the log_prob function
    log_prob_mod = lambda theta, x, y, yerr : log_prob(theta, x, y, yerr, priors, mod)
    # construct the sampler
    sampler = emcee.EnsembleSampler(nwalk, ndim, log_prob_mod, args=(x, y, yerr))
    # and return it
    return pos, sampler, mod, labels

def runSampler(pos, sampler, labels, samples=1000):
    #with Pool() as pool:
        # assign the pool to the sampler
    #    sampler.pool = pool
        # run the sampler
    sampler.run_mcmc(pos, samples, progress=True)
    # fetch some samples, discarding 10% and thining to 50% of the sample count
    flat_samples = sampler.get_chain(discard=int(samples*0.1), flat=True)
    # plot a corner plot
    #fig = corner.corner(flat_samples, labels=labels)
    #plt.show()

def MCMCFit(x, y):
    # experimental model
    def model(x, magnitude, period, phase):
       return magnitude * np.sin(x * 2 * np.pi / period + phase)

    # construct model
    with pm.Model() as transitModel:
        # data
        data = pm.Data("data", y)
        # parameters
        mag = pm.Uniform("mag", lower=0, upper=1000)
        per = pm.Uniform("per", lower=0, upper=1000)
        phase = pm.Uniform("phase", lower=0, upper=2*np.pi)
        # expected outcome
        mu = model(x, mag, per, phase) #mag * pm.math.sin(x * (2 * np.pi / per) + phase)
        # observation
        obs = pm.Normal("y", mu=mu, sigma=1, observed=y)

    with transitModel:
        trace = pm.sample(1000, init="adapt_diag", return_inferencedata=False, chains=4)
        az.plot_trace(trace, var_names=["mag", "per", "phase"])
        print(pm.summary(trace))
