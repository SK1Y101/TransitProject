# python modules
from multiprocessing import Pool
import matplotlib.pylab as plt
import astropy.constants as c
import scipy.optimize as opt
import astropy.units as u
import numpy as np
import emcee

''' < TTV Models for the fitting step > '''

def _extraModelParam_(*bodies):
    # group bodies by param length (a, mu, ecc, inc, arg, t0)
    bodies = np.reshape(bodies, (-1, 6))
    # compute periods
    P = 2*np.pi*np.sqrt((bodies[:,0]*u.m)**3 / (c.G * c.M_sun * bodies[0][1] * (1 + bodies[:, 1])))
    # return the period and body array
    return bodies, P

def _trueFromMean_(t, P, bodies):
    # compute mean anomaly
    M = [2*np.pi*(t-t0*u.d) / p.to(u.d) for (_,_,_,_,_,t0),p in zip(bodies, P)]
    # fetch eccentricity
    e = bodies[:,2]
    # fourier expansion for v
    v = [Mp + (2*ep-0.25*ep**3)*np.sin(Mp*u.radian) + 1.25*ep**2*np.sin(3*Mp*u.radian)+(13/12)*ep**3*np.sin(3*Mp*u.radian) for Mp,ep in zip(M, e)]
    # true anomaly is always in the range of 0 to 2pi
    return [v%(2*np.pi) for v in v]

# 'n' interior, coplanar, non-interacting, zero-eccentricity, perturbing planets
def model1(t, *bodies):
    # fetch period and 2d bodies array
    bodies, P = _extraModelParam_(*bodies)
    # compute the barycentre shift
    r_b = [-a * mu * np.sin(2*np.pi*(t-t0*u.d)*u.radian / p.to(u.d)) for (a,mu,_,_,_,t0),p in zip(bodies[2:], P[2:])]
    # compute the TTV offset
    TTV = (P[1] / (2*np.pi * bodies[1][0])) * sum(r_b)
    # return
    return TTV / u.s

# 'n' interior, coplanar, non-interacting, non-zero-eccentricity, perturbing planets
def model2(t, *bodies):
    # fetch period and 2d bodies array
    bodies, P = _extraModelParam_(*bodies)
    # fetch true anomalies
    V = _trueFromMean_(t, P, bodies)
    # compute distance at current point in orbit
    R = [a * (1-e**2) / (1 + e*np.cos(v*u.radian)) for (a,mu,e,_,arg,t0),v in zip(bodies, V)]
    # compute the barycentre shift
    r_b = [-r * mu * np.sin(v*u.radian+arg*u.radian) for (a,mu,e,_,arg,t0),p,v,r in zip(bodies[2:], P[2:], V[2:], R[2:])]
    # compute the velocity offset due to eccentricity
    eccoff = np.sqrt(1-bodies[1][2]**2) / np.sqrt(2 + 2*bodies[1][2]*np.cos(V[1]*u.radian+bodies[1][4]*u.radian) - (1-bodies[1][2]**2))
    # compute the TTV offset
    TTV = (P[1] / (2*np.pi * bodies[1][0])) * eccoff * sum(r_b)
    # return
    return TTV / u.s

def modeln(t, *bodies):
    # fetch period and 2d bodies array
    bodies, P = _extraModelParam_(*bodies)
    # and true anomalies
    v = _trueFromMean_(t, P, bodies)
    M = (t-bodies[1][5]*u.d) / P[1].to(u.d)
    # TTV
    TTV = bodies[2][1]*(1-bodies[2][2]**-1.5)/(2*np.pi*bodies[0][1])*(P[1]**2 / P[2])*(v[1]+bodies[2][2]*np.sin(v[1]*u.radian))
    return TTV / u.s

''' < TTV model fitting and graphing methods > '''

def plotModels(x, y, yerr, models, bodies, xlim=None, xlimz=(0,2), fname="TTVModel"):
    # fetch size
    m = len(models)
    plotnum = 1+m if m>1 else 1
    # fetch transiting body period
    p = _extraModelParam_(bodies)[1][1].to(u.d)
    # convert x to transit numbers
    x = np.around(x/p)
    # intermediary x
    x1 = np.linspace(min(x), max(x), 100000)
    # figure
    plt.figure(figsize=(20, (10/3) * plotnum if plotnum > 1 else 10))
    # plot everything
    plt.style.use('seaborn-colorblind')
    plt.subplot(plotnum, 1, 1)
    for model in models:
        # convert transit number to time for the model to run
        plt.plot(x1, model(x1*p, *bodies), label="Analytical TTV, {}".format(model.__name__.capitalize()))
    plt.errorbar(x, y, yerr=yerr, fmt="o", label="Simulated TTV", c="k")
    plt.legend()
    plt.title("Comparison between Numerical and Analytical TTV methods", fontsize="xx-large")
    plt.ylabel("TTV [s]")
    plt.xlabel("Transit number")
    if xlim:
        plt.xlim(xlim)
    # for each model
    if plotnum>1:
        for i, model in enumerate(models):
            plt.subplot(plotnum, 1, i+2)
            plt.style.use('seaborn-colorblind')
            # offset colour
            for j in range(i):
                plt.plot([-10],[0])
                # convert transit number to time for the model to run
            plt.plot(x1, model(x1*p, *bodies), label="Analytical TTV, {}".format(model.__name__.capitalize()))
            plt.errorbar(x, y, yerr=yerr, fmt="o", label="Simulated TTV", c="k")
            plt.legend()
            plt.title("Zoomed Comparison for {}".format(model.__name__.capitalize()), fontsize="xx-large")
            plt.ylabel("TTV [s]")
            plt.xlabel("Transit number")
            if xlimz:
                plt.xlim(xlimz)
    plt.tight_layout()
    plt.savefig(f"{fname}.png", transparent=False, bbox_inches='tight')
    plt.savefig(f"{fname}_transparent.png", transparent=True, bbox_inches='tight')
    plt.show()

''' < Statistical programming to make the whole fitting thing work > '''

def log_like(theta, x, y, yerr, model):
    # params and logarithm of the underestimation
    params, log_f = theta[:-1], theta[-1]
    # expected value
    exp = model(x, *params)
    # variance
    sigma2 = yerr ** 2 + exp ** 2 * np.exp(2 * log_f)
    # return likelihood
    return -0.5 * np.sum((y - exp) ** 2 / sigma2 + np.log(sigma2))

def log_prior(theta, priors):
    # determine whether values are within prior bound
    if all(p[0] <= t <= p[1] for t, p in zip(theta, priors)):
        # return zero if they are (e**0 == 1)
        return 0.0
    # return minus infinity if they aren't (e**-inf == 0)
    return -np.inf

def log_prob(theta, x, y, yerr, priors, model):
    # fetch the log prior
    lp = log_prior(theta, priors)
    # if the log prior is infinite
    if not np.isfinite(lp):
        return -np.inf
    # return log likelihood
    return lp + log_like(theta, x, y, yerr, model)

def BIC(k, n, L):
    ''' compute the Bayesian information criterion of the fit.
        k: Number of free parameters.
        n: Number of data points.
        L: Maximum observed value of the log_likelihood function.'''
    return k*np.log(n) - 2*L

def chooseRange(arr, choices=1, resolution=1000):
    # parameter space from upper and lower bounds
    pSpace = np.linspace(arr[:,0], arr[:,1], resolution)
    # ensure it is shaped correctly
    pSpace = np.reshape(pSpace, pSpace.shape[:2])
    # select random indices of the parameter space
    ind = np.random.rand(*pSpace.shape).argsort(axis=0)[:choices]
    # select the elements from the parameter space
    return np.take_along_axis(pSpace, ind, axis=0).T

def parameterSearch(x, y, yerr, priors, model, initial=None):
    # negative log likelihood
    def nll(*args):
        return -log_like(*args)
    # if we were not supplied initial parameters
    if not isinstance(initial, (np.ndarray, list, tuple)):
        # create parameter guesses
        initial = chooseRange(np.array(priors))
    # add the underestimation factor too
    initial = np.array(list(initial)+[0])
    priors = np.array(list(priors)+[[None, None]])
    # fetch the solution
    solution = opt.minimize(nll, initial, args=(x, y, yerr, model), bounds=priors)
    #solution = opt.basinhopping(nll, initial, minimizer_kwargs={"args":(x, y, yerr, model)})
    # return the parameters (including the underestimation factor)
    return solution.x

def parameterMCMC(x, y, yerr, sol, priors, model, samples=5000, walkers=32, poolWorkers=None):
    # ensure we have at least 4 walkers per parameter
    walkers = max(4*len(sol), walkers)
    # create walkers spread around the solution
    pos = sol + 1e-4 * np.random.randn(walkers, len(sol))
    # fetch the shape of parameter-walker space
    nwalk, ndim = pos.shape
    if poolWorkers:
        # start the multiprocessing pool
        with Pool(poolWorkers) as pool:
            # construct the emcee sampler
            sampler = emcee.EnsembleSampler(nwalk, ndim, log_prob, args=(x, y, yerr, priors, model), pool=pool)
            # and run to find parameters
            sampler.run_mcmc(pos, samples, progress=True)
    else:
        # construct the emcee sampler
        sampler = emcee.EnsembleSampler(nwalk, ndim, log_prob, args=(x, y, yerr, priors, model))
        # and run to find parameters
        sampler.run_mcmc(pos, samples, progress=True)
    # return the sampler
    return sampler
