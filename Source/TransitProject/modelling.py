# python modules
from multiprocessing import Pool
import matplotlib.pylab as plt
import astropy.constants as c
import scipy.optimize as opt
import astropy.units as u
from scipy import stats
from tqdm import tqdm
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
    # compute the TTV offset
    TTV = (P[1] / (2*np.pi * bodies[1][0])) * sum(r_b)
    # return
    return TTV / u.s

# 'n' interior, coplanar, non-interacting, non-zero-eccentricity, perturbing planets, including the transiting planets eccentricity
def model3(t, *bodies):
    # fetch period and 2d bodies array
    bodies, P = _extraModelParam_(*bodies)
    # fetch true anomalies
    V = _trueFromMean_(t, P, bodies)
    # compute distance at current point in orbit
    R = [a * (1-e**2) / (1 + e*np.cos(v*u.radian)) for (a,mu,e,_,arg,t0),v in zip(bodies, V)]
    # compute the barycentre shift
    r_b = [-r * mu * np.sin(v*u.radian+arg*u.radian) for (a,mu,e,_,arg,t0),p,v,r in zip(bodies[2:], P[2:], V[2:], R[2:])]
    # compute the velocity offset due to eccentricity
    eccoff = np.sqrt(1-bodies[1][2]**2) / np.sqrt(1 + 2*bodies[1][2]*np.cos(V[1]*u.radian+bodies[1][4]*u.radian) + bodies[1][2]**2)
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

def _intersectingHillSphere_(bodies):
    # mass, semi-major axis, and eccentricity
    M = bodies[:, 1]
    a = bodies[:, 0]
    e = bodies[:, 2]
    # inner and outer regions of stability
    inner = a[1:] * (1-e[1:]) * (1 - (M[1:] / (3*M[0]))**(1/3))
    outer = a[1:] * (1+e[1:]) * (1 + (M[1:] / (3*M[0]))**(1/3))
    # if any planet lies within the hill sphere of another
    intersection = [np.logical_and(np.delete(inner, i) <= body[0], body[0] <= np.delete(outer, i)) for i, body in enumerate(bodies[1:])]
    # return as an array
    return np.array(intersection)

''' < Statistical programming to make the whole fitting thing work > '''

def log_like(theta, x, y, yerr, model, additionalCheck=lambda x: False):
    # params and logarithm of the underestimation
    params, log_f = theta[:-1], theta[-1]
    # expected value
    exp = model(x, *params)
    # variance
    sigma2 = yerr ** 2 + exp ** 2 * np.exp(2 * log_f)
    # if we have an additional parameter check that is not a set of boundaries
    if additionalCheck(params):
        # (The function must return true to fail these parameters)
        return -np.inf
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

def chooseRange(arr, choices=1, resolution=1000):
    # ensure we have a numpy array
    arr = np.array(arr)
    # parameter space from upper and lower bounds
    pSpace = np.linspace(arr[:,0], arr[:,1], resolution)
    # ensure it is shaped correctly
    pSpace = np.reshape(pSpace, pSpace.shape[:2])
    # select random indices of the parameter space
    ind = np.random.rand(*pSpace.shape).argsort(axis=0)[:choices]
    # select the elements from the parameter space
    return np.take_along_axis(pSpace, ind, axis=0).T

def basin(func, p0, x, y, yerr, model, priors, niter=1000):
    ''' Define the basin hopping method for fitting a model to data. '''
    with tqdm(total=niter, smoothing=0, desc="Basin hopping") as bar:
        def update(*a):
            bar.update(1)
        return opt.basinhopping(func, p0, minimizer_kwargs={"args":(x, y, yerr, model), "bounds":priors}, callback=update, niter=niter)

def difevo(func, p0, x, y, yerr, model, priors, niter=10000):
    ''' Define the differential evolution method for fitting a model to data. '''
    with tqdm(total=niter, smoothing=0, desc="Differential Evolution") as bar:
        def update(*a, convergence=0):
            bar.set_postfix_str(f"convergence: {100*min(1,convergence):0.2f}%", refresh=False)
            bar.update(1)
        return opt.differential_evolution(func, priors, x0=p0, args=(x, y, yerr, model), callback=update, polish=True, maxiter=niter)

def dual(func, p0, x, y, yerr, model, priors, niter=10000):
    ''' Define the dual annealing method for fitting a model to data. '''
    with tqdm(total=niter, smoothing=0, desc="Dual annealing") as bar:
        def update(*a):
            bar.update(1)
        return opt.dual_annealing(func, priors, x0=p0, args=(x, y, yerr, model), callback=update, maxiter=niter)

def parameterSearch(x, y, yerr, priors, model, initial=None, method="dif", additionalCheck=lambda x: False):
    # negative log likelihood
    def nll(*args):
        return -log_like(*args, additionalCheck=additionalCheck)
    # if we were not supplied initial parameters
    if not isinstance(initial, (np.ndarray, list, tuple)):
        # create parameter guesses
        initial = chooseRange(priors)
    # add the underestimation factor too
    initial = np.array(list(initial)+[0])
    priors = np.array(list(priors)+[[-100, 100]])
    # choose the minimiser method
    if method == "lsq":
        # Least Squares Regression (Local)
        solution = opt.least_squares(nll, initial, args=(x, y, yerr, model))
    elif method == "basin":
        # Basin Hopping (Global)
        solution = basin(nll, initial, x, y, yerr, model, priors)
    elif method == "dif":
        # Differential Evolution (Global)
        solution = difevo(nll, initial, x, y, yerr, model, priors)
    elif method == "dual":
        # Dual Annealing method (Global)
        solution = dual(nll, initial, x, y, yerr, model, priors)
    elif method == "shgo":
        # Simplical Homology Global Optimisation (Global)
        solution = opt.shgo(nll, priors, args=(x, y, yerr, model))
    else:
        # Standard minimiser (Local)
        solution = opt.minimize(nll, initial, args=(x, y, yerr, model), bounds=priors)
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

''' < Graphical Outputs > '''

def pSpaceToReal(theta):
    val = theta.reshape((-1, 6))
    # convert from AU to m, and years to days
    val[:, 0] = val[:, 0] * (1*u.AU).to(u.m).value
    val[:, 5] = val[:, 5] * (1*u.yr).to(u.d).value
    return val

def plotModels(x, y, yerr, periods, models, params, xlim=None, xlimz=(0,2), fname=None, showCombined=False):
    ''' Plot the data and predicted output of given models. Will also show each model output on seperate subplots.
        x, y, yerr:   The x, y, and y-error of the data.
        periods:      The orbital period (in days) of all bodies
        models:       A list of models to compare
        params:       An array of model parameters to use for the plotting.
        xlim:         The x-limit of the combined model/data view
        xlimz:        The x-limit of the subplots for each model individually.
        fname:        The filename to save the final plot as. Leaving this as None will prevent the script from saving the output.'''
    # fetch size
    m = len(models)
    plotnum = m + int(showCombined)
    # shift plot to zero
    x -= min(x)
    # fetch transiting period
    p = periods[1]
    # convert x to transit numbers
    x = np.around(x/p)
    # remove datapoints outside of the limits if xlim is given
    if xlim != None:
        keep = np.where(np.logical_and(x >= np.floor(xlim[0]), x <= np.ceil(xlim[1])))[0]
        x = x[keep]
        y = y[keep]
        if isinstance(yerr, np.ndarray):
            if yerr.shape[0] == 2:
                yerr = [err[keep] for err in yerr]
            else:
                yerr = yerr[keep]
    # intermediary x, 100 values per orbit of the smallest period
    x1 = np.arange(min(x), max(x), 0.01 * min(periods[1:]) / p)
    # shape Nx6 for a single parameter array of N planets, or MxNx6 for M different parameter spaces
    pshape = np.array(params).shape
    if pshape[0] != len(models):
        # reshape to MxNxP
        params = np.reshape(params, (-1, pshape[-2], pshape[-1]))
        # if M == 1
        if params.shape[0] == 1:
            # repeat the parameters for each model
            params = params.repeat(len(models), axis=0)
        # ensure M is not mismatched
        assert params.shape[0] == len(models)
    # figure
    plt.figure(figsize=(10, (5/3) * plotnum))
    # plot everything
    plt.style.use('seaborn-colorblind')
    # show the first plot if desired
    if showCombined:
        plt.subplot(plotnum, 1, 1)
        for model, param in zip(models, params):
            param = param if param.size % 6 == 0 else param[:-1]
            # convert transit number to time for the model to run
            plt.plot(x1, model(x1*p, *param), label="Analytical TTV, {}".format(model.__name__.capitalize()))
        plt.errorbar(x, y, yerr=yerr, fmt="o", label="Simulated TTV", c="k")
        plt.title("Comparison between Numerical and Analytical TTV methods", fontsize="xx-large")
        plt.ylabel("TTV [s]")
        plt.xlabel("Transit number")
        # move the x axis up slightly
        ax = plt.gca()
        ax.xaxis.labelpad = -2
        if xlim:
            plt.xlim(xlim)
    # remove datapoints outside of the limits if xlim is given
    if xlimz != None:
        keep = np.where(np.logical_and(x >= np.floor(xlimz[0]), x <= np.ceil(xlimz[1])))[0]
        x = x[keep]
        y = y[keep]
        if isinstance(yerr, np.ndarray):
            if yerr.shape[0] == 2:
                yerr = [err[keep] for err in yerr]
            else:
                yerr = yerr[keep]
        x1 = np.arange(min(x), max(x), 0.01 * min(periods[1:]) / p)
    # for each model-parameter combination
    for i, modpar in enumerate(zip(models, params)):
        # decompose
        model, param = modpar
        param = param if param.size % 6 == 0 else param[:-1]
        plt.subplot(plotnum, 1, i+1+int(showCombined))
        plt.style.use('seaborn-colorblind')
        # offset colour
        for j in range(i):
            plt.plot([-10],[0])
            # convert transit number to time for the model to run
        plt.plot(x1, model(x1*p, *param), label="Analytical TTV, {}".format(model.__name__.capitalize()))
        plt.errorbar(x, y, yerr=yerr, fmt="o", label="Simulated TTV", c="k")
        plt.title("Zoomed Comparison for {}".format(model.__name__.capitalize()), fontsize="xx-large")
        plt.ylabel("TTV [s]")
        plt.xlabel("Transit number")
        # move the x axis up slightly
        ax = plt.gca()
        ax.xaxis.labelpad = -2
        if xlimz:
            plt.xlim(xlimz)
    plt.tight_layout()
    # if we were given a filename, save the plot
    if fname:
        plt.savefig(f"{fname}.pdf", transparent=False, bbox_inches='tight')
        plt.savefig(f"{fname}_transparent.pdf", transparent=True, bbox_inches='tight')
    plt.show()

''' < Model fitting and selection > '''

def BIC(k, n, L):
    ''' Compute the Bayesian information criterion of the fit.
        k: Number of free parameters.
        n: Number of data points.
        L: Maximum observed value of the log_likelihood function.'''
    return k*np.log(n) - 2*L

def AIC(k, L):
    ''' Compute the Akaike information criterion of the fit.
        k: Number of free parameters.
        L: Maximum observed value of the log_likelihood function.'''
    return 2*k - 2*L

def HQC(k, n, L):
    ''' Compute the Hannan-Quinn information criterion of the fit.
        k: Number of free parameters.
        n: Number of data points.
        L: Maximum observed value of the log_likelihood function.'''
    return 2*k*np.log(np.log(n)) - 2*L

def setup_model(model, system, perturbers=1):
    ''' Setup a model with the unchanging system parameters and any number of perturbing planets.
        model: the TTV model to setup.
        system: the system parameters for objects that are not to be fit for (ie: transiting planet & star).
        perturbers: the number of new perturbing planets to introduce to the system.'''
    #define the unchanging parameters
    initial, priors, labels = [], [], []
    # estimate the radius of the star, to define the innermost allowed orbit
    r_star = max(system[0][1]**0.57, system[0][1]**0.8) * (695700000 / 149597870700)
    # create the parameter space for changing values
    for p in range(perturbers):
        # setup the standard priors
        a  = [r_star * 10, system[1][0] / 149597870700] #AU
        mu = [1e-6, 0.05]
        e  = [0, 0.5]
        t0 = [0, 200] #BJD [2400000, 2500000] in years
        i  = [-1e-5, 1e-5]
        arg= [-2*np.pi, 2*np.pi]
        # and combine
        priors += [a, mu, e, i, arg, t0]
        # generate randomised initial parameters
        init = chooseRange([a, mu, e, i, arg, t0])
        # start with zero eccentricity
        init[2] = 0
        # flatten, convert to list, and append
        initial += list(init.flatten())
        # add labels to the output
        labels += [f"a_{p}", f"mu_{p}", f"e_{p}", f"i_{p}", f"ω_{p}", f"t0_{p}"]
    # define the default parameters (star & target planet)
    default = system[:2].flatten()
    # use default parameters with the model
    def modelHandler(x, *theta):
        # append the default values
        theta = np.append(default, theta).reshape((-1, 6))
        # convert from AU to m, and years to days
        theta[:, 0] = theta[:, 0] * (1*u.AU).to(u.m).value
        theta[:, 5] = theta[:, 5] * (1*u.yr).to(u.d).value
        # return the model with the default values included
        return model(x, *theta.flatten())
    # return the found values and model function
    return priors, initial, labels, modelHandler

def randomError(x, y, yerr):
    ''' Determine the likelihood that the observed TTV are due to an error of one kind or another.
        x, y, yerr: the observed data. '''
    # ephemerides
    def model(x, *theta):
        m, b = theta
        return x*(m/u.d) + b
    plt.show()
    priors = [[-1000, 1000], [-1000, 1000]]
    sol = parameterSearch(x, y, yerr, priors, model, initial=[0,0])
    return sol, model

''' < This is the main bit that is ran, where all the power comes > '''

def free_params(x, model, initial):
    ''' determine, systematically, the number of free parameters a model has.
        x:      Some data to provide to the model
        model:  The model to determine the free parameters of.
        initial: A set of initial parameters to use with the model.
    '''
    # determine the base values of the model with all parameters given
    base = model(x, initial)
    k = 0
    # for each of the parameters given
    for i, val in enumerate(initial):
        # construct a new parameter set with this one missing
        param = np.copy(initial)
        param[i] = 0
        # increment the free parameters if changing this value gives a different result
        if any(base != model(x, param)):
            k += 1
    return k

def optimiser(x, y, yerr, models, system, methods=["dif", "dual"], p=[1, 2, 3]):
    ''' Run the optimisation pipeline on any number of models and perturbers.
        x, y, yerr: The observed data.
        models:     A list of models to fit data to.
        system:     A list of system parameters to determine starting conditions from.
        methods:    A list of optimisation methods to be used in finiding the global solution.
        p:          A list, or range of the number of perturbers to introduce for each model.
    '''
    # method names
    m2name = {"dif":"differential evolution", "dual":"dual annealing", "lsq":"least squares regression", "basin":"basin hopping", "shgo":"simplicial homology global optimization"}
    defmethod = "bounded limited memory Broyden–Fletcher–Goldfarb–Shanno"
    # output dictionary
    solution_dict = {"models":[], "priors":[], "labels":[], "freeParams":[], "solutions":[]}
    # for each model/perturber combination
    for model in models:
        print(f"\nFitting model {model.__name__}")
        for per in p:
            for method in methods:
                # setup the model
                priors, initial, labels, modelHandler = setup_model(model, system, per)
                # determine the number of true free-parameters, includeing the underestimation factor too
                freeParams = free_params(x, modelHandler, initial) + 1
                # determine this solution, ensuring parameters that intersect hill spheres are dissallowed
                solution = parameterSearch(x, y, yerr, priors, modelHandler, initial=initial, method=method)#, additionalCheck=_intersectingHillSphere_)
                # also carry over the name
                modelHandler.__name__ = f"{model.__name__}: {per} perturbers, {m2name[method] if method in m2name else defmethod} fitting"
                # append to the dictionary
                solution_dict["models"].append(modelHandler)
                solution_dict["priors"].append(priors)
                solution_dict["labels"].append(labels)
                solution_dict["freeParams"].append(freeParams)
                solution_dict["solutions"].append(solution)
    return solution_dict

def determineUncertainties(x, y, yerr, solution):
    # output dictionary
    output = {}
    # for each setup
    for idx, model in enumerate(solution["models"]):
        priors = solution["priors"][idx]
        labels = solution["labels"][idx]
        freeParams = solution["freeParams"][idx]
        soln = solution["solutions"][idx]

        # run MCMC analysis
        sampler = parameterMCMC(x, y, yerr, soln, priors=priors, model=model, samples=10000)

        # fetch the samples
        flat_samples = sampler.get_chain(discard=100, flat=True)

        best, besterror = [], []

        # for each parameter
        for idx, val in enumerate(solution)[:-1]:
            mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
            q = np.diff(mcmc)
            best.append(mcmc[1])
            besterror.append([[q[0], q[1]]])

        print(model.__name__)
        print(best)
        from time import sleep as wait
        wait(10000000)
