if __name__ == "__main__":
    # fetch the program arguments
    args = fetchArgs()

    # fetch the system data
    systemdf = ts.fetchParams(args.planet)[0]

    # fetch the parameters for each body in the standard format
    bodies = dfToParams(systemdf)

    # if we passed a real planet:
    if ".csv" not in args.planet:
        # fetch the planetary data
        transitdf = fetchMidTransitTimes(args.planet, args.dataSource)
        # convert to useful units
        x = tp.DatetoHJD(transitdf["date"]).to_numpy() * u.d
        y = (transitdf["oc"].to_numpy() * u.min).to(u.s).value
        yerr = (0.5*(transitdf["ocel"]+transitdf["oceu"]).to_numpy() * u.min).to(u.s).value
    else:
        # compute TTV from simulation instead
        x, y, yerr = simulateMidTransitTimes(args.planet)

    # define the models used
    models = [tm.model1]#, tm.model2, tm.model3]
    # defint the number of free parameters introduced per body
    freeParams = [3, 5, 5]

    # plot the initial states of the models used
    xlim, xlimz = [0.5, 15.5], [5.5, 10.5]
    P = tm._extraModelParam_(tm.pSpaceToReal(bodies))[1].to(u.d)
    #tm.plotModels(x, y, yerr, P, models, bodies, xlim=xlim, xlimz=xlimz, fname="TTVTestModelComparison")

    # determine the optimal soltion of the fit of 'n' models with 'p' parameters and combined optimisation methods.
    solutions = tm.optimiser(x, y, yerr, models, bodies, p=[1])
    # save the solutions so they can be refered too
    solDf = pd.DataFrame.from_dict(solutions)
    solDf["models"] = [model.__name__ for model in solDf["models"]]
    tp.saveDataFrame(solDf, "TTVTestModelFittingParameters")
    # and graph them too
    tm.plotModels(x, y, yerr, P, solutions["models"], solutions["solutions"], xlim=xlim, xlimz=xlimz, fname="TTVTestModelFitting")

    # fetch the MCMC values for the models, as well as evaluating their information criterion
    MCMCVal = tm.determineUncertainties(x, y, yerr, solutions)

    print(MCMCVal)
    import time
    time.sleep(100000)

    sampleParams = sampler.chain[:, 100:, :].reshape((-1, len(solution)))

    '''# plot the errors
    import matplotlib.ticker as mtick
    plt.figure(figsize=(20, 10))
    # fetch values
    r = best/initial
    idx = np.isinf(r)|np.isnan(r)
    yerrval = np.array([b[~idx] for b in beste.T])/initial[~idx]
    # plot the zero points
    plt.plot([-0.5, len(initial)+0.5], [1, 1], "--", c="k", alpha=0.3)
    # divide up the area based on body type
    plt.text(0, max(r[~idx]+yerrval[1])/0.99, "Central Star", va="top", fontsize="x-large")
    plt.plot([5.5, 5.5], [0, 10], "--", c="k", alpha=0.2)
    plt.text(6, max(r[~idx]+yerrval[1])/0.99, "Transiting Planet", va="top", fontsize="x-large")
    plt.plot([11.5, 11.5], [0, 10], "--", c="k", alpha=0.2)
    plt.text(12, max(r[~idx]+yerrval[1])/0.99, "Perturbing Planets", va="top", fontsize="x-large")
    # limits
    plt.xlim([-0.5, len(initial)+0.5])
    plt.ylim([min(r[~idx]-yerrval[0])*0.98, max(r[~idx]+yerrval[1])/0.98])
    # plot the relative amount
    plt.errorbar(np.arange(len(best))[~idx], best[~idx]/initial[~idx], \
                 yerr=yerrval, fmt="o", label="Parameters with initial values")
    # plot anything missing
    xval = list(np.arange(len(best))[idx])+[len(tm.model2best)]
    plt.errorbar(xval, np.ones(len(xval)), fmt="o", label="Parameters without initial values")
    # legend
    plt.legend(loc="lower left")
    # labels
    plt.ylabel("Ratio of best fit parameters to true parameters", fontsize="x-large")
    plt.xlabel("Model Parameters", fontsize="x-large")
    plt.title("Ratio of best fit parameters to true parameters", fontsize="xx-large")
    # axis ticks
    plt.gca().yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0))
    plt.xticks(range(len(lab)), lab, size='small')
    #plt.savefig("TTVBestFitParameters.pdf", transparent=False, bbox_inches='tight')
    #plt.savefig("TTVBestFitParameters_transparent.pdf", transparent=True, bbox_inches='tight')
    plt.show()'''

    # fit data to models
    #tm.MCMCFit(dateAsFloat, transitdf["oc"])
    '''m1pos, m1samp, m1mod, m1labels = tm.modelToMCMC(x, y, yerr, tm.model1, star, planets[0], 2)
    tm.runSampler(m1pos, m1samp, m1labels)

    # compute BIC
    k = m1pos.shape[1]
    n = len(x)

    # fetch best fit params
    theta = []
    flat_samples = m1samp.get_chain(discard=100, thin=15, flat=True)
    for i in range(k):
        mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
        theta.append(mcmc[1])
    L = tm.log_likelihood(tuple(theta), x, y, yerr, model=m1mod)
    print(planets[1])
    print(tm.toSI(theta))
    print(tm.BIC(k, n, L))

    plt.plot(x1, m1mod(x1, tm.fromSI(planets[1:])))
    plt.plot(x1, m1mod(x1, theta))
    plt.errorbar(x, y, yerr=yerr, fmt=".")
    plt.show()'''

    #print(m1pos, m1samp)

    '''# fit a sine function to the data as a test
    dateAsFloat = tp.DatetoHJD(transitdf["date"])
    x, y = np.array(dateAsFloat), transitdf["oc"].to_numpy()
    # average errors so they are symmetric
    yerr = 0.5*(transitdf["ocel"].to_numpy() + transitdf["oceu"].to_numpy())
     fit x, y, yerr to a model.

    import corner
    print(systemdf.iloc[2]["per"], systemdf.iloc[2]["per_e1"])
    samples = np.vstack([trace[k] for k in ["mag", "per", "phase"]]).T
    corner.corner(samples, labels=["mag", "per", "phase"], show_titles=True)


    print("Model estimates")
    params = []
    for i, item in enumerate(["Magnitude", "Period", "Phase"]):
        mcmc = np.percentile(samples[:, i], [16, 50, 84])
        q = np.diff(mcmc)
        print("-{:12}: {:0.3f} ± {:0.3f}".format(item, mcmc[1], q[0], q[1]))
        params.append(mcmc[1])

    fitfunc2 = lambda x: model(x, *params)
    fitfunc2.__name__ = "MCMC fit"'''

    # plot the data with the fitting residuals
    #plotMidtransits(transitdf, [fitfunc2], args.addSim)
