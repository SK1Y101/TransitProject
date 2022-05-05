""" Use Juliet to fetch TESS ligthcurve data """
# python modules

from astroquery.mast import Observations
from matplotlib import pylab as plt #testing purposes
from astropy.table import join
from astropy.io import fits
import juliet, re, os
import pandas as pd
import numpy as np

# my modules
from . import *

def _mastQuery_(target):
    # start with no target
    target_name = ""
    # check if the object passed was a tess object
    match = re.match("^(tess|tic) ?(\d+)$", target)
    # if it was
    if match:
        target_name = f"{match.group(2).zfill(9)}"
    # if we have a target to search for
    if target_name:
        # fetch the observation
        obs = Observations.query_criteria(target_name=target_name)
        # return the observation
        return obs

def _downloadlc_(obs, fluxType="PDCSAP_FLUX"):
    # fetch the fits url for this object
    url = "https://mast.stsci.edu/api/v0.1/Download/file?uri={}".format(obs["dataURI"])
    # load the data
    data = fits.getdata(url)
    # ensure we have the data we want
    if fluxType not in data.columns.names:
        return
    # remove zero or nan flux locations
    mask = np.where((data[fluxType]!=0)&(~np.isnan(data[fluxType])))[0]
    # fetch the times
    time = data["TIME"][mask]
    # fetch the flux and flux error, normalised by the median
    flux = data[fluxType][mask] / np.median(data[fluxType][mask])
    ferr = data[fluxType+"_ERR"][mask] / np.median(data[fluxType][mask])
    # compute the units for each column we"re searching
    cols = [np.where(np.array(data.columns.names) == reqCol)[0][0] for reqCol in ["TIME", fluxType, fluxType+"_ERR"]]
    units = [data.columns.units[col] for col in cols]
    # return the data
    return time, flux, ferr, units

def _lightcurves_(target, loc="/raw_data/tess_data/", stale_time=7, throwError=False):
    # filename for this target
    lcfname = loc+target.replace(" ","_").lower()+"/raw_lc.csv"
    # check if we have local and non-stale data to use
    if checkForFile(lcfname, stale_time):
        # load the dataframe and return
        return loadDataFrame(lcfname)
    # fetch mast observations
    obs = _mastQuery_(target)
    # if we had any
    if obs:
        # fetch the data products given the observations
        prod = Observations.get_product_list(obs)
        # ensure we are only using science data
        prod = Observations.filter_products(prod, productType=["SCIENCE"])
        # join the two results into a single table
        res = join(obs, prod, keys="obs_id", join_type="right", table_names=["", "products"])
        # sort by obj id
        res.sort("obs_id")
        # ensure mast knows we are using tess data
        res["obs_collection"] = "TESS"
        # fetch all the lightcurves that aren"t the fast-lc
        mask = np.array([fname.endswith("lc.fits") and not fname.endswith("fast-lc.fits") for fname in res["productFilename"]])
        res = res[mask]
        # output arrays (and default bjd offset)
        times, fluxes, ferrs, unit = np.array([]), np.array([]), np.array([]), ("BJD - 2457000",)
        # for each file found
        for i, item in enumerate(res):
            # ensure we are using the tess data, and the lc, not the
            if item["obs_collection_"] == "TESS":
                # attempt to fetch the lightcurve data
                lcdata = _downloadlc_(item)
                # if it was successful
                if lcdata:
                    # unpack times, flux, flux errors, and units
                    t, f, ferr, unit = lcdata
                    # append to respective arrays
                    times = np.append(times, t)
                    fluxes = np.append(fluxes, f)
                    ferrs = np.append(ferrs, ferr)
                    units = unit
        # create dataframe from the data
        df = pd.DataFrame()
        df["time"], df["flux"], df["flux_e"] = times, fluxes, ferrs
        # convert our times from bjd - 2457000 to bjd
        bjdoffset = findFloats(unit[0])[0].split(" ")[1]
        df["time"] += int(bjdoffset)
        # ensure the folder exists
        makeFolder(lcfname)
        # save the dataframe
        saveDataFrame(df, lcfname)
        # and return it too
        return df
    # if there was no observation data
    if throwError:
        raise Exception("No LightCurve Data for {}".format(target))

def fetchStandardJulietPriors(df, t):
    # define the standard priors that are consistent between models
    param = ["q1_TESS", "q2_TESS", "rho", "mdilution_TESS", "mflux_TESS", "sigma_w_TESS", "GP_sigma_TESS", "GP_rho_TESS"]
    dists = ["uniform", "uniform", "loguniform", "fixed", "normal", "loguniform", "loguniform", "loguniform"]
    hypes = [[0., 1.], [0., 1.], [100., 10000.], 1.0, [0.,0.1], [0.1, 1000.], [1e-6, 1e6], [1e-3, 1e3]]
    # sort by increasing semimajor axis (or increasing period)
    df = df.sort_values(by="pl_orbsmax")
    # compute the earliest time in the times array
    initial_t = np.min(t)
    # for each planet in the system
    for i, idx in enumerate(df.index):
        # fetch this planet
        pl = df.iloc[i]
        # if the planet is non-transiting, we don't need to fit for it
        if not pl["tran_flag"]:
            continue
        # fetch the orbital period and initial epoch
        per, t0 = pl["pl_orbper"], pl["pl_tranmid"]
        # fetch the size of a standard deviation by considering the distance between the upper and lower errors
        persig, t0sig = 0.5*(pl["pl_orbpererr1"]-pl["pl_orbpererr2"]), 0.5*(pl["pl_tranmiderr1"]-pl["pl_tranmiderr2"])
        # compute how many orbits the planet made between t0 and the start of the dataset
        ndif = (t0 - initial_t) // per
        # rewind t0 as far as possible, and compund the sigma terms
        t0, t0sig = t0-ndif*per, t0sig*ndif
        # fetch the priors on its values
        param += "P_p{0}, t0_p{0}, r1_p{0}, r2_p{0}, ecc_p{0}, omega_p{0}".format(i+1).split(",")
        dists += ["normal", "normal", "uniform", "uniform", "fixed", "fixed"]
        hypes += [[per, persig], [t0, t0sig], [0., 1.], [0., 1.], 0., 90.]
    # build the priors dictionary
    return juliet.utils.generate_priors(param, dists, hypes)

def _lcData_(target, lc, planetdf, loc="/raw_data/tess_data/"):
    # fetch the directory for this target
    folder = loc+target.replace(" ", "_")+"/Transit"
    # decompose the lc dataframe into time, flux, and error
    t, f, ferr = np.array(lc["time"]), np.array(lc["flux"]), np.array(lc["flux_e"])
    # if the juliet fit already exists
    if os.path.exists(folder):
        # load from the folder
        return juliet.load(input_folder = folder)
    else:
        # ensure the output folder exists
        makeFolder(folder)
    # fetch the priors for the juliet model
    priors = fetchStandardJulietPriors(planetdf, t)
    # Define all of the dictionaries
    times, fluxes, fluxes_error = {}, {}, {}
    # setup the dictionaries for juliet
    times["TESS"], fluxes["TESS"], fluxes_error["TESS"] = t, f, ferr
    # and return the juliet stuff
    return juliet.load(priors = priors, t_lc = times, y_lc = fluxes, \
                       yerr_lc = fluxes_error, GP_regressors_lc = times, \
                       out_folder = folder)

def tessMidTransits(target, planetdf, loc="/raw_data/tess_data/", stale_time=7):
    # lowercase the target name, and filepath prepended
    target, loc = target.lower(), os.getcwd()+loc if loc[0]=="/" else loc
    # fetch the lightcurve data
    df = _lightcurves_(target, loc, stale_time)
    # fetch the juliet lightcurve data
    lcdata = _lcData_(target, df, planetdf, loc)
    # fetch the phases
    #P,t0 =  4.7423729 ,  2457376.68539
    #phases = juliet.utils.get_phases(df["time"], P, t0)
    # fetch the results
    results = lcdata.fit(sampler="dynesty")#, nthreads=4)
    # Extract full model:
    transit_plus_GP_model = results.lc.evaluate("TESS")

    # Deterministic part of the model (in our case transit divided by mflux):
    transit_model = results.lc.model["TESS"]["deterministic"]

    # GP part of the model:
    gp_model = np.array(results.lc.model["TESS"]["GP"])

    # Now plot. First preambles:
    import matplotlib.gridspec as gridspec
    fig = plt.figure(figsize=(12,4))
    gs = gridspec.GridSpec(1, 2, width_ratios=[2,1])
    ax1 = plt.subplot(gs[0])

    # Plot data
    ax1.errorbar(lcdata.times_lc["TESS"], lcdata.data_lc["TESS"], yerr = lcdata.errors_lc["TESS"], fmt = ".", alpha = 0.1)

    # Plot the (full, transit + GP) model:
    ax1.plot(lcdata.times_lc["TESS"], transit_plus_GP_model, color="black",zorder=10)

    ax1.set_xlim([2457000 + 1328, 2457000 + 1350])
    ax1.set_ylim([0.96, 1.04])
    ax1.set_xlabel("Time (BJD)")
    ax1.set_ylabel("Relative flux")

    #ax2 = plt.subplot(gs[1])

    # Now plot phase-folded lightcurve but with the GP part removed
    #removed_gp = np.array(lcdata.data_lc["TESS"]) - np.array(gp_model)
    #ax2.errorbar(np.array(phases), np.array(lcdata.data_lc["TESS"] - gp_model), yerr = np.array(lcdata.errors_lc["TESS"]), fmt = ".", alpha = 0.3)
    #ax2.errorbar(phases, lcdata.data_lc["TESS"], yerr = lcdata.errors_lc["TESS"], fmt = ".", alpha = 0.3)

    # Plot transit-only (divided by mflux) model:
    #idx = np.argsort(phases)
    #ax2.plot(phases[idx], transit_model[idx], color="black",zorder=10)
    #ax2.yaxis.set_major_formatter(plt.NullFormatter())
    #ax2.set_xlabel("Phases")
    #ax2.set_xlim([-0.03, 0.03])
    #ax2.set_ylim([0.96, 1.04])
    # fit transits to the lightcurve data using juliet
    #midTransitTimes = fitTransits(df["Times"])
    #plt.errorbar(t, f, yerr=ferr, fmt=".", alpha=0.2)
    #plt.xlabel("Time (BJD)")
    #plt.ylabel("Relative flux")
    #plt.xlim((2457000 + 1325, 2457000+1355))
    #plt.ylim((0.94, 1.06))
    plt.show()
