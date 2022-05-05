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
    else:
        return pd.DataFrame()

def fetchExpectedEpochs(planet, times):
    # compute the earliest time in the times array
    initial_t = np.min(times)
    # fetch the orbital period and initial epoch
    per, t0 = planet["pl_orbper"], planet["pl_tranmid"]
    # compute how many orbits the planet made between t0 and the start of the dataset
    ndif = (t0 - initial_t) // per
    # for all possible transits that could have been observed
    t_exp = np.arange(t0-ndif*per, np.max(times), per)
    # compute the minimum distance between the transit and the array of time datapoint
    t_dist = [np.abs(times - t).min() for t in t_exp]
    # find only the transits that are within the data
    transits = np.where(t_dist < (times[1]-times[0]))[0]
    # return the epoch times and epoch number
    return t_exp, transits

def fetchJulietPriors(df, times):
    # define the standard priors that are consistent between models
    param = ["q1_TESS", "q2_TESS", "rho", "mdilution_TESS", "mflux_TESS", "sigma_w_TESS", "GP_sigma_TESS", "GP_rho_TESS"]
    dists = ["uniform", "uniform", "loguniform", "fixed", "normal", "loguniform", "loguniform", "loguniform"]
    hypes = [[0., 1.], [0., 1.], [100., 10000.], 1.0, [0.,0.1], [0.1, 1000.], [1e-6, 1e6], [1e-3, 1e3]]
    # sort by increasing semimajor axis (or increasing period)
    df = df.sort_values(by="pl_orbsmax")
    # for each planet in the system
    for i, idx in enumerate(df.index):
        # fetch this planet
        pl = df.iloc[i]
        # if the planet is non-transiting, we don't need to fit for it
        if not pl["tran_flag"]:
            continue
        # compile the required priors
        param += "r1_p{0},r2_p{0},ecc_p{0},omega_p{0}".format(i+1).split(",")
        dists += ["uniform", "uniform", "fixed", "fixed"]
        hypes += [[0., 1.], [0., 1.], 0., 90.]
        # find the expected transits that are within the data
        t_exp, transits = fetchExpectedEpochs(pl, times)
        # for each transit, add the predicted transit time to the prior parameters
        param += list(["T_p{0}_TESS_{1}".format(i+1, t) for t in transits])
        hypes += list([[round(t_exp[t],1), .2] for t in transits])
        dists += list(["normal" for t in transits])
    # build the priors dictionary
    return juliet.utils.generate_priors(param, dists, hypes)

def _lcData_(target, lc, planetdf, loc="/raw_data/tess_data/"):
    # fetch the directory for this target
    folder = loc+target.replace(" ", "_")+"/Transit"
    # decompose the lc dataframe into time, flux, and error
    t, f, ferr = np.array(lc["time"]), np.array(lc["flux"]), np.array(lc["flux_e"])
    # if the juliet fit already exists
    if os.path.exists(folder):
        # load from the folder. prepend the syspath as juliet uses absolute filepaths
        return juliet.load(input_folder = os.getcwd()+folder)
    else:
        # ensure the output folder exists
        makeFolder(folder)
    # fetch the priors for the juliet model
    priors = fetchJulietPriors(planetdf, t)
    # Define all of the dictionaries
    times, fluxes, fluxes_error = {}, {}, {}
    # setup the dictionaries for juliet
    times["TESS"], fluxes["TESS"], fluxes_error["TESS"] = t, f, ferr
    # and return the juliet stuff. prepend the syspath as juliet uses absolute filepaths
    return juliet.load(priors = priors, t_lc = times, y_lc = fluxes, \
                       yerr_lc = fluxes_error, GP_regressors_lc = times, \
                       out_folder = os.getcwd()+folder)

def plotLC(target, lcdata, model):
    # figure size
    fig = plt.figure(figsize=(12,4))
    # Plot data
    plt.errorbar(lcdata.times_lc["TESS"], lcdata.data_lc["TESS"], yerr=lcdata.errors_lc["TESS"], \
                 fmt=".", alpha=0.1, label="TESS Data")
    # Plot the model curve
    plt.plot(lcdata.times_lc["TESS"], model, color="k", zorder=10, label="Lightcurve fit")
    # fetch the maximum deviation of the model
    lim = max(abs(model-1))
    ylim = float(lim + np.average(abs(lcdata.errors_lc["TESS"])))
    # set the limits
    plt.ylim([1-ylim, 1+ylim])
    plt.xlabel("Time (BJD)")
    plt.ylabel("Relative flux")
    plt.title("TESS Lightcurve data for {}".format(target))
    plt.legend()
    plt.show()

def tessLCData(target, planetdf, loc="/raw_data/tess_data/", stale_time=7, returnData=False):
    # lowercase the target name, and filepath prepended
    target = target.lower()
    # fetch the lightcurve data
    df = _lightcurves_(target, loc, stale_time)
    # if there wasn't any for this planetary system, or we don't want to return the data
    if (df.empty) or (not returnData):
        # return empty
        return pd.DataFrame(), pd.DataFrame()
    # fetch the juliet lightcurve data
    lcdata = _lcData_(target, df, planetdf, loc)
    # return the lc data
    return df, lcdata

def tessMidTransits(target, planetdf, loc="/raw_data/tess_data/", stale_time=7):
    # lowercase the target name, and filepath prepended
    target = target.lower()
    # fetch the Lightcurve data
    df, lcdata = tessLCData(target, planetdf, loc, stale_time, returnData=True)
    # if there wasn't any for this planetary system
    if df.empty:
        return
    # fetch the results
    results = lcdata.fit(sampler="dynesty")#, nthreads=4)
    # fetch the transits
    print(results.posteriors)
    # fetch the transit model
    model = results.lc.evaluate("TESS")
    # plot the fit
    plotLC(target, lcdata, model)
