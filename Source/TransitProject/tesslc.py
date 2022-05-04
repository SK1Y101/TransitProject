''' Use Juliet to fetch TESS ligthcurve data '''
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

def findFloats(txt):
    return re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", txt)

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
    # compute the units for each column we're searching
    cols = [np.where(np.array(data.columns.names) == reqCol)[0][0] for reqCol in ["TIME", fluxType, fluxType+"_ERR"]]
    units = [data.columns.units[col] for col in cols]
    # return the data
    return time, flux, ferr, units

def _lightcurves_(target, loc="/raw_data/tess_data/", stale_time=7):
    # filename for this target
    lcfname = loc+target.replace(" ","_")+"_lc.csv"
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
        # fetch all the lightcurves that aren't the fast-lc
        mask = np.array([fname.endswith("lc.fits") and not fname.endswith("fast-lc.fits") for fname in res["productFilename"]])
        res = res[mask]
        # output arrays
        times, fluxes, ferrs, units = np.array([]), np.array([]), np.array([]), ()
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
        # save the dataframe
        saveDataFrame(df, lcfname)
        # and return it too
        return df
    # if there was no observation data
    raise Exception("No LightCurve Data for {}".format(target))

def tessMidTransits(target, loc, stale_time):
    # lowercase the target name
    target = target.lower()
    # ensure the tess lightcurve folder exists
    makeFolder(loc)
    # fetch the lightcurve data
    df = _lightcurves_(target, loc, stale_time)
    print("LC Found")
    # decompose into time, flux, and error
    t, f, ferr = np.array(df["time"]), np.array(df["flux"]), np.array(df["flux_e"])
    # Period and t0 from Anderson et al. (201X):
    P,t0 =  4.7423729 ,  2457376.68539 - 2457000
    # Get phases --- identify out-of-transit (oot) times by phasing the data
    # and selecting all points at absolute phases larger than 0.02:
    phases = juliet.utils.get_phases(t, P, t0)
    idx_oot = np.where(np.abs(phases)>0.02)[0]
    # Save the out-of-transit data into dictionaries so we can feed them to juliet:
    times, fluxes, fluxes_error = {},{},{}
    times['TESS'], fluxes['TESS'], fluxes_error['TESS'] = t[idx_oot], f[idx_oot], ferr[idx_oot]
    # First define the priors:
    priors = {}

    # Same priors as for the transit-only fit, but we now add the GP priors:
    params = ['P_p1','t0_p1','r1_p1','r2_p1','q1_TESS','q2_TESS','ecc_p1','omega_p1',\
              'rho', 'mdilution_TESS', 'mflux_TESS', 'sigma_w_TESS', \
              'GP_sigma_TESS', 'GP_rho_TESS']

    dists = ['normal','normal','uniform','uniform','uniform','uniform','fixed','fixed',\
             'loguniform', 'fixed', 'normal', 'loguniform', \
             'loguniform', 'loguniform']

    hyperps = [[4.7,0.1], [1329.9,0.1], [0.,1], [0.,1.], [0., 1.], [0., 1.], 0.0, 90.,\
               [100., 10000.], 1.0, [0.,0.1], [0.1, 1000.], \
               [1e-6, 1e6], [1e-3, 1e3]]

    # Populate the priors dictionary:
    for param, dist, hyperp in zip(params, dists, hyperps):
        priors[param] = {}
        priors[param]['distribution'], priors[param]['hyperparameters'] = dist, hyperp

    print("Params set")

    times['TESS'], fluxes['TESS'], fluxes_error['TESS'] = t,f,ferr
    dataset = juliet.load(priors=priors, t_lc = times, y_lc = fluxes, \
                          yerr_lc = fluxes_error, GP_regressors_lc = times, \
                          out_folder = '{}_transitGP'.format(target.replace(" ","_")), verbose = True)
    print("Dataset saved")

    results = dataset.fit()
    # Extract full model:
    transit_plus_GP_model = results.lc.evaluate('TESS')

    # Deterministic part of the model (in our case transit divided by mflux):
    transit_model = results.lc.model['TESS']['deterministic']

    # GP part of the model:
    gp_model = results.lc.model['TESS']['GP']

    # Now plot. First preambles:
    fig = plt.figure(figsize=(12,4))
    gs = gridspec.GridSpec(1, 2, width_ratios=[2,1])
    ax1 = plt.subplot(gs[0])

    # Plot data
    ax1.errorbar(dataset.times_lc['TESS'], dataset.data_lc['TESS'], \
                 yerr = dataset.errors_lc['TESS'], fmt = '.', alpha = 0.1)

    # Plot the (full, transit + GP) model:
    ax1.plot(dataset.times_lc['TESS'], transit_plus_GP_model, color='black',zorder=10)

    ax1.set_xlim([1328,1350])
    ax1.set_ylim([0.96,1.04])
    ax1.set_xlabel('Time (BJD - 2457000)')
    ax1.set_ylabel('Relative flux')

    ax2 = plt.subplot(gs[1])

    # Now plot phase-folded lightcurve but with the GP part removed:
    ax2.errorbar(phases, dataset.data_lc['TESS'] - gp_model, \
                 yerr = dataset.errors_lc['TESS'], fmt = '.', alpha = 0.3)

    # Plot transit-only (divided by mflux) model:
    idx = np.argsort(phases)
    ax2.plot(phases[idx],transit_model[idx], color='black',zorder=10)
    ax2.yaxis.set_major_formatter(plt.NullFormatter())
    ax2.set_xlabel('Phases')
    ax2.set_xlim([-0.03,0.03])
    ax2.set_ylim([0.96,1.04])
    # fit transits to the lightcurve data using juliet
    #midTransitTimes = fitTransits(df["Times"])
    #plt.errorbar(t, f, yerr=ferr, fmt='.', alpha=0.2)
    #plt.xlabel('Time (BJD)')
    #plt.ylabel('Relative flux')
    #plt.xlim((2457000 + 1325, 2457000+1355))
    #plt.ylim((0.94, 1.06))
    plt.show()
