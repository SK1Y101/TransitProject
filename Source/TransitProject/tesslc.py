''' Use Juliet to fetch TESS ligthcurve data

The light curve and target pixel files have the following file naming convention:

hlsp_tess-spoc_tess_phot_<tid>-<sector>_tess_v1_<type>.fits

<tid> = The full, 16-digit, zero-padded TIC ID.
<sector> = The Sector represented as a 4-digit, zero-padded string, preceded by an 's', e.g., 's0026' for Sector 26.
<type> = 'lc' for the light curve and 'tp' for the target pixel files.'''

from astroquery.mast import Observations
from matplotlib import pylab as plt
from astropy.table import join
from astropy import units as u
from astropy.io import fits
import juliet, re, os
import numpy as np

def _mastQuery_(target):
    # start with no target
    target_name=""
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
    # return the data
    return time, flux, ferr

def lightcurves(target):
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
        # fetch all the data products that end with the lightcurve filetype
        mask = np.array([fname.endswith("lc.fits") for fname in res["productFilename"]])
        res = res[mask]
        # for each file found
        for i, item in enumerate(res):
            # ensure we are using the tess data
            if item["obs_collection_"] == "TESS":
                # attempt to fetch the lightcurve data
                lcdata = _downloadlc_(item)
                # if it was successful
                if lcdata:
                    # unpack times, flux, and flux errors
                    t, f, ferr = lcdata
                    # plot temporarily
                    plt.errorbar(t, f, yerr=ferr, fmt=".")
                    plt.show()

lightcurves("tic 281541555")
