import juliet
import numpy as np
import matplotlib.pylab as plt

tid = 9725627#261136679#20096620#281541555

for sector in range(1, 100):
    print("Searching sector {:04}".format(sector))
    try:
        url = 'https://archive.stsci.edu/hlsps/tess-data-alerts/hlsp_tess-data-'+\
          'alerts_tess_phot_{:011}-s{:02}_tess_v1_lc.fits'.format(tid, sector)

        t, f, ferr  = juliet.get_TESS_data(url)
        # Put data arrays into dictionaries so we can fit it with juliet:
        times, fluxes, fluxes_error = {},{},{}
        times['TESS'], fluxes['TESS'], fluxes_error['TESS'] = t,f,ferr

        # Plot data:
        plt.errorbar(t, f, yerr=ferr, fmt='.')
        plt.xlim([np.min(t),np.max(t)])
        plt.show()

        break
    except:
        pass
