# Fetch required exoplanet data. See: Readme/#Programming steps

# Python modules
import pandas as pd
from tqdm import tqdm, trange

# my modules
import TransitProject as tp
import TransitProject.webScraping as ws

'''
Goals: Collect midtransit times from ETD and ExoClock
       Collect System information from exoplanet archive
       Collect ephemerides from Exoclock and exoplanet archive.
       Only download data occasionally (minimise network overhead)
'''

'''
# ETD (Lightcurves)
ws.fetchETDData("http://var2.astro.cz/ETD/index.php?lang=en", sourceName="ETD")

# Exoclock (Lightcurves & system info)
# Access system information
ws.saveDataDict("https://www.exoclock.space/database/planets_json", "/raw_data/exoClockEphemerides.csv")
# fetch midtransit data (Old style)
ws.fetchExoClockData("https://www.exoclock.space/database/observations", sourceName="ExoClock observation")
ws.fetchExoClockData("https://www.exoclock.space/database/literature", sourceName="ExoClock literature")
ws.fetchExoClockData("https://www.exoclock.space/database/etd", sourceName="ExoClock ETD")
# fetch midtransit data (New style)

# Exoplanet Archive
# Composite information
ws.fetchExoplanetArchive("pscomppars")
# all information
ws.fetchExoplanetArchive("ps")'''

# Transiting Exoplanets Survey Satelite
ws.fetchTESSData()
