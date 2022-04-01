# Fetch required exoplanet data. See: Readme/#Programming steps

# Python modules
import pandas as pd
from tqdm import tqdm, trange

# my modules
import TransitProject as tp
import TransitProject.webScraping as ws

# Examples for fetching

# ETD (Lightcurves)
ws.fetchETDData("http://var2.astro.cz/ETD/index.php?lang=en", sourceName="ETD")

# Exoclock (Lightcurves & system info)
# Access system information
ws.saveDataDict("https://www.exoclock.space/database/planets_json", "/raw_data/exoClockEphemerides.csv")
# fetch midtransit data
ws.fetchExoClockData("https://www.exoclock.space/database/observations", sourceName="ExoClock observation")
ws.fetchExoClockData("https://www.exoclock.space/database/literature", sourceName="ExoClock literature")
ws.fetchExoClockData("https://www.exoclock.space/database/etd", sourceName="ExoClock ETD")

# Exoplanet Archive
# Composite information
ws.fetchExoplanetArchive("pscomppars")
# all information
ws.fetchExoplanetArchive("ps")
