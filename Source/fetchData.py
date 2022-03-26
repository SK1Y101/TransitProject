# Fetch required exoplanet data. See: Readme/#Programming steps

# Python modules
import pandas as pd

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
wait(1000)
ws.fetchExoClockData("https://www.exoclock.space/database/literature", sourceName="ExoClock literature")
ws.fetchExoClockData("https://www.exoclock.space/database/etd", sourceName="ExoClock ETD")

# Exoplanet Archive (System info)

# Exoplanet Catalogue (System info & digital renders for pretty graphs)
