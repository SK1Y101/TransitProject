# Fetch required exoplanet data. See: Readme/#Programming steps

# Python modules
import pandas as pd

# my modules
import TransitProject as tp
import TransitProject.webScraping as ws

# Examples for fetching

# ETD (Lightcurves)

# Exoclock (Lightcurves & system info)
# Access system information

ws.fetchExoClockData("https://www.exoclock.space/database/observations", sourceName="ExoClock observation")
ws.fetchExoClockData("https://www.exoclock.space/database/literature", sourceName="ExoClock literature")
ws.fetchExoClockData("https://www.exoclock.space/database/etd", sourceName="ExoClock ETD")

# Exoplanet Archive (System info)

# Exoplanet Catalogue (System info & digital renders for pretty graphs)
