# Fetch required exoplanet data. See: Readme/#Programming steps
import TransitProject as tp
import TransitProject.webScraping as ws

# Examples for fetching

# ETD (Lightcurves)

# Exoclock (Lightcurves & system info)
# Access system information

def fetch_from_ExoClock(url=""):
    # fetch the page source code
    page = ws.bs4("https://www.exoclock.space/database/observations")
    # locate the exoplanet rows from the main table
    rows = page.find("tbody").find_all("tr")
    # search the table for each exoplanet
    exoplanets, this_data, this_exoplanet = {}, [], ""
    for x in rows:
        # if the row has an id tag, its the title of an exoplanet
        if x.has_attr("id"):
            # if we had data from the previous exoplanet
            if this_data:
                # store it
                exoplanets[this_exoplanet] = this_data
            # clear the temporary exoplanet data, and fetch the new name
            this_data, this_exoplanet = [], x.find("button").text
            continue
        try:
            # fetch the columns
            cols = x.find_all("td")
            # fetch the o-c time, and observation date
            obs, oc, ocerror = cols[0].text, *cols[-2].text.split(" ± ")
            this_data.append([obs, oc, ocerror])
        except:
            pass
    print(exoplanets)

fetch_from_ExoClock("https://www.exoclock.space/database/observations")

# Exoplanet Archive (System info)

# Exoplanet Catalogue (System info & digital renders for pretty graphs)
