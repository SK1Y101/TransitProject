# Fetch required exoplanet data. See: Readme/#Programming steps
import TransitProject as tp
import TransitProject.webScraping as ws
from time import sleep as wait

# replace <br/> tags with newlines in a bs4 object
def replace_breaks(elem):
    for br in elem.find_all("br"):
        br.replace_with("\n")

# Examples for fetching

# ETD (Lightcurves)

# Exoclock (Lightcurves & system info)
# Access system information

def fetchExoClockData(url, stale_time=7):
    # this function is currently useless, but will eventually be used to handle local storage
    # fetch data
    return fetch_from_ExoClock(url, stale_time)

def fetch_from_ExoClock(url="", stale_time=7):
    # fetch the page source code
    page = ws.bs4(url, stale_time=stale_time)
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
        # fetch the columns
        cols = x.find_all("td")
        # if we have data this column
        if len(cols) > 2:
            thiscol = {}
            # fetch the o-c time
            thiscol["date"] = cols[0].text.strip()
            # observer and observatory
            replace_breaks(cols[1])
            thiscol["observer"], thiscol["observatory"] = cols[1].text.strip().split("\n")
            # Telescope, camera, filter, exposure
            replace_breaks(cols[2])
            cols2 = cols[2].text.split("\n")
            thiscol["telescope"] = cols2[0].strip()
            thiscol["camera"], thiscol["filter"], thiscol["exposure"] = cols2[1].replace("\\n", "").strip().split(" / ")[-3:]
            # oc and oc error
            thiscol["oc"], thiscol["ocerror"] = cols[3].text.strip().split(" ± ")
            # add this data to our exoplanet list
            this_data.append(thiscol)
    # return the dictionary
    return exoplanets

dict = fetchExoClockData("https://www.exoclock.space/database/observations")
#print(dict)
print(dict["HAT-P-13b"])

# Exoplanet Archive (System info)

# Exoplanet Catalogue (System info & digital renders for pretty graphs)
