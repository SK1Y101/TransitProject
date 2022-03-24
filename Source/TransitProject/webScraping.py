# import python modules
from urllib.request import urlopen
from bs4 import BeautifulSoup

# import my modules
from . import *

# -----< Methods >-----

def fetchUrlSource(url="", stale_time=7):
    ''' Will fetch the source code from a given URL.
        url: The URL to fetch the source code from.
        stale_time: the time in days to consider the local source for consideration. '''
    # convert the url to a string for filestorage
    file_dir = "/raw_data/webpage_source/{}.html".format(url.replace("https://", "").replace("/", "_"))
    # if a local version exists, and is not considered stale
    if checkForFile(file_dir, stale_time=stale_time):
        # then return the contents of that
        return readFromFile(file_dir)
    # otherwise, fetch the source
    source = str(urlopen(url).read())
    # ensure the folder for the file exists
    makeFolder(file_dir)
    # store the fetched source, and return it
    return writeToFile(file_dir, source)

def replace_breaks(elem):
    ''' replace <br/> tags with newlines in a bs4 object.
        elem: the bs4 element to replace newline tags in.'''
    for br in elem.find_all("br"):
        br.replace_with("\n")

def remove_val(elem, val=""):
    ''' remove value from bs4 object.
        elem: the bs4 element to remove things from.
        val: the thing to remove. '''
    for br in elem.find_all(val):
        br.replace_with("")

def bs4(url="", features="html5lib", stale_time=7):
    ''' Will fetch the source code from a given URL as a BeautifulSoup4 object
        url: The URL to fetch the source code from.
        features: The bs4 method for loading the source.
        stale_time: the time in days to consider the local source for consideration. '''
    # fetch the page source code
    source = fetchUrlSource(url, stale_time)
    # convert to BS4 object
    page_source = BeautifulSoup(source, features=features)
    # and return
    return page_source

def fetchExoClockData(url, loc="/raw_data/midTransitTimes/", stale_time=7, sourceName=""):
    # check for local storage file
    makeFolder(loc)
    # fetch data
    data_dict = fetch_from_ExoClock(url, stale_time)
    # itterate on data
    for x in data_dict:
        # add the source name to every value
        for y in data_dict[x]:
            y["source"] = sourceName
        # convert to a DataFrame
        df = pd.DataFrame(data_dict[x])
        # check if one has been stored locally already
        if checkForFile(loc+x+".csv"):
            # fetch the old data
            old_df = loadDataFrame(loc+x+".csv")
            # compare with the collected data, adding any new references
            old_df = pd.concat([df, old_df], ignore_index=True)#.drop_duplicates(keep=False)
            # and reoder by the date column
            old_df.sort_values(by="date", inplace=True, ascending=False)
            old_df.reset_index(inplace=True)
            # if the dataframe hasn't changed, skip to the next reference
            if df.equals(old_df):
                continue
            # otherwise, update our df reference to store
            df = old_df
        # remove unnamed columns, and any of the index columns
        df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
        df = df.loc[:, ~df.columns.str.contains('index')]
        # store the dataframe
        saveDataFrame(df, loc+x+".csv")

def fetch_from_ExoClock(url="", stale_time=7):
    # fetch the page source code
    page = bs4(url, stale_time=stale_time)
    # locate the exoplanet rows from the main table
    rows = page.find("tbody").find_all("tr")
    # fetch which column the observer information is held
    observer_col = -3 if "literature" in url else 1
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
            # observation date
            thiscol["date"] = cols[0].text.strip()
            # observer
            replace_breaks(cols[observer_col])
            thiscol["observer"] = cols[observer_col].text.strip().split("\n")[0]
            # oc and oc error
            thiscol["oc"], thiscol["oce"] = [x.strip() for x in findFloats(cols[-2].text)]
            # add this data to our exoplanet list
            this_data.append(thiscol)
    # return the dictionary
    return exoplanets
