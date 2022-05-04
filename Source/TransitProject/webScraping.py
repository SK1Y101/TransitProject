# import python modules
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive as archive
from urllib.request import urlopen
from bs4 import BeautifulSoup
from tqdm import tqdm
import pandas as pd
import json, io

# import my modules
from . import *

# -----< Methods >-----

def saveDataDict(data, loc, stale_time=7):
    ''' Take a dictionary of data and save as a pandas dataframe for future usage.
        data: a dictionary of data to store, or a json url to source the dictionary data from.
        loc: the filepath within which to store the data.'''
    # ensure the location exists
    makeFolder(loc)
    # continue if we don't have stored data available
    if checkForFile(loc, stale_time):
        return
    # if we were given a url string
    if isinstance(data, str):
        # source the data
        data = loadJsonURL(data)
    # convert to a pandas dataframe
    df = pd.DataFrame.from_dict(data, orient="index")
    # and store at the given location
    saveDataFrame(df, loc)

def loadJsonURL(url="", stale_time=7, loc="/raw_data/webpage_source/"):
    ''' Load JSON Data from a remote url.
        url: the url to fetch data from.
        stale_time: the time in days to use a locally stored copy instead. '''
    source = fetchUrlSource(url, stale_time, loc=loc)
    return json.loads(source)

def fetchUrlSource(url="", stale_time=7, loc="/raw_data/webpage_source/"):
    ''' Will fetch the source code from a given URL.
        url: The URL to fetch the source code from.
        stale_time: the time in days to consider the local source for consideration. '''
    # remove spaces from the url if any are given
    url = url.replace(" ","")
    # convert the url to a string for filestorage
    file_dir = "{}{}.html".format(loc, url.replace("https://", "").replace("/", "_"))
    # if a local version exists, and is not considered stale
    if checkForFile(file_dir, stale_time=stale_time):
        # then return the contents of that
        return readFromFile(file_dir)
    # otherwise, fetch the source
    source = urlopen(url).read()
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

def bs4(url="", features="html5lib", stale_time=7, loc="/raw_data/webpage_source/"):
    ''' Will fetch the source code from a given URL as a BeautifulSoup4 object
        url: The URL to fetch the source code from.
        features: The bs4 method for loading the source.
        stale_time: the time in days to consider the local source for consideration. '''
    # fetch the page source code
    source = fetchUrlSource(url, stale_time, loc=loc)
    # convert to BS4 object
    page_source = BeautifulSoup(source, features=features)
    # and return
    return page_source

def saveTransitTimes(df, loc="/raw_data/midTransitTimes/", name="", stale_time=7):
    ''' Save the midtransit time for an exoplanet.
        loc: The location to save the exoplanet details at.
        stale_time: how long to consider any stored data as accesible.
        df: the dataframe of information to store '''
    # append the exoplanet name
    if "/" in name:
        name = name[name.index("/")+1:]
    # ensure the name has no spaces, and add to the location string
    name = name.capitalize().replace(" ", "")
    loc += "/"+name+".csv"
    # ensure the location exists
    makeFolder(loc)
    # add the given data to the old dataframe
    appendDataFrame(df, loc, sortby="date")
    # add the exoplanet name to the list of exoplanets we're keeping track of
    appendDataFrame(pd.DataFrame({"name":name}, index=[0]), "/raw_data/exoplanetList.csv", ascending=True)

def fetchExoClockData(url, loc="/raw_data/midTransitTimes/", stale_time=7, sourceName=""):
    ''' Scrape the ExoClock website for exoplanet o-c data, and store as locally accesible csv data.
        url:        The webpage to scrape.
        loc:        The local location to store the data
        stale_time: The time in days to use a locally stored copy of the page instead.
        sourceName: The title to log in the source column of the data. '''
    # check for local storage file
    makeFolder(loc)
    # fetch the page source code
    page = bs4(url, stale_time=stale_time)
    # locate the exoplanet rows from the main table
    rows = page.find("tbody").find_all("tr")
    # fetch which column the observer information is held
    observer_col = -3 if "literature" in url else 1
    # search the table for each exoplanet
    exoplanets, this_data, this_exoplanet = {}, [], ""
    for x in tqdm(rows, desc="Fetching {}".format(sourceName), miniters=0):
        # if the row has an id tag, its the title of an exoplanet
        if x.has_attr("id"):
            # if we had data from the previous exoplanet
            if this_data:
            # convert to a dataframe
                this_data = pd.DataFrame(this_data)
                # add the source name
                this_data["source"] = sourceName
                # and store the midtransit times
                saveTransitTimes(this_data, name = this_exoplanet)
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

def rowToData(elem):
    ''' fetch the data from a table row, ensuring it is capitalised correctly,
        and anything after a close bracket is ignored. Also remove latinate non-
        breaking spaces if any are included. '''
    return [x.text.capitalize().replace(u'\xa0', u' ') for x in elem.find_all("td")]

def decomposeTable(source, selection_target=0):
    # find the table element
    table = source.find_all("tbody")[selection_target].find_all("tr")
    # decompose the table into a python readbale form and return the found data
    return [rowToData(row) for row in table]

def fetchETDData(url, sourceName="", stale_time=7):
    # fetch the mainurl for the webpage
    suburl = url[:url.find("index.php")]
    # fetch the websource
    source = bs4(url, stale_time=stale_time)
    # fetch the rows on the exoplanet table
    rows = source.find("div", {"class":"centerlike"}).find("tbody").find_all("tr")
    # the column name for the dataframe are first row
    cols = rowToData(rows[0])
    # split the middle column, as it's incorrect
    cols[3] = "Time span from"
    cols.insert(4, "Time span till")
    # initialise the dataframe
    df = pd.DataFrame(columns = cols)
    # for each exoplanet
    for row in tqdm(rows[1:], desc="Fetching ETD Data", miniters=0):
        # fetch this data in a nice form
        thisExoplanet = rowToData(row)
        # add this exoplanet to storage dataframe
        df.loc[len(df.index)] = thisExoplanet
        # if there are no observations, skip
        if thisExoplanet[3] == "0":
            continue
        # locate the link for each exoplanet
        exoUrl = suburl+row.find_all("td")[1].find("a")["href"]
        # navigate to that page
        subpage = bs4(exoUrl, stale_time=stale_time, loc="/raw_data/webpage_source/ETD/")
        # source the tabular data
        transit_data = sourceFromETDTable(subpage, thisExoplanet[1])
        if transit_data.empty:
            continue
        #transit_data = sourceFromETDAscii(subpage, suburl)
        # add the source name
        transit_data["source"]=sourceName
        # convert the date from HJD to YYYY-MM-DD, and the OC from days to minutes
        transit_data["date"] = HJDtoDate(transit_data["date"]).dt.strftime('%Y-%m-%d %H:%M%S')
        transit_data["oc"] = round(pd.to_numeric(transit_data["oc"])*1440,4)
        transit_data["oce"] = round(pd.to_numeric(transit_data["oce"])*1440,4)
        # some rows in the dataset were coppied incorrectly, ie: HAT-P-13b has a BJD error of 113 days,
        # looking at the literateure, this should actually be 113 seconds.
        transit_data.loc[transit_data["oce"] > 86400, "oce"] /= 86400
        # save the data
        saveTransitTimes(transit_data, name=thisExoplanet[1])
    # remove the index column from ETD
    df.drop(" ", axis=1, inplace=True)
    # save the dataframe
    saveDataFrame(df, "/raw_data/{}Ephemerides".format(sourceName))

def sourceFromETDTable(page, name=""):
    ''' Fetch midtransit times from the main ETD Table on the page.
        page:   the bs4 page source.'''
    # fetch the table
    replace_breaks(page)
    table = decomposeTable(page, -1)
    # if we have no data, return nothing
    if len(table) <= 1:
        return pd.DataFrame(columns=["date","observer","oc","oce"])
    table = pd.DataFrame(table[1:], columns=table[0])
    # split the HJDMid time into its two parts
    try:
        table[["Hjd mid (2400000 +)", "pm", "HJDmid Error"]] = table["Hjd mid (2400000 +)"].str.split(expand=True)
    except:
        # if there is not an error value given
        table[["Hjd mid (2400000 +)", "pm"]] = table["Hjd mid (2400000 +)"].str.split(expand=True)
        table["HJDmid Error"]=""
    # split the author and reference
    table[["Author", "reference"]] = table["Author & reference"].str.split("\n", expand=True)
    # and fetch only the useful columns
    table = table.loc[:, ["Hjd mid (2400000 +)", "Author", "O-c (d)", "HJDmid Error"]]
    # rename the columns
    table.columns = ["date","observer","oc","oce"]
    # lower any clearly wrong dates and convert to floats
    table["date"] = pd.to_numeric(table["date"])
    if (table["date"]>2400000).any():
        table = table[table["date"] < 2400000]
    # and return
    return table

def sourceFromETDAscii(subpage, suburl, stale_time=7):
    ''' Fetch midtransit data from the ETD ASCII table service.
        suburl:     The url of the homepage.
        subpage:    the bs4 element of the current page.
        stale_time: As everything else. '''
    # fetch the url of the ascii page
    dataUrl = suburl+subpage.find("a", string="Show data as ASCII table separated by semicolon")["href"]
    # go to the data url
    data = bs4(dataUrl, stale_time=stale_time, loc="/raw_data/webpage_source/ETD/").find("body").text
    # the data begins at the octothorp
    data = data[data.index("#"):-1]#.replace(";", ",")
    # convert to a dataframe
    data = pd.read_csv(io.StringIO(data), sep=";")
    # pull out the useful information
    transit_data = data.loc[:, ["HJDmid", "Observer", "O-C", "HJDmid Error"]]
    # rename the columns
    transit_data.columns = ["date","observer","oc","oce"]
    # return the data
    return transit_data

def fetchExoplanetArchive(table="pscomppars", loc="/raw_data/", stale_time=7):
    ''' Fetch exoplanet data from the Exoplanet Archive.
        table:      specifies which of the archive tables to use (Defaults to pscomppars)
        loc:        The directory to store the table.
        stale_time: How long to consider the local data fit for use. '''
    # fetch the location to search
    loc += table+".csv"
    # if we don't have any local data,
    if not checkForFile(loc, stale_time):
        # the function to fetch the data
        def fetchData():
            # fetch the tabular data
            rtable = archive.query_criteria(table)
            # convert to dataframe
            rtable = rtable.to_pandas()
            # and sort by planet name
            rtable.sort_values(by="pl_name", inplace=True, ascending=True)
            # store locally
            saveDataFrame(rtable, loc)
        # call the function with an animated loading output
        animated_loading_function(fetchData, name="Fetching Exoplanet Archive {} Data".format(table.capitalize()))

def fetchTESSData(loc="/raw_data/", stale_time=7):
    ''' Fetch exoplanet midtransit data from TESS.
        loc:        The directory to store the table.
        stale_time: How long to consider the local data fit for use. '''
