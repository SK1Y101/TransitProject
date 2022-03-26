# import python modules
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
        name = name[name.index("/"):]
    loc += name.capitalize().replace(" ", "")+".csv"
    # ensure the location exists
    makeFolder(loc)
    # check for any old data
    if checkForFile(loc):
        old = loadDataFrame(loc)
        old = pd.concat([df, old], ignore_index=True)
        old.sort_values(by="date", inplace=True, ascending=False)
        old.reset_index(inplace=True)
        if df.equals(old):
            return
        # otherwise, update our df reference to store
        df = old
    # save the data
    saveDataFrame(df, loc)

def fetchExoClockData(url, loc="/raw_data/midTransitTimes/", stale_time=7, sourceName=""):
    ''' Scrape the ExoClock website for exoplanet o-c data, and store as locally accesible csv data.
        url: The webpage to scrape.
        loc: the local location to store the data
        stale_time: the time in days to use a locally stored copy of the page instead.
        sourceName: The title to log in the source column of the data. '''
    # check for local storage file
    makeFolder(loc)
    # fetch data
    data_dict = fetch_from_ExoClock(url, stale_time)
    # itterate on data
    for x in tqdm(data_dict, desc="Fetching {}".format(sourceName), miniters=0):
        # convert to a DataFrame
        df = pd.DataFrame(data_dict[x])
        # and save
        saveTransitTimes(df)

def fetch_from_ExoClock(url="", sttransit_dateale_time=7):
    ''' Scrape the ExoClock website for exoplanet o-c data.
        url: The webpage to scrape.
        stale_time: the time in days to use a locally stored copy of the page instead. '''
    # fetch the page source code
    page = bs4(url, stale_time=stale_time)
    # locate the exoplanet rows from the main table
    rows = page.find("tbody").find_all("tr")
    # fetch which column the observer information is held
    observer_col = -3 if "literature" in url else 1
    # search the table for each exoplanet
    exoplanets, this_data, this_exoplanet = {}, [], ""
    for x in tqdm(rows):
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
        transit_data = sourceFromETDTable(subpage)
        if transit_data.empty:
            continue
        #transit_data = sourceFromETDAscii(subpage, suburl)
        # add the source name
        transit_data["source"]=sourceName
        # convert the date from HJD to YYYY-MM-DD, and the OC from days to minutes
        transit_data["date"] = pd.to_datetime(transit_data["date"], unit="D", \
                               origin=pd.Timestamp("1858-11-16 12:00")).dt.strftime('%Y-%m-%d')
        transit_data["oc"] = round(pd.to_numeric(transit_data["oc"])*1440,4)
        transit_data["oce"] = round(pd.to_numeric(transit_data["oce"])*1440,4)
        # save the data
        saveTransitTimes(transit_data, name=thisExoplanet[1])
    # remove the index column from ETD
    df.drop(" ", axis=1, inplace=True)
    # save the dataframe
    saveDataFrame(df, "/raw_data/{}Ephemerides".format(sourceName))

def sourceFromETDTable(page):
    ''' Fetch midtransit times from the main ETD Table on the page.
        page: the bs4 page source.'''
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
        suburl: The url of the homepage.
        subpage: the bs4 element of the current page.
        stale_time: as everything else. '''
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
