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
