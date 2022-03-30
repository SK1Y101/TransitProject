# python modules
from datetime import datetime
import pandas as pd
import os, time, re

# -----< Methods >-----

def loadDataFrame(fpath):
    ''' Load a pandas.DataFrame object to a given file.
        fpath: The file path to load the dataframe from.
        *Note*  This will infer the load type from any given extension,
                else will default to loading with csv.

        Compatible extensions:
        .csv .h5 .html .json .pkl .xlsx'''
    # if the path begins with a slash, append the current working directory
    if fpath[0] == "/":
        fpath = os.getcwd()+fpath
    # fetch the extension
    ext = os.path.splitext(fpath)[1]
    # if there isn't one, assume csv
    if not ext:
        ext, fpath = ".csv", fpath+".csv"
    # table of allowed extensions
    exts = {".csv":pd.read_csv, ".h5":pd.read_hdf, ".html":pd.read_html,
            ".json":pd.read_json, ".pkl":pd.read_pickle, ".xlsx":pd.read_excel}
    # check for extension
    if ext in exts:
        # load the dataframe from the file
        df = exts[ext](fpath)
        # and return
        return df

def saveDataFrame(df, fpath):
    ''' Save a pandas.DataFrame object to a given file.
    df: The dataframe to store.
        fpath: The file path to store the dataframe as.
        *Note*  If a compatible file extension is given (ie: h5, json) this will
                store in that format, else will default to csv.

        Compatible extensions:
        .csv .h5 .html .json .pkl .xlsx'''
    # if the path begins with a slash, append the current working directory
    if fpath[0] == "/":
        fpath = os.getcwd()+fpath
    # fetch the extension
    ext = os.path.splitext(fpath)[1]
    # if there isn't one, assume csv
    if not ext:
        ext, fpath = ".csv", fpath+".csv"
    # remove unnamed columns, and any of the index columns
    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
    df = df.loc[:, ~df.columns.str.contains('index')]
    # don't store duplaicates
    df = df.loc[df.astype(str).drop_duplicates().index]
    # table of allowed extensions
    exts = {".csv":df.to_csv, ".h5":lambda a : df.to_hdf(a, key="df", mode="w"),
            ".html":df.to_html, ".json":df.to_json, ".pkl":df.to_pickle, ".xlsx":df.to_excel}
    # check for extension
    if ext in exts:
        # load the dataframe into the file
        exts[ext](fpath, index=False)

def appendDataFrame(df, loc, sortby=None, ascending=False):
    ''' Save a pandas.DataFrame object to a given file.
        df: The dataframe to store.
        loc: The location of the saved dataFrame to append to.
        sortby: The column name to sort the dataset by. will use the first column by default'''
    if not sortby:
        sortby = df.columns[0]
    # check for any old data
    if checkForFile(loc):
        old = loadDataFrame(loc)
        # concatenate the new and old dataframe
        old = pd.concat([df, old], ignore_index=True)
        # sort the dataset and remove duplicates
        old.sort_values(by=sortby, inplace=True, ascending=ascending)
        old.reset_index(inplace=True)
        old = old.loc[old.astype(str).drop_duplicates().index]
        # if the dataframe hasn't changed, we're fine
        if df.equals(old):
            return
        # otherwise, update our df reference to store
        df = old
    # save the data
    saveDataFrame(df, loc)

def readFromFile(fileName):
    ''' Read the contenmt of a file using the standard python functions.
        fileName: the file to read from. '''
    # if the path begins with a slash, append the current working directory
    if fileName[0] == "/":
        fileName = os.getcwd()+fileName
    # opne the file and return the contents
    try:
        with open(fileName, "r") as f:
            return f.read()
    except:
        with open(fileName, "rb") as f:
            return f.read()

def writeToFile(fileName, data):
    ''' Write the specified data to a file using the standard python functions.
        fileName: the file to write to.
        data: the data to write to the file. '''
    # if the path begins with a slash, append the current working directory
    if fileName[0] == "/":
        fileName = os.getcwd()+fileName
    # write to the file
    try:
        with open(fileName, "w") as f:
            f.write(data)
    except:
        # if it failed, attempt to write bytes
        with open(fileName, "wb") as f:
            f.write(data)
    # and return the written data
    return data

def checkForFile(fileName, stale_time=7):
    ''' Check if a file exists, and whether it will be considered stale.
        fileName: The filepath to the file.
        stale_time: The number of days this file should be younger than to not be considered stale.
        useLocalStorage: Will prepend the local storage directory to the filename if set to true. '''
    # if the path begins with a slash, append the current working directory
    if fileName[0] == "/":
        fileName = os.getcwd()+fileName
    # return nothing if the file doesnt exist
    if not os.path.exists(fileName):
        return False
    # return nothing if the file is empty:
    if not os.path.getsize(fileName):
        return False
    # return nothing if the file is stale
    if (time.time() - os.path.getmtime(fileName)) > stale_time*86400:
        return False
    # otherwise, return true
    return True

def makeFolder(loc):
    ''' Create any folders in the "loc" path that don't exist yet. '''
    # remove the filename from loc if it was given
    loc = os.path.split(loc)[0]
    # if the path begins with a slash, append the current working directory
    if loc[0] == "/":
        loc = os.getcwd()+loc
    # if the folder dosent exist
    if not os.path.exists(loc):
        # split the filepath by directory
        dirs = loc.split("/")
        # itterate backwards over the directory
        for x in range(len(dirs))[::-1]:
            # if this directory exists, break here
            if os.path.exists("/".join(dirs[:x])):
                break
        # and now create all of the directories that didn't exist
        for y in range(x, len(dirs)):
            try:
                os.mkdir("/".join(dirs[:y+1]))
            except:
                pass
    # return the folder existing
    return os.path.exists(os.path.split(loc)[0])

''' An absolutely magical regex function to get all floats in a string. '''
def findFloats(txt=""):
    return re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", txt)
    # Do I understand how this works? Vaguely.
    # But it was sourced from an excellent stackoverflow answer that does explain it anyway
    # https://stackoverflow.com/a/45001796
