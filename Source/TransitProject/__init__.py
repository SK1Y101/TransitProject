# python modules
from datetime import datetime
from tqdm import tqdm
import pandas as pd
import os, time, re
import numpy as np

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

def inDataFrame(loc, target, col=None):
    ''' check if a target value exists in a dataframe.
        loc: the dataframe file to load and check.
        target: the target value to search for.
        col: the column to search in the dataframe for. '''
    # fetch the dataframe
    df = loadDataFrame(loc)
    # ensure the column is valid
    if not col:
        col = df.columns[0]
    # search for the target in the column
    return target in list(df[col])

def moveRowToTop(df, targetIdx):
    ''' Move a given row to the top of a dataframe.
        df: The dataframe to manipulate
        targetIdx: The row index to move to the top '''
    # append an empty row
    df.loc[len(df.index)] = np.nan
    # shift down
    df = df.shift(1)
    # copy target to first row
    df.iloc[0] = df.iloc[targetIdx+1]
    # remove the old value and return
    return df[df.index != targetIdx+1]

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

def noneIfInf(x):
    ''' return None if a value is infinity, otherwise return the value. '''
    if np.isinf(x):
        return None
    return x

def animated_loading_function(func, *args, name="Waiting"):
    ''' Execute a function with an animation showing the loading progress.
        Usefull for long execution that does not have an output of it's own.
        func:   The function to execute, passed as a lambda. '''
    # import threading and Queue
    from multiprocessing.pool import ThreadPool as Pool
    from time import sleep as wait
    from itertools import cycle
    # function to print the current output.
    def __showElapsed__(name, elapsed, spinner=" ", fin=False):
        print("{}: {:02.0f}:{:02.0f} {}".format(name, *divmod(elapsed, 60), spinner),
              end="\n"*fin+"\r"*(not fin))
    # start the pool
    with Pool(1) as p:
        # setup the counting variables
        el, wt, cyc = 0, .1, cycle(["⢿", "⣻", "⣽", "⣾", "⣷", "⣯", "⣟", "⡿"])
        # execute the function
        r = p.apply_async(func, args)
        # wait until it is complete
        while not r.ready():
            wait(wt)
            el += wt
            # show the ongoing output
            __showElapsed__(name, el, next(cyc))
    # show the completion output
    __showElapsed__(name, el, fin=True)
    # return any of the outputs
    return r.get()

def toItterableInput(*inputs, onlyUnique=True, keep=(None,)):
    ''' Will decompose any number of inputs into a list of arguments for a
        function call with multiprocessing.pool().starmap
        *inputs: any number of inputs to the function to multiprocess.
        onlyUnique: Whether the output array will only contain unique arguments.
        keep: the index of any inputs that should be kept as their tuple.'''
    # the types of inputs we want to expand
    tupleTypes = (tuple, list, np.ndarray)
    # and ensure our input is a list, rather than a tuple
    inputs = list(inputs)
    # if desired, sort to only unique inputs
    if onlyUnique:
        # for each input
        for x in range(len(inputs)):
            # try to condense down to only unique values
            try:
                # if we aren't keeping this input as is
                if x not in keep:
                    inputs[x] = np.unique(inputs[x], axis=0)
            # otherwise, keep as is
            except:
                pass
    # return if an input should be split or not
    def splitInput(inputs, x, keep):
        return (x not in keep) and isinstance(inputs[x], tupleTypes)
    # length of input tuple
    inputNum = range(len(inputs))
    # compute required number of inputs
    num = int(np.product([len(inputs[x]) for x in inputNum if splitInput(inputs, x, keep)]))
    # create the output by itteration of the input
    out = [[inputs[y][x%len(inputs[y])] if splitInput(inputs, y, keep) else inputs[y] \
            for y in inputNum] for x in range(num)]
    # and return the array of inputs
    return out

def parallelJob(function, inputs, outType=None, workers=None, reserveCpu=True, tqdmLeave=True):
    ''' Will execute a given function as a multiprocessing pool.
        function: The function to parallelise.
        inputs: The inputs to the function.
        workers: The number of worker processes to spawn. Defaults to half of the number of cpu cores
        reserveCpu: If true, the multiprocessing pool will ensure at least one CPU core is idle for other OS tasks.
        outType: The type that the final output should be. Defaults to list. Accepts: list, tuple, dict, np.array.
        tqdmLeave: Whether to leave the overarchive tqdm bar when completed.'''
    # imports required for this to work
    from time import sleep as wait
    import multiprocessing as mp
    # if worker number isn't give, set it to the job length or half cpu count.
    # if it was given, set it to the given number not exceeding the maximum cpu count
    workers = min(len(inputs), mp.cpu_count()//2) if not workers else min(workers, mp.cpu_count()-reserveCpu)
    # the number of completed jobs
    comp_0, comp_1 = 0, 0

    # create the tqdm bar for the output
    with tqdm(total=len(inputs), desc="Parallel {}".format(function.__name__), \
              leave=tqdmLeave, smoothing=0) as bar:
        # create the multiprocessing pool
        with mp.Pool(workers) as p:
            # create the worker tasks
            results = p.starmap_async(function, inputs, chunksize=1)
            # show the tqdm output while waiting for things to finish
            while not results.ready():
                # update the number of completed jobs
                bar.n = bar.total-results._number_left
                bar.update(0)
        # ensure the bar is completed when the job is done
        bar.n = bar.total
        # update the bar
        bar.refresh()
    # fetch the results of the pool
    res = results.get()
    # if the outType is defined, and is applicable
    if outType and outType in [list, tuple, dict, np.array]:
        return outType(res)
    # return the output
    return res
