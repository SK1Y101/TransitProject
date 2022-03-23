import os, time

# -----< Methods >-----

def readFromFile(fileName):
    ''' Read the contenmt of a file using the standard python functions.
        fileName: the file to read from. '''
    # if the path begins with a slash, append the current working directory
    if fileName[0] == "/":
        fileName = os.getcwd()+fileName
    # opne the file and return the contents
    with open(fileName, "r") as f:
        return f.read()

def writeToFile(fileName, data):
    ''' Write the specified data to a file using the standard python functions.
        fileName: the file to write to.
        data: the data to write to the file. '''
    # if the path begins with a slash, append the current working directory
    if fileName[0] == "/":
        fileName = os.getcwd()+fileName
    # write to the file
    with open(fileName, "w") as f:
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
