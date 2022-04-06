from multiprocessing import Pool
from time import sleep as wait
from tqdm import tqdm, trange
import multiprocessing
import numpy as np

def testfunc(sim, tstep, transits=1000):
    workerId = multiprocessing.current_process()._identity[0]
    for x in trange(10+workerId, leave=None, position=workerId):
        wait(0.1)
    return sim

def toItterableInput(*inputs, onlyUnique=True):
    # the types of inputs we want to itterate
    tupleTypes = (tuple, list, np.ndarray)
    # if desired, sort to only unique inputs
    if onlyUnique:
        # for each input
        for x in range(len(inputs)):
            # try to condense down to only unique values
            try:
                inputs[x] = np.unique(inputs[x], axis=0)
            # otherwise, keep as is
            except:
                pass
    # compute required number of inputs
    num = int(np.product([len(x) for x in inputs if isinstance(x, tupleTypes)]))
    # create the output by itteration of the input
    out = [[y[x%len(y)] if isinstance(y, tupleTypes) else y \
            for y in inputs] for x in range(num)]
    # and return the array of inputs
    return out

class a:
    pass

sims = [a(), a(), a(), a(), a()]
itterable = toItterableInput(sims, [0.1, 0.1, 0.1], 1000)
maxjobs = len(itterable)
_complete=0

with tqdm(total = maxjobs, leave=False, desc="Multiprocessing pool", smoothing=0) as bar:
    with Pool(8) as p:
        r=p.starmap_async(testfunc, itterable, chunksize=1)
        while not r.ready():
            completed = maxjobs-r._number_left
            bar.update((completed)-_complete)
            _complete=completed
    bar.n = bar.total
    bar.refresh()
    print(r.get())
