from multiprocessing import Pool
from time import sleep as wait
from tqdm import tqdm, trange
import multiprocessing
import numpy as np

def testfunc(sim, tstep, transits=1000):
    workerId = multiprocessing.current_process()._identity[0]
    for x in trange(10+workerId, leave=None, position=workerId):
        wait(0.1)
    return workerId

sim = np.random.rand(100, 10)

itterable = [(sim[x], 100, 1000) for x in range(len(sim))]
maxjobs = len(itterable)
_complete=0

with tqdm(total = maxjobs, leave=False, desc="Multiprocessing pool", smoothing=-.3) as bar:
    with Pool(8) as p:
        r=p.starmap_async(testfunc, itterable, chunksize=1)
        while not r.ready():
            completed = maxjobs-r._number_left
            bar.update((completed)-_complete)
            _complete=completed
    bar.n = bar.total
    bar.refresh()
    print(r.get())