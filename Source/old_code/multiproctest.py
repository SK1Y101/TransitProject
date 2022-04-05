from multiprocessing import Pool
from time import sleep as wait
import numpy as np

def testfunc(sim, tstep, transits=1000, i=0):
    wait(1.5)
    return i

sim = np.random.rand(5, 10)

itterable = [(sim[x], 100, 1000, x) for x in range(len(sim))]

with Pool() as p:
    r=p.starmap_async(testfunc, itterable)
    while not r.ready():
        print(r._number_left)
        wait(0.5)
    print(r.get())
