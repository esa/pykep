q#!/usr/bin/env python
# coding: utf-8

# Performing a preliminary conjunction analysis
# We build in parallel the distribution of close events


# Core imports
import pykep
import scipy
import numpy as np
import multiprocessing as mp
import pickle as pk

def perform_N(N, sats, t0, years, satellite, seed, q):
    np.random.seed(seed)
    # This will contain the results
    mins = []
    # We loop
    for i in range(N):
        epoch = pykep.epoch(t0 + np.random.rand() * years * 365.25, "mjd")
        pos_satellite, _ = np.array(satellite.eph(epoch)) / pykep.EARTH_RADIUS
        
        # All the population
        xs = [0]*len(sats)
        ys = [0]*len(sats)
        zs = [0]*len(sats)
        for i,s in enumerate(sats):
            try:
                r, v = s.eph(epoch)
                xs[i] = r[0]/pykep.EARTH_RADIUS
                ys[i] = r[1]/pykep.EARTH_RADIUS
                zs[i] = r[2]/pykep.EARTH_RADIUS
            except:
                xs[i] = 0
                ys[i] = 0
                zs[i] = 0
                continue
        pos = np.array((xs,ys,zs)).T
        D = scipy.spatial.distance.cdist([pos_satellite],pos)
        sat_id = np.argmin(D[0])
        D[0,sat_id] = 3
        mins.append(min(D[0]))
        with open("seed_"+str(seed)+".pk", "wb") as file:
            pk.dump(mins, file)
    q.put(mins)

if __name__ == "__main__":  # confirms that the code is under main function
    # Avoids to use the buggy fork
    ctx = mp.get_context('spawn')
    # Here we import all the objects in the catalogue all.tle into pykep.planet objects
    sats = pykep.util.read_tle(tle_file="all.tle", with_name=True)

    # We define the satellite Sentinel 3B from its tle.
    sentinel3b_l1="1 43437U 18039A   20145.62789730 +.00000000 +00000-0 +18124-4 0  9996"
    sentinel3b_l2="2 43437 098.6292 212.6693 0000789 104.4098 255.7170 14.26736461108341"
    sentinel3b = pykep.planet.tle(sentinel3b_l1, sentinel3b_l2)

    # We make the Monte-Carlo simulation
    # Number of runs, should be at least 100000 (hours to compute on a single CPU)
    N = 100
    # This is the lower bound for the epoch
    t0 = 58991.90384230018
    # This is the number of years to sample
    years = 3

    procs = []
    q = mp.Queue()

    for i in range(mp.cpu_count()):
        seed = np.random.randint(10000000,100000000)
        proc = ctx.Process(target=perform_N, args=(N, sats, t0, years, sentinel3b, seed, q))
        procs.append(proc)
        proc.start()
        print("Process started: ", proc.pid)
    # complete the processes
    for proc in procs:
        proc.join()

    res = []
    size = q.qsize()
    for i in range(size):
        res = res + q.get()

    with open("res.pk", "wb") as file:
        pk.dump(res, file)

