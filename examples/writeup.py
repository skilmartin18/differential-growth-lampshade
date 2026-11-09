import matplotlib.pyplot as plt
import sys
import logging
import numpy as np
from scipy import signal
from datetime import datetime
import os
import pickle
from copy import deepcopy
import open3d as o3d

sys.path.append('../')
import particle
import sim as s

log = logging.getLogger("logger")
log.setLevel(logging.INFO)

# this create a vase with a consistent cross section, based on a single 
# layer from the the whole vase
# they may be basic duplications of layers, or modulated forms (with twisting and changing radii)

LOAD_PATH = '../saves/med-dense/layers.pkl'
SAVE_ID = '-med-dense'

SHOW = False


### FIGURE 1 ###
## create clump

if False:
    cl = particle.clump()

    # create sim using clump
    sim = s.growth_sim(cl, threeD=False, saveToFile=True)

    # load clump from file
    sim.load_clump_from_layers(LOAD_PATH)
    sim.clump = sim.layers[len(sim.layers)-1]

    # sim.draw_clump(colour='olive')
    # sim.show()

    sim.run_clump_smoothing(10, smoothing_function=sim.gaussian_smooth, smoothingKwargs={'points': 3})
    sim.convert_to_polar()
    sim.normalise_to_circle()
    sim.convert_to_cartesian()
    # sim.polar_update_clump()
    sim.draw_clump(colour='olive')
    sim.show()

### CLUMP PARAMS ###
ITERATIONS = 250
_perturbation_func = particle.perturbation(prob=0.5, magnitude=2)
_general_force = particle.force(rDistance=75, aDistance=100, k=30, neighbors=2, clamp = 1e999)

SIZE = 6
RADIUS = 40
INSERT_DIST = 15

MASS = 1
FRICTION = 0.8
LIFE = 99999

particleParams = {
    'mass': MASS,
    'friction': FRICTION,
    'life': LIFE
}

growingParams = {
    "growthZoneIterations": 9999,
    "life": 50,
    "growth": 2,
    "coverage": 0,
}
### FIGURE 2 ###
if True:
    

    ### RUN SIMULATION ###

    log.setLevel(logging.INFO)

    ## create clump
    cl = particle.clump(size=SIZE, radius=RADIUS, insertDistance=INSERT_DIST, _perturbation=_perturbation_func, _force=_general_force, particleKwargs= particleParams)

    # create sim using clump
    sim = s.growth_sim(cl, threeD=False, growingKwargs=growingParams, saveToFile=True)
    sim.grow_clump(ITERATIONS,showProgress=False)
    # sim.draw_clump(colour='olive',showParticles=True)
    sim.show()


if False:
    ## create clump
    cl = particle.clump(size=4, radius=RADIUS, insertDistance=INSERT_DIST, _perturbation=_perturbation_func, _force=_general_force, particleKwargs= particleParams)

    cl.insert()

    # create sim using clump
    sim = s.growth_sim(cl, threeD=False, saveToFile=True)

    sim.draw_clump(showParticles=True, colour='olive')
    sim.show()