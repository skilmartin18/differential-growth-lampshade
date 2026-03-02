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
import sim

log = logging.getLogger("logger")


### WIDE VASE ###
#
#
# A quick script to generate a wide vase, could be good as a light sconce?
#
# Bottom half of these is always wack, so use some other software to bin those away lol
#
# The terrible meshing algorithm im using generates some actually nice patterns here (even more emergent design..)
#
#


### CLUMP PARAMS ###
ITERATIONS = 300
_perturbation_func = particle.perturbation(prob=0.5, magnitude=2)
_general_force = particle.force(rDistance=75, aDistance=100, k=40, neighbors=1, clamp = 1e999)

SIZE = 12
RADIUS = 40
INSERT_DIST = 15

MASS = 1
FRICTION = 0.5
LIFE = 100

particleParams = {
    'mass': MASS,
    'friction': FRICTION,
    'life': LIFE
}

growingParams = {
    "growthZoneIterations": 999,
    "life": 100,
    "growth": 1,
    "coverage": 0,
}

### RUN SIMULATION ###

log.setLevel(logging.INFO)

## create clump
cl = particle.clump(size=SIZE, radius=RADIUS, insertDistance=INSERT_DIST, _perturbation=_perturbation_func, _force=_general_force, particleKwargs= particleParams)

# create sim using clump
sim = sim.growth_sim(cl, threeD=True, growingKwargs=growingParams, saveToFile=True)

# load clump from file
sim.load_clump_from_layers("/home/skilmartin/code/differential-growth/saves/very_dense_symm/layers.pkl")

# smoothing is what makes this one wide
sim.run_vase_smoothing(60, sim.gaussian_smooth, smoothingKwargs={'points': 11}, displayIts=999)
sim.trim_vase(150)
sim.draw_vase()
sim.convert_layers_to_mesh(showPlotly=True)
sim.export_mesh_stl()

