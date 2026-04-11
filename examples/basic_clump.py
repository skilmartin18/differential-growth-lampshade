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

LOAD_PATH = '../saves/very_dense_symm/layers.pkl'
SAVE_ID = '-med-dense'


SHOW = False


## create clump
cl = particle.clump()

# create sim using clump
sim = s.growth_sim(cl, threeD=True, saveToFile=True)

# load clump from file
sim.load_clump_from_layers(LOAD_PATH)
sim.clump = sim.layers[len(sim.layers)-1]
sim.draw_clump(colour='olive')
sim.show()
