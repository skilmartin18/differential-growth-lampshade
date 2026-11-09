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

LOAD_PATH = '../saves/offcenter_growth/layers.pkl'
SAVE_ID = '-med-dense'

SMOOTH_IT = 2
SMOOTH_PTS = 5

SHOW = True

_rBounds = [1,160,15]



## create clump
cl = particle.clump()

# create sim using clump
sim = s.growth_sim(cl, threeD=True, saveToFile=True)

# load clump from file
sim.load_clump_from_layers(LOAD_PATH)
sim.convert_to_polar()
sim.normalise_to_circle()
sim.polar_update_clump()

# create vase using a single layer from previous sim
sim.clump = sim.layers[len(sim.layers)-1]
sim.clump.insertDistance = 5
sim.clump.insert()
sim.layers[len(sim.layers)-1] = sim.clump
sim.run_vase_smoothing(SMOOTH_IT, sim.gaussian_smooth, smoothingKwargs={'points': SMOOTH_PTS}, displayIts=50)

_rBounds = [50,160,15]

sim.modulated_vase(
    layers=20,
    layer=len(sim.layers)-1,
    period=0,
    zInc=3,
    rBounds=_rBounds,
    taper=0.0,
    rotation=0,
    showProgress=True
    )

sim.convert_layers_to_mesh(showPlotly=True)
sim.export_mesh_stl(filename="../stls/offcenter_narrowers_circ"+SAVE_ID)



if SHOW:
    sim.show()
# sim.convert_layers_to_mesh(showPlotly=True)
# sim.export_mesh_stl(filename="../stls/prism_taper"+SAVE_ID)
