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

SMOOTH_IT = 5
SMOOTH_PTS = 5

SHOW = True

_rBounds = [60  ,160,15]


## Varying width WALLS and rotation and custom taper function for lamp base
if False:
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

    sim.run_vase_smoothing(SMOOTH_IT, sim.gaussian_smooth, smoothingKwargs={'points': SMOOTH_PTS}, displayIts=50)


    def taper_func(i):
        print(i)
        if i < 25:
        #    print(np.sin( (i/29 * 0.6 * 2 *3.14) + 2.8*3.14/2) * 80)
            
            return np.sin( (i/25 * 0.8 * 2 *3.14) + 2.8*3.14/2) * 40 + (i-25)*2
        else:
            return np.sin( (25/25 * 0.8 * 2 *3.14) + 2.8*3.14/2) * 40 + (i-25)*0.7
       
    
    sim.modulated_vase(
        layers=130,
        layer=len(sim.layers)-1,
        period=0,
        zInc=8,
        rBounds=_rBounds,
        taper=taper_func,
        rotation=80,
        )

    if SHOW:
        sim.show()
    sim.convert_layers_to_mesh(showPlotly=True)
    sim.export_mesh_stl(filename="../stls/normallookinglamp_base"+SAVE_ID)


_rBounds = [60  ,160,15]
## Varying width WALLS and rotation and custom taper function for lamp base
if True:
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

    sim.run_vase_smoothing(SMOOTH_IT, sim.gaussian_smooth, smoothingKwargs={'points': SMOOTH_PTS}, displayIts=50)

    
    sim.modulated_vase(
        layers=130,
        layer=len(sim.layers)-1,
        period=0,
        zInc=8,
        rBounds=_rBounds,
        taper=0.0,
        rotation=60,
        )

    if SHOW:
        sim.show()
    sim.convert_layers_to_mesh(showPlotly=True)
    sim.export_mesh_stl(filename="../stls/normallookinglamp_shade"+SAVE_ID)