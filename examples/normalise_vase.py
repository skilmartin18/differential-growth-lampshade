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

SMOOTH_IT = 2
SMOOTH_PTS = 5

SHOW = True

_rBounds = [60  ,160,15]


## STRAIGHT WALLS
if False:
    ## create clump
    cl = particle.clump()

    # create sim using clump
    sim = s.growth_sim(cl, threeD=True, saveToFile=True)

    # load clump from file
    sim.load_clump_from_layers(LOAD_PATH)
   

    # create vase using a single layer from previous sim

    sim.run_vase_smoothing(SMOOTH_IT, sim.gaussian_smooth, smoothingKwargs={'points': SMOOTH_PTS}, displayIts=50)

    sim.modulated_vase(
        layers=100,
        layer=len(sim.layers)-1,
        period=0,
        zInc=3,
        rBounds=_rBounds,
        taper=0.0,
        rotation=0,
        )

    if SHOW:
        sim.show()
    sim.convert_layers_to_mesh(showPlotly=True)
    sim.export_mesh_stl(filename="../stls/prism"+SAVE_ID)

## PLAIN TAPER WALLS
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

    sim.modulated_vase(
        layers=100,
        layer=len(sim.layers)-1,
        period=0,
        zInc=3,
        rBounds=_rBounds,
        taper=0.7,
        rotation=0,
        )

    if SHOW:
        sim.show()
    sim.convert_layers_to_mesh(showPlotly=True)
    sim.export_mesh_stl(filename="../stls/prism_taper"+SAVE_ID)

## Varying width WALLS
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

    sim.modulated_vase(
        layers=100,
        layer=len(sim.layers)-1,
        period=2.5,
        zInc=5,
        rBounds=_rBounds,
        taper=0.0,
        rotation=0,
        )

    if SHOW:
        sim.show()
    sim.convert_layers_to_mesh(showPlotly=True)
    sim.export_mesh_stl(filename="../stls/prism_width"+SAVE_ID)

## Varying width WALLS and rotation
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

    sim.modulated_vase(
        layers=100,
        layer=len(sim.layers)-1,
        period=1.2,
        zInc=5,
        rBounds=_rBounds,
        taper=0.2,
        rotation=90,
        )

    if SHOW:
        sim.show()
    sim.convert_layers_to_mesh(showPlotly=True)
    sim.export_mesh_stl(filename="../stls/prism_width_rotation"+SAVE_ID)

## Varying width WALLS and rotation and taper
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
        layers=100,
        layer=len(sim.layers)-1,
        period=0,
        zInc=5,
        rBounds=_rBounds,
        taper=0.8,
        rotation=80,
        )

    if SHOW:
        sim.show()
    sim.convert_layers_to_mesh(showPlotly=True)
    sim.export_mesh_stl(filename="../stls/prism_width_rotation_taper"+SAVE_ID)

## Varying width WALLS and rotation and custom taper function
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
        if i < 30:
        #    print(np.sin( (i/29 * 0.6 * 2 *3.14) + 2.8*3.14/2) * 80)
           return np.sin( (i/30 * 0.6 * 2 *3.14) + 2.8*3.14/2) * 80
        elif i >= 30 and i < 50:
            return  76.18623736083171 - 0.9*(i-30)
        else:
            return  np.sin( (0.6 * 2 *3.14) + 2.8*3.14/2) * 80 + 0.2*(i-50) - 0.9*30
       
    
    sim.modulated_vase(
        layers=200,
        layer=len(sim.layers)-1,
        period=0,
        zInc=5,
        rBounds=_rBounds,
        taper=taper_func,
        rotation=80,
        )

    if SHOW:
        sim.show()
    sim.convert_layers_to_mesh(showPlotly=True)
    sim.export_mesh_stl(filename="../stls/prism_width_rotation_custom_taper"+SAVE_ID)