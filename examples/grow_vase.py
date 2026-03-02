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
import logging
sys.path.append('../')
import particle
import sim

### CLUMP PARAMS ###

# this will be a basic clump

if False:
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
    log = logging.getLogger("logger")
    log.setLevel(logging.INFO)

    ## create clump
    cl = particle.clump(size=SIZE, radius=RADIUS, insertDistance=INSERT_DIST, _perturbation=_perturbation_func, _force=_general_force, particleKwargs= particleParams)

    # create sim using clump
    sim = sim.growth_sim(cl, threeD=True, growingKwargs=growingParams, saveToFile=True)

    # simulate vase growth
    sim.grow_vase(ITERATIONS, layers=ITERATIONS,saveToFile=True, draw=False)
    sim.draw_vase(showProgress=True)
    sim.show()

if False: 

    ### CLUMP PARAMS ###

    # faster growing
    ITERATIONS = 200
    _perturbation_func = particle.perturbation(prob=0.5, magnitude=2)
    _general_force = particle.force(rDistance=75, aDistance=100, k=40, neighbors=1, clamp = 1e999)

    SIZE = 12
    RADIUS = 60
    INSERT_DIST = 15

    MASS = 1
    FRICTION = 0.65
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
    log = logging.getLogger("logger")
    log.setLevel(logging.INFO)

    ## create clump
    cl = particle.clump(size=SIZE, radius=RADIUS, insertDistance=INSERT_DIST, _perturbation=_perturbation_func, _force=_general_force, particleKwargs= particleParams)

    # create sim using clump
    sim = sim.growth_sim(cl, threeD=True, growingKwargs=growingParams, saveToFile=True)

    # simulate vase growth
    sim.grow_vase(ITERATIONS, layers=ITERATIONS,saveToFile=True, draw=True)
    sim.draw_vase(showProgress=True)
    sim.show()
       
if False: 

    ### CLUMP PARAMS ###

    # faster growing
    ITERATIONS = 200
    _perturbation_func = particle.perturbation(prob=0.5, magnitude=2)
    _general_force = particle.force(rDistance=75, aDistance=100, k=40, neighbors=1, clamp = 1e999)

    SIZE = 12
    RADIUS = 30
    INSERT_DIST = 15

    MASS = 1
    FRICTION = 0.65
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
    log = logging.getLogger("logger")
    log.setLevel(logging.INFO)

    ## create clump
    cl = particle.clump(size=SIZE, radius=RADIUS, insertDistance=INSERT_DIST, _perturbation=_perturbation_func, _force=_general_force, particleKwargs= particleParams)

    # create sim using clump
    sim = sim.growth_sim(cl, threeD=True, growingKwargs=growingParams, saveToFile=True)

    # simulate vase growth
    sim.grow_vase(ITERATIONS, layers=ITERATIONS,saveToFile=True, draw=True)
    sim.draw_vase(showProgress=True)
    sim.show()
        
if False: 

    ### CLUMP PARAMS ###

    # faster growing
    ITERATIONS = 200
    _perturbation_func = particle.perturbation(prob=0.5, magnitude=2)
    _general_force = particle.force(rDistance=75, aDistance=100, k=40, neighbors=1, clamp = 1e999)

    SIZE = 6
    RADIUS = 30
    INSERT_DIST = 15

    MASS = 1
    FRICTION = 0.6
    LIFE = 30

    particleParams = {
        'mass': MASS,
        'friction': FRICTION,
        'life': LIFE
    }

    growingParams = {
        "growthZoneIterations": 30,
        "life": 100,
        "growth": 1,
        "coverage": 0.3,
    }

    ### RUN SIMULATION ###
    log = logging.getLogger("logger")
    log.setLevel(logging.INFO)

    ## create clump
    cl = particle.clump(size=SIZE, radius=RADIUS, insertDistance=INSERT_DIST, _perturbation=_perturbation_func, _force=_general_force, particleKwargs= particleParams)

    # create sim using clump
    sim = sim.growth_sim(cl, threeD=True, growingKwargs=growingParams, saveToFile=True)

    # simulate vase growth
    sim.grow_vase(ITERATIONS, layers=ITERATIONS,saveToFile=True, draw=True)
    sim.draw_vase(showProgress=True)
    sim.show()
    
       
if True: 

    ### CLUMP PARAMS ###

    # faster growing
    ITERATIONS = 200
    _perturbation_func = particle.perturbation(prob=0.5, magnitude=2)
    _general_force = particle.force(rDistance=75, aDistance=100, k=40, neighbors=1, clamp = 1e999)

    SIZE = 6
    RADIUS = 30
    INSERT_DIST = 15

    MASS = 1
    FRICTION = 0.6
    LIFE = 20

    particleParams = {
        'mass': MASS,
        'friction': FRICTION,
        'life': LIFE
    }

    growingParams = {
        "growthZoneIterations": 30,
        "life": 100,
        "growth": 1,
        "coverage": 0.3,
    }

    ### RUN SIMULATION ###
    log = logging.getLogger("logger")
    log.setLevel(logging.INFO)

    ## create clump
    cl = particle.clump(size=SIZE, radius=RADIUS, insertDistance=INSERT_DIST, _perturbation=_perturbation_func, _force=_general_force, particleKwargs= particleParams)

    # create sim using clump
    sim = sim.growth_sim(cl, threeD=True, growingKwargs=growingParams, saveToFile=True)

    # simulate vase growth
    sim.grow_vase(ITERATIONS, layers=ITERATIONS,saveToFile=True, draw=True)
    sim.draw_vase(showProgress=True)
    sim.show()