
import matplotlib.pyplot as plt
import particle
import logging
import numpy as np
from scipy import signal
from datetime import datetime
import os
import pickle
from copy import deepcopy
import open3d as o3d

log = logging.getLogger("logger")

class growth_sim():
    
    def __init__(self, clump: particle.clump | str, threeD=False, saveToFile = False, growingKwargs = None):
        
        
        # load clump from object or pickle path
        if type(clump) == str:
            self.clump = self.load_clump_from_layers(clump,0)
        if type(clump) == particle.clump:
            self.clump = clump
            
        # basic plot details
        
        # make it interactive to show animations
        plt.ion() 
        self.fig = plt.figure(figsize=(12, 8))
        
        # for 3d (mainly used when looking at vases)
        self.threeD = threeD
        if threeD:
            self.ax = self.fig.add_subplot(projection='3d')
        else:
            self.ax = self.fig.add_subplot()
        
        self.ax.set_autoscale_on(True)

        # pretty colours for when progress is being saved
        self.colourMap = plt.cm.viridis 
       
        # kwargs for making growth zones
        self.growingKwargs = growingKwargs
        
        self.layers = []
        self.polarClump = []
        self.stop = False
    
    ### clump ###
    
    def grow_clump(self, iterations, animate = True, showProgress = False, saveToFile = False):
        
        # grow the input clump for a certain number of iterations 
        for i in range(iterations):
            
            self.clump.run_sim_iteration(incrementZ=False)
            self.insert_growing_zone(i, **self.growingKwargs)
            
            if animate:
                self.draw_clump(showProgress=showProgress, colour=self.colourMap(i/iterations))
            
            if self.stop:
                break
        
        # save an object and a png at the end of the run
        if saveToFile:
            self.layers.append(deepcopy(self.clump))
            self.save_clump()
      
    def insert_growing_zone(self, iteration, **kwargs):
        
        # get kwargs with defaults
        iterations = kwargs.get("growthZoneIterations", 1e99)
        coverage = kwargs.get('coverage', 0)
        life = kwargs.get('life', 100)
        growth = kwargs.get('growth',1)
        
        # dont trigger on first it, trigger on subsequent milestones
        if (iteration+1)%iterations == 0:
            log.info("creating zone")
            self.clump.create_growing_zone(coverage, life, growth)
    
    def draw_clump(self, showProgress = False, showParticles = False, colour = 'red'):
        
        # keep old lines drawn
        if not showProgress:
            self.ax.cla()
        
        positions =  [particle.positionVec for particle in self.clump.particles]
        positions.append(positions[0])
        
        xPos = [position[0] for position in positions]
        yPos = [position[1] for position in positions]
        
        if self.threeD:
            zPos = [position[2] for position in positions]
        else:
            zPos = []
        
        if showParticles:
            marker = 'o'
        else:
            marker = ""
            
        self.ax.plot(xPos,yPos,zPos, marker=marker, color = colour)
        
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()
        
        pass
    
    ### vase ###
    
    def grow_vase(self, iterations, layers = 20, startingIt = 0, saveToFile = False):
        
        itsperlayer = int(np.ceil(iterations/layers))
        layer = 0
        
        # the same as growing a clump, but this time we increase Z values to make a 2.5D vase
        for i in range(iterations):
            
            self.clump.run_sim_iteration(incrementZ=False)
            self.insert_growing_zone(i, **self.growingKwargs)
            
            # dont always want a layer per sim iteration
            if i % itsperlayer == 0 and i > startingIt:
            
                # increment Z in all clump particles
                for particle in self.clump.particles:
                    particle.positionVec[2] = layer
                
                layer = layer + 1
                self.draw_clump(showProgress=True, colour=self.colourMap(i/iterations))
                
                # inside loop only save to object, dont write to file
                if saveToFile:
                    self.layers.append(deepcopy(self.clump))
            
            if self.stop:
                break
        
        if saveToFile:
            self.save_clump()
  
    def draw_vase(self, showProgress = False):
        
        # two kinds of show progress, we want to retain each 2d plot but dont want 
        # to retain multiple vases
        if not showProgress:
            self.ax.cla()
            
        for i,layer in enumerate(self.layers):
            self.clump = layer
            self.draw_clump(showProgress=True, colour=self.colourMap(i/len(self.layers)))
        
    
    def draw_vase_from_layer(self, layers = 20, layer = 0, showProgress = False):
        
        # two kinds of show progress, we want to retain each 2d plot but dont want 
        # to retain multiple vases
        if not showProgress:
            self.ax.cla()
            
        firstLayer = deepcopy(self.layers[layer])
        self.layers = []
        
        for i in range(layers):
            
            clump = firstLayer
    
            # increment Z in all clump particles
            for particle in clump.particles:
                particle.positionVec[2] = i
            
            self.layers.append(deepcopy(clump))
            self.clump = clump
            self.draw_clump(showProgress=True, colour=self.colourMap(i/len(self.layers)))
    
    ### smoothing ###
    
    def smooth_vase(self, smoothing_function, smoothingKwargs = {}):
        
        for layer in self.layers:
            self.clump = layer
            smoothing_function(**smoothingKwargs)
    
    def run_vase_smoothing(self, iterations, smoothingFunction, smoothingKwargs = {}, displayIts = 10):
        
        log.info(f"running vase smoothing for {iterations} its")
        
        for i in range(iterations):
            
            self.smooth_vase(smoothingFunction, smoothingKwargs)
            
            if i % displayIts == 0:
                self.draw_vase()
     
    
    def moving_average_smooth(self, **kwargs):
        
        points = kwargs.get("points", 3)
        
        log.info(f'running moving average smoothing with {points} points')
        
        # create kernel
        blurKernelMovingAverage = np.repeat(1/points,points)
        blurKernelMovingAverage = np.array([np.repeat(0,points),blurKernelMovingAverage,np.repeat(0,points)])
        log.debug(blurKernelMovingAverage)
        
        particlesPositions = []
        particlesPositions.append([particle.positionVec for particle in self.clump.particles])
        particlesPositions = np.array(particlesPositions[0])
        log.debug(particlesPositions)
        
        smooth = signal.convolve2d(particlesPositions, blurKernelMovingAverage.T, mode='same', boundary='wrap')
        
        for i, p in enumerate(self.clump.particles):
            
            p.positionVec = smooth[i]
    
    def gaussian_smooth(self, **kwargs):
        
        points = kwargs.get("points", 5)
        
        log.info(f'running gaussian smoothing kernel size of {points}')
        
        def gaussian(x, mu, sig):
            return (
                1.0 / (np.sqrt(2.0 * np.pi) * sig) * np.exp(-np.power((x - mu) / sig, 2.0) / 2)
            )

        mu = np.floor(points/2)
        sig = points/6
        xVals = np.linspace(0,points-1,points)
        gaussianKernel = [gaussian(xVal,mu,sig) for xVal in xVals]
        gaussianKernel = np.array([np.repeat(0,points),gaussianKernel,np.repeat(0,points)])
        log.debug(gaussianKernel)
        
        particlesPositions = []
        particlesPositions.append([particle.positionVec for particle in self.clump.particles])
        particlesPositions = np.array(particlesPositions[0])
        log.debug(particlesPositions)
        
        smooth = signal.convolve2d(particlesPositions, gaussianKernel.T, mode='same', boundary='wrap')
        
        for i, p in enumerate(self.clump.particles):
            
            p.positionVec = smooth[i]
        
    
    def run_clump_smoothing(self, iterations, smoothing_function, smoothingKwargs, displayIts = 10):
        
        for i in range(iterations):
            log.info(f"running smoothing iteration {i}")
            smoothing_function(**smoothingKwargs)
            self.draw_clump()

    def convert_to_polar(self):
        
        particlesPositions = []
        particlesPositions.append([particle.positionVec for particle in self.clump.particles])
        
        # this is terrible and unreadable but it works
        particlesPositions[0].append(particlesPositions[0][0])
        particlesPositions = np.array(particlesPositions[0])
        
        
        X = np.copy(particlesPositions[:, 0])
        Y = np.copy(particlesPositions[:,1])
    
        R = np.sqrt(X**2+Y**2)
        T = np.arctan2(Y,X)
        
        self.polarClump = [T,R]
        return [T,R]
    
    def convert_to_cartesian(self):
 
        x = self.polarClump[1] * np.cos(self.polarClump[0])
        y = self.polarClump[1] * np.sin(self.polarClump[0])
        return [x,y]
    
    def polar_update_clump(self):
        
        [x,y] = self.convert_to_cartesian()
        
        for i,particle in enumerate(self.clump.particles):
            
            particle.positionVec[0] = x[i]
            particle.positionVec[1] = y[i]
    
    def normalise_to_circle(self, min = 80, max = 160):
        
        # ingest the clump in POLAR format
        
        # we normalise the values to [0 1]
        minMag = np.min(self.polarClump[1])
        maxMag = np.max(self.polarClump[1])
        
        magNorms = []
        for mag in self.polarClump[1]:
            mag = (mag-minMag)/(maxMag-minMag)
            magNorms.append(mag)
        
        magNorms = np.array(magNorms)
        
        # we map that to our new range [min max]
        newRange = max-min
        magNorms = magNorms*newRange
        magNorms = magNorms+min
        
      
        self.polarClump = [self.polarClump[0], magNorms]
    
    def draw_polar(self):   
        # figures for polar normalisation
        self.polarfig = plt.figure(figsize=(6, 6))
        self.polarax = self.polarfig.add_subplot(111, projection='polar') 
        self.polarax.plot(self.polarClump[0], self.polarClump[1], color='b', linewidth=2)
        
    ### utility ###

    def init_interactive_sim(self):
        
        # make room for sliders
        self.fig.subplots_adjust(left=0.25, bottom=0.25)

        # SIZE
        axsize = self.fig.add_axes([0.15, 0.05, 0.3, 0.03])
        self.size_slider = plt.Slider(
            ax=axsize,
            label='Initial Size',
            valmin=4,
            valmax=60,
            valinit=8,
            valstep=1
        )

        # RADIUS
        axradius = self.fig.add_axes([0.15, 0.1, 0.3, 0.03])
        self.radius_slider = plt.Slider(
            ax=axradius,
            label='Initial Circle Radius',
            valmin=1,
            valmax=60,
            valinit=20,
            valstep=1
        )
        
        # INSERT
        axinsert = self.fig.add_axes([0.15, 0.15, 0.3, 0.03])
        self.insert_slider = plt.Slider(
            ax=axinsert,
            label='Insert Distance',
            valmin=5,
            valmax=60,
            valinit=15,
            valstep=1
        )
        
        # FRICTION
        axfriction = self.fig.add_axes([0.15, 0.2, 0.3, 0.03])
        self.friction_slider = plt.Slider(
            ax=axfriction,
            label='Friction',
            valmin=0,
            valmax=1,
            valinit=0.7,
        )
        
        # iterations
        axiterations = self.fig.add_axes([0.6, 0.05, 0.3, 0.03])
        self.iterations_slider = plt.Slider(
            ax=axiterations,
            label='iterations',
            valmin=1,
            valmax=1000,
            valinit=100,
            valstep=1,
        )
        
        # life
        axlife = self.fig.add_axes([0.6, 0.1, 0.3, 0.03])
        self.life_slider = plt.Slider(
            ax=axlife,
            label='life',
            valmin=0,
            valmax=1000,
            valinit=100,
            valstep=1
        )
        
        # mass
        axmass = self.fig.add_axes([0.6, 0.15, 0.3, 0.03])
        self.mass_slider = plt.Slider(
            ax=axmass,
            label='mass',
            valmin=0,
            valmax=5,
            valinit=1,
        )

        # Whenever we update the slider, we reset to simulation, and run the iterations again
        def update(val):
            self.stop = True
            self.clump.stop = True
            del self.clump
            self.clump = particle.clump(size=self.size_slider.val, radius=self.radius_slider.val,insertDistance=self.insert_slider.val,_perturbation=_perturbation_func,_force=_general_force, particleKwargs= particleParams)
            self.stop = False
            self.grow_clump(self.iterations_slider.val)
            self.fig.canvas.draw_idle()


        self.size_slider.on_changed(update)
        self.radius_slider.on_changed(update)
        self.insert_slider.on_changed(update)
        self.friction_slider.on_changed(update)
        self.iterations_slider.on_changed(update)
        self.life_slider.on_changed(update)
        self.mass_slider.on_changed(update)
        
    def show(self):
        plt.show(block=True)
    
    def save_clump(self, path = None, filename = None):
        
        if path == None:
            now = datetime.now()
            timestamp_str = now.strftime("%Y-%m-%d_%H-%M-%S")
            path = "./saves/"+timestamp_str+"/"
        
        os.makedirs(path)

        # save the actual layer details
        with open(path+"layers.pkl", 'wb') as handle:
            pickle.dump(self.layers, handle, protocol=pickle.HIGHEST_PROTOCOL)

        # save a jpeg of the matplotlib window? overwrite so only one per run
        self.fig.savefig(path+"preview.png")
        
        log.info(f"sim saved to {path}")
        
    def load_clump_from_layers(self, path, layer: str | int = 'last'):
        
        with open(path, 'rb') as handle:
            layers = pickle.load(handle)
            
        if layer == 'last':
            layer = len(layers)-1
        
        
        
        self.layers = layers
        self.clump = layers[layer]
        
        log.info(f"loaded layers with {len(layers)} layers")
        log.info(f"loaded layer {layer} with {len(self.clump.particles)} vertices")
        
        return layers[layer]
    
if __name__ == "__main__":
    
    
    ### CLUMP PARAMS ###
    ITERATIONS = 500
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
    sim = growth_sim(cl, threeD=True, growingKwargs=growingParams, saveToFile=True)
    
    # create sim using path to pickle
    # sim = growth_sim("/home/skilmartin/code/differential-growth/saves/saves/example_save/layers.pkl", threeD=True, growingKwargs=growingParams)
    
    # simulate clump growth
    # sim.grow_clump(ITERATIONS, showProgress=True)
    
    # simulate vase growth
    # sim.grow_vase(ITERATIONS, layers=50,saveToFile=True)
    
    # load clump from file
    sim.load_clump_from_layers("/home/skilmartin/code/differential-growth/saves/example_save/layers.pkl",layer=48)
    # sim.draw_vase()
    
    # sim.run_vase_smoothing(100, sim.gaussian_smooth, smoothingKwargs={'points': 11}, displayIts=20)
    
    # open3d testing
    
    # 1. Create a NumPy array of 3D points (N, 3 shape)
    # This example creates 100 random points for demonstration
    
    # WILL PROBABLY NEED TO SUBDIVIDE SO ALL LAYERS HAVE THE SAME AMOUNT OF VERTICES
    
    # find the amount of additional vertices required
    
    # order the current list of vertices pairs, by the magnitude between them, highest to lowest
    # i.e for each pair of vertices, calculate the magnitude between them
    # label the pairs 
    # sort the magnitude (and labels) highest to lowest
    # starting with the pair with greatest magnitude, insert a point between them
    # run through until the magnitude labels list is exhausted 
    # repeat until correct number of vertices
    
    # repeat for all layers
    
    # final layer is always the largest
    maxVertices = len(sim.layers[-1].particles)
    
    print(maxVertices)
    
    allPositions = []
    for i,layer in enumerate(sim.layers):
        
        positions =  np.array([particle.positionVec for particle in layer.particles])
        newPositions = deepcopy(positions)
        
        # shift all positions one to the left (vert x+1 now in position x)
        positionsNextVert = np.roll(positions,-1, axis=0)
        
        # difference between next vertice and current one normed
        # gives an array of magnitudes, where each index corresponds to a colinear pair of vertices
        # i.e index 1 has the distance between 1->2, index 2 has 2->3
        nextVertMag = np.linalg.norm(positionsNextVert - positions, axis=1)

        # sorted indices
        sortedVertMagI = np.flip(np.argsort(nextVertMag))
        
        print(f"num pos before {len(newPositions)}")
        for pair in sortedVertMagI:
            
            print("inserting")
            newParticlePos = positions[pair] + positionsNextVert[pair]*0.5
            newPositions = np.insert(newPositions, pair, newParticlePos, axis=0)
            # we now have the pair of points with the highest distance between them
            # pair X means vertex X->X+1
            # we may have to insert into a new array as we are iterating through the same array that we are editing
            pass
        
        print(f"num pos after {len(newPositions)}")
        
        break
        
        # # print(positions)
        # allPositions.append(positions)
    print(newPositions)
    xPos = [position[0] for position in newPositions]
    yPos = [position[1] for position in newPositions]
    
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot()
        
    ax.set_autoscale_on(True)
    ax.plot(xPos,yPos, color = 'red')
    
    print(len(xPos))
    print(len(yPos))
    plt.ioff()
    plt.show()
    
    xPos = [position[0] for position in positions]
    yPos = [position[1] for position in positions]
    
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot()
        
    ax.set_autoscale_on(True)
    ax.plot(xPos,yPos, color = 'red')
    
    print(len(xPos))
    print(len(yPos))
    plt.ioff()
    plt.show()
    
    # 2. Create a PointCloud object
    # pcd = o3d.geometry.PointCloud()
    # # 3. Assign the NumPy array to the point cloud
    # pcd.points = o3d.utility.Vector3dVector(allPositions)
    
    # pcd.estimate_normals()

    
    # radii = [1,2,3,4,5,6,7,8,9,10,15,20,25]
    # rec_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_ball_pivoting(
    # pcd, o3d.utility.DoubleVector(radii))
    # # Optional: Visualize the point cloud
    # o3d.visualization.draw_plotly([rec_mesh])
        
    
    
    
    
    
    
    # sim.convert_to_polar()
    # sim.normalise_to_circle()
    # sim.polar_update_clump()
    # sim.draw_clump()
    # # create vase using a single layer from previous sim
    # sim.draw_vase_from_layer(layer=48)
    # sim.run_vase_smoothing(20, sim.gaussian_smooth, smoothingKwargs={'points': 7}, displayIts=10)
    # sim.draw_polar()
    # sim.show()
    