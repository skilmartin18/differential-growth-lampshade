
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
from collections.abc import Callable

log = logging.getLogger("differential-growth")

"""growth_sim is a utility module used to help generate lampshades/vases from simulated clumps
It has slowly grown as required, and unfortunately is one of the worst designed pieces of software
i have ever written.
"""

class growth_sim():
    
    def __init__(self, clump: particle.clump | str, threeD: bool =False, saveToFile: bool = False, growingKwargs: dict = None):
        """initialise a vase simulation

        Args:
            clump (particle.clump | str): either pass in a clump to grow, or the path to a pickled layers file
            threeD (bool, optional): changes the plotting mode to 3d. Defaults to False.
            saveToFile (bool, optional): save the resultant vase to a file. Defaults to False.
            growingKwargs (dict, optional): growthZone kwargs to be passed to clump. Defaults to None.
        """
        
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
        self.saveToFile = saveToFile
        self.layers = []
        self.polarClump = []
        self.stop = False
        self.fignum = 0
    ### clump ###
    
    def grow_clump(self, iterations: int, animate: bool = True, showProgress: bool = False):
        """grows self.clump for a set number of iterations

        Args:
            iterations (int): number of iterations
            animate (bool, optional): should the clump be drawn each step. Defaults to True.
            showProgress (bool, optional): should canvas be wiped or retained per iteration. Defaults to False.
        """
        for i in range(iterations):
            self.fignum = i
            self.clump.run_sim_iteration()
            self.insert_growing_zone(i, **self.growingKwargs)
            
            if animate:
                self.draw_clump(showProgress=showProgress, colour=self.colourMap(i/iterations))
            
            if self.stop:
                break
        
        # save an object and a png at the end of the run
        if self.saveToFile:
            self.layers.append(deepcopy(self.clump))
            self.save_clump()
      
    def insert_growing_zone(self, iteration: int, **kwargs):
        """insert a growing zone in the clump based on growingKwargs

        Args:
            iteration (int): current sim iteration
        """
        # get kwargs with defaults
        # how many iterations before a zone is added
        iterations = kwargs.get("growthZoneIterations", 1e99)
        coverage = kwargs.get('coverage', 0)
        life = kwargs.get('life', 100)
        growth = kwargs.get('growth',1)
        
        # dont trigger on first it, trigger on subsequent milestones
        if (iteration+1)%iterations == 0:
            log.info("creating zone")
            self.clump.create_growing_zone(coverage, life, growth)
    
    def draw_clump(self, showProgress: bool = False, showParticles: bool = False, colour: str = 'red'):
        """draw a 2d clump

        Args:
            showProgress (bool, optional): should canvas be wiped. Defaults to False.
            showParticles (bool, optional): show particles as well as individual lines. Defaults to False.
            colour (str, optional): matplotlib compatible colour. Defaults to 'red'.
        """
        
        # keep old lines drawn?
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
        self.fig.savefig(f"./saves/{self.fignum}.png")
        self.fig.canvas.flush_events()
        
    
    ### vase ###
    
    def grow_vase(self, iterations: int, layers: int = 20, startingIt: int = 0, draw: bool = True):
        """grows a vase (a 2.5D clump)

        Args:
            iterations (int): how many iterations to grow for
            layers (int, optional): how many layers tall to generate. Defaults to 20.
            startingIt (int, optional): how many sim iterations to run before starting vase creation. Defaults to 0.
            draw (bool, optional): draw the vase?. Defaults to True.
        """
        
        itsperlayer = int(np.ceil(iterations/layers))
        layer = 0
        
        # the same as growing a clump, but this time we increase Z values to make a 2.5D vase
        for i in range(iterations):
            
            self.clump.run_sim_iteration()
            self.insert_growing_zone(i, **self.growingKwargs)
            
            # dont always want a layer per sim iteration
            if i % itsperlayer == 0 and i > startingIt:
            
                # increment Z in all clump particles
                for particle in self.clump.particles:
                    particle.positionVec[2] = layer
                
                layer = layer + 1
                
                if draw:
                    self.draw_clump(showProgress=True, colour=self.colourMap(i/iterations))
                
                # inside loop only save to object, dont write to file
                if self.saveToFile:
                    self.layers.append(deepcopy(self.clump))
            
            if self.stop:
                break
        
        if self.saveToFile:
            self.save_clump()
  
    def draw_vase(self, showProgress: bool = False):
        """draw a whole '2.75D' vase

        Args:
            showProgress (bool, optional): show progress... Defaults to False.
        """
        
        # two kinds of show progress, we want to retain each 2d plot but dont want 
        # to retain multiple vases
        
        # this draws a 2.75d vase... although not using 3D differential growth, we use 
        # new simulations per layer giving a tree/coral like view
        
        if not showProgress:
            self.ax.cla()
            
        for i,layer in enumerate(self.layers):
            self.clump = layer
            self.draw_clump(showProgress=True, colour=self.colourMap(i/len(self.layers)))
        
    
    def draw_vase_from_layer(self, layers: int = 20, layer: int = 0, showProgress: bool = False):
        """Draws a 2.5D vase (a 3d vase which is extruded from a single layer)

        Args:
            layers (int, optional): how many layers high to extrude. Defaults to 20.
            layer (int, optional): what layer to extrude. Defaults to 0.
            showProgress (bool, optional): show progress... Defaults to False.
        """
        
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
    
    def modulated_vase(self,layers = 20, layer = 0, rBounds = (80,160,15), period=1, showProgress=False, zInc = 1, taper: float | Callable = 0, rotation = 0):
        """generate a modulated vase- this is the core function for generating cool vases from clump perimeters. 

        Args:
            layers (int, optional): how many layers tall to make the vase. Defaults to 20.
            layer (int, optional): which layer to use to create the vase. Defaults to 0.
            rBounds (tuple, optional): for normalising the clump to a circle- tuple contains (minimum radius, maximum radius, mod)
            where mod is the magnitude of the sinusoidal radius variations. Defaults to (80,160,15).
            period (int, optional): period of sinusoidal variations. Defaults to 1.
            showProgress (bool, optional): False wipes figure after each iteration. Defaults to False.
            zInc (int, optional): Z-axis distance between layers. Defaults to 1.
            taper (float | Callable, optional): float for linear per-iteration reduction in radius, else function to determine 
            taper amount based on iteration. Defaults to 0.
            rotation (int, optional): angular rotation (degrees) between consecutive vase layers. Defaults to 0.
        """
        # two kinds of show progress, we want to retain each 2d plot but dont want 
        # to retain multiple vases
        if not showProgress:
            self.ax.cla()
            
        firstLayer = deepcopy(self.layers[layer])
        self.layers = []
        
        for i in range(layers):
            
            self.fignum = i
            rangeFactor = np.sin( (i/layers * period * 2 *3.14) ) * rBounds[2]
      
            # for each layer we get a clump (functions operate on self so make ach layer hte main clump)
            self.clump = deepcopy(firstLayer)
            
            # we then convert it to polar
            self.convert_to_polar()
            
            if type(taper) == float:
                t = taper*i
            elif callable(taper):
                t = taper(i)
                
            # normalise, with bounds alternating based on a sin wave
            self.normalise_to_circle(min=rBounds[0]+rangeFactor-t,max=rBounds[1]+rangeFactor-t)
            self.rotate_layer(np.deg2rad( (i/layers)*rotation))
            self.polar_update_clump()
    
            # increment Z in all clump particles
            for particle in self.clump.particles:
                particle.positionVec[2] = i*zInc
            
            self.layers.append(deepcopy(self.clump))
            self.draw_clump(showProgress=showProgress, colour=self.colourMap(i/layers))
        
    def rotate_layer(self,radians: float = 0):
        """rotate the current sims polar clump

        Args:
            radians (float, optional): radians to rotate by. Defaults to 0.
        """
        for i in range(len(self.polarClump[0])):
            self.polarClump[0][i] = (self.polarClump[0][i]+radians)
   
    ### smoothing ###
    
    def smooth_vase(self, smoothing_function: Callable, smoothingKwargs = {}):
        """apply smoothing to all layers of a vase

        Args:
            smoothing_function (Callable): some smoothing function
            smoothingKwargs (dict, optional): kawrgs for smoothing. Defaults to {}.
        """
        for layer in self.layers:
            self.clump = layer
            smoothing_function(**smoothingKwargs)
    
    def run_vase_smoothing(self, iterations: int, smoothingFunction: Callable, smoothingKwargs = {}, displayIts: int = 10):
        """run multiple iterations of vase smoothin and display

        Args:
            iterations (int): how many iterations to smooth for
            smoothingFunction (Callable): smoothing function
            smoothingKwargs (dict, optional): smoothing kwargs. Defaults to {}.
            displayIts (int, optional): draw vase every displayIts. Defaults to 10.
        """
        log.info(f"running vase smoothing for {iterations} its")
        
        for i in range(iterations):
            
            self.smooth_vase(smoothingFunction, smoothingKwargs)
            
            if i+1 % displayIts == 0:
                self.draw_vase()
     
    
    def moving_average_smooth(self, **kwargs):
        """smoothing function using an X point MA
        """
        
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
        """gaussian smooth using an X point gaussian dist.
        """
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
        
    
    def run_clump_smoothing(self, iterations: int, smoothing_function: Callable, smoothingKwargs = {}, displayIts = 10):
        """run iteratons of smoothing on a clump

        Args:
            iterations (int): how many smoothing iterations to run  
            smoothing_function (Callable): smoothing function
            smoothingKwargs (dict, optional): smoothing kwargs. Defaults to {}.
            displayIts (int, optional): unused lol. Defaults to 10.
        """
        for i in range(iterations):
            log.info(f"running smoothing iteration {i}")
            smoothing_function(**smoothingKwargs)
            self.draw_clump()

    def convert_to_polar(self):
        """converts the current clump particles into polar representation (angle, mag)
        """
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
        """converts current polar clump into cartensian form
        """
        x = self.polarClump[1] * np.cos(self.polarClump[0])
        y = self.polarClump[1] * np.sin(self.polarClump[0])
        return [x,y]
    
    def polar_update_clump(self):
        """uses convert_to_cartesian to convert current polar clump to cartesian 
        andd overwrites current cartesian clump w it
        """
        [x,y] = self.convert_to_cartesian()
        
        for i,particle in enumerate(self.clump.particles):
            
            particle.positionVec[0] = x[i]
            particle.positionVec[1] = y[i]
    
    def normalise_to_circle(self, min: float = 80, max: float = 160):
        """normalises the coordinates of the current polar clump to sit linearly between a new min and max 
        bound. 

        Args:
            min (float, optional): minimum magnitude. Defaults to 80.
            max (float, optional): maximum magnitude. Defaults to 160.
        """
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
        """draws figure from polar clump on a polar axis
        """
        # figures for polar normalisation
        self.polarfig = plt.figure(figsize=(6, 6))
        self.polarax = self.polarfig.add_subplot(111, projection='polar') 
        self.polarax.plot(self.polarClump[0], self.polarClump[1], color='b', linewidth=2)
        
    ### utility ###

    def init_interactive_sim(self):
        """never used 
        """
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
        """show AND block
        """
        plt.show(block=True)
    
    def save_clump(self, path: str = None, filename: str = None): 
        """save the raw clump object and a png to a file

        Args:
            path (str, optional): path to saved file. Defaults to ./saves/{timestamp}/.
            filename (str, optional): filename for saved file. Defaults to layers.pkl and preview.png.
        """
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
        
    def load_clump_from_layers(self, path: str, layer: str | int = 'last') -> particle.clump:
        """loads a clump from a layer file 

        Args:
            path (str): path to layers pkl
            layer (str | int, optional): index of layer. Defaults to 'last'.

        Returns:
            particle.clump: _description_
        """
        with open(path, 'rb') as handle:
            layers = pickle.load(handle)
            
        if layer == 'last':
            layer = len(layers)-1
        
        
        
        self.layers = layers
        self.clump = layers[layer]
        
        log.info(f"loaded layers with {len(layers)} layers")
        log.info(f"loaded layer {layer} with {len(self.clump.particles)} vertices")
        
        return layers[layer]
    
    def trim_vase(self,layersToTrim: int):
        """trims a vase from a point until end of vase

        Args:
            layersToTrim (int): trims vase from this until end
        """
        self.layers = self.layers[layersToTrim:]
    
  
    def convert_layers_to_mesh(self,meshTopBottom = False, showPlotly=False):
        """converts a set of vase layers into a mesh and displays it in the browser

        Args:
            meshTopBottom (bool, optional): who knows. Defaults to False.
            showPlotly (bool, optional): display output in browser. Defaults to False.

        Returns:
            o3d.geometry.TriangleMesh: mesh
        """
        # this is possibly THE WORST way to do this.
        # i have no 3d modelling experience, so this is a brute force method of creating a mesh
        
        # final layer is always the largest (this is mathematically true... not a guess)
        maxVertices = len(self.layers[-1].particles)
        log.info(f"all layers need to be remeshed to {maxVertices} vertices")
        
        newLayers = []
        for i,layer in enumerate(self.layers):
            
            positions =  np.array([particle.positionVec for particle in layer.particles])
            log.info(f"adding {maxVertices-len(positions)} vertices to layer {i}")
            
            while len(positions) < maxVertices:
                
                # shift all positions one to the right (vert x-1 now in position x)
                # (particles are orderer anti clockwise?)
                positionsNextVert = np.roll(positions,-1, axis=0)

                # difference between next vertice and current one normed
                # gives an array of magnitudes, where each index corresponds to a colinear pair of vertices
                # i.e index 1 has the distance between 1->2, index 2 has 2->3
                nextVertMag = np.linalg.norm(positionsNextVert - positions, axis=1)
            
                # get pair with highest magnitude
                # previous code used the whole list and iterated though it, but adding data to a list while adding to it 
                # causes issues that are easy to fix with brute force instead of thinking
                furthestPair = np.flip(np.argsort(nextVertMag))[0]
            
                # get a new vertex halfway between current and neighbor
                diff = positions[(furthestPair+1)%len(positions)] - positions[furthestPair]
                newParticlePos = positions[furthestPair] + diff*0.5
                
                positions = np.insert(positions, furthestPair+1, newParticlePos, axis=0)
                # we now have the pair of points with the highest distance between them
                # pair X means vertex X->X+1
                # we may have to insert into a new array as we are iterating through the same array that we are editing
            
            newLayers.append(positions)
            
        # convert array with seperate layers into one long array
        allLayers = np.concatenate(newLayers)
        
        # to index a point at vertex i in layer j we use this formula (modulo so we can wrap in x axis)
        # allLayers[ (i%maxVertices) + j*maxVertices ]
        def ij_to_index(i,j):
            return (i%maxVertices) + j*maxVertices
        
        triangles = []
        # iterate through all layers, starting at the bottom, skip the last layer as it does not connect to anything
        for j in range(len(newLayers)-1):
            
            for i in range(maxVertices):
                
                # first triangle goes from
                # *  *  *    
                # |\
                # | \
                # *--*  *    
                # i,j -> i+1,j -> i,j+1
                # second triangle goes from
                # *--*  *    
                # |\ |
                # | \|
                # *--*  * 
                # i,j+1 -> i+1,j+1 -> i+1, j
                
                triangle1 = np.array([
                    ij_to_index(i,j),
                    ij_to_index(i+1,j),
                    ij_to_index(i,j+1)
                    
                ])
                triangle2 = np.array([
                    ij_to_index(i,j+1),
                    ij_to_index(i+1,j+1),
                    ij_to_index(i+1,j)
                ])
                
                triangles.append(triangle1)
                triangles.append(triangle2)
                
                
        mesh = o3d.geometry.TriangleMesh()
        mesh.vertices = o3d.utility.Vector3dVector(allLayers)
        mesh.triangles = o3d.utility.Vector3iVector(triangles)
        mesh.compute_vertex_normals()

        if showPlotly:
            o3d.visualization.draw_plotly([mesh])

        self.mesh = mesh
        
        return mesh
    
    def export_mesh_stl(self, filename = 'out'):
        """generate an stl from o3d mesh

        Args:
            filename (str, optional): filename to save as. Defaults to 'out'.
        """
        o3d.io.write_triangle_mesh(f"{filename}.stl", self.mesh)
        
if __name__ == "__main__":
    
    
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
    sim = growth_sim(cl, threeD=True, growingKwargs=growingParams, saveToFile=True)
    
    # load clump from file
    sim.load_clump_from_layers("/home/skilmartin/code/differential-growth/saves/example_save/layers.pkl",layer=48)
    sim.run_vase_smoothing(10, sim.gaussian_smooth, smoothingKwargs={'points': 11}, displayIts=20)
    sim.trim_vase(180)
    sim.convert_layers_to_mesh(showPlotly=True)
    
  
    
   
   
        
    
    
    
    
    
    
   
    