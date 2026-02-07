
import matplotlib.pyplot as plt
import particle
import logging
log = logging.getLogger("logger")

class growth_sim():
    
    def __init__(self, clump, threeD=False):
        
        plt.ion() 
        self.fig = plt.figure(figsize=(10, 8))
        self.threeD = threeD
        if threeD:
            self.ax = self.fig.add_subplot(projection='3d')
        else:
            self.ax = self.fig.add_subplot()
            
        self.colourMap = plt.cm.viridis 
        self.clump = clump
        pass
    
        
    def grow_clump(self, iterations, animate = True, showProgress = False, saveToFile = False):
        
        for i in range(iterations):
    
            self.clump.run_sim_iteration(incrementZ=False)
            self.draw_clump(showProgress=showProgress, colour=self.colourMap(i/iterations))
        
        plt.ioff()
        plt.show()
    
    def grow_vase(self):
        
        pass
    
    def draw_clump(self, showProgress = False, showParticles = False, colour = 'red'):
        
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
        self.ax.set_xlim(-400,400)
        self.ax.set_ylim(-400,400)
        
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()
        
        pass
    
    

if __name__ == "__main__":
    
    
    ### CLUMP PARAMS ###
    ITERATIONS = 300
    _perturbation_func = particle.perturbation(prob=1, magnitude=2)
    _general_force = particle.force(rDistance=75, aDistance=100, k=30, neighbors=1, clamp = 1e999)

    SIZE = 8
    RADIUS = 40
    INSERT_DIST = 15

    MASS = 1
    FRICTION = 0.5
    LIFE = 60

    particleParams = {
        'mass': MASS,
        'friction': FRICTION,
        'life': LIFE
    }

    ### RUN SIMULATION ###

    ## create clump
    log.setLevel(logging.WARNING)

    
    cl = particle.clump(size=SIZE, radius=RADIUS, insertDistance=INSERT_DIST, _perturbation=_perturbation_func, _force=_general_force, particleKwargs= particleParams)
    
    sim = growth_sim(cl, threeD=False)
    
    sim.grow_clump(ITERATIONS, showProgress=True)
    
    