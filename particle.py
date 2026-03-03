import numpy as np
import matplotlib.pyplot as plt
import time
from typing import Self
import logging
import multiprocessing
import random

logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s - %(levelname)s - %(message)s',
)

log = logging.getLogger("differential-growth")


"""
A collection of objects used to simulate 2d differential growth

This module supplies force, perturbation, particle, and clump objects. These
allow the simple simulation of differential growth. 

"""


class force:
    """The force function is used in clump simulation to apply forces between particles.
    The general force descried here is effectivley a 'spring' force: F = -kx. In this case the
    force of attraction and repulsion are calculated the same, in some implementations they can have 
    different implementations
    """
    def __init__(self, rDistance = 75, aDistance = 100, k = 20, neighbors = 1, clamp = 1e99):
        """generic force for application to particles. 

        Args:
            rDistance (int, optional): distance within which particles will begin to repel eachother. Defaults to 75.
            aDistance (int, optional): distance beyond which which particles will attract . Defaults to 100.
            k (int, optional): sptring constant. Defaults to 20.
            clamp (_type_, optional): max force that can be applied. Defaults to 1e99.
        """
        
        self.repulsionDistance = rDistance
        self.attractionDistance = aDistance
        self.springConstant = k
        self.clamp = clamp
        self.neighbors = neighbors
        
    
    def calc_force_repulse(self, norm: float) -> float:
        """calculates repulsive force

        Args:
            norm (float): norm of distance vector between two points (how far away are two points)

        Returns:
            float: force
        """
        if norm == 0:
            return 0
        force = (self.springConstant/(norm*norm)) 
        
        if force > self.clamp:
            force = self.clamp
            
        return -force
    
    def calc_force_attract(self, norm: float) -> float:
        """calculates attractive force

        Args:
            norm (float): norm of distance vector between two points (how far away are two points)

        Returns:
            float: force
        """
        if norm == 0:
            return 0
        force = (self.springConstant/(norm*norm)) 
        
        if force > self.clamp:
            force = self.clamp
            
        return force
    
class perturbation():
    """perturbation is used to randomly place starting points 
    This heavily effects clump growth patterns
    """
    
    def __init__(self, prob = 0.75, magnitude=1):
        """initialise a perturbation function

        Args:
            prob (float, optional): probability that a point will be perturbed. Defaults to 0.75.
            magnitude (int, optional): multiplier for perturbation amount which is 0-1. Defaults to 1.
        """
        self.magnitude = magnitude
        self.prob = prob
        
    def perturb(self, vector: np.array) -> np.array:
        """moves a vector by a small amount

        Args:
            vector (np.array): position vector 3d

        Returns:
            np.array: position vector after perturbation
        """
        if random.random() < self.prob:
            amount = random.random()
            return vector + np.array([1,1,0])*amount*self.magnitude
        else:
            return vector
            
class particle:
    """A single particle in a clump
    """
    
    def __init__(self, initialPos: list[float] = [0,0,0], **kwargs):
        """initialise a particle at a position

        Args:
            initialPos (list[float], optional): initial position. Defaults to [0,0,0].
        Kwargs:
            mass: mass of the particle, used in force-> velocity calculations
            friction: particle friction- velocity loss per iteration as a percentage
            life: total amount of iterations that a particle will remain active for
            growthRate: multiplier for particle velocity->movement calculations
        """
        self.positionVec = np.array(initialPos)
        self.velocityVec = np.array([0,0,0])
        self.mass = kwargs.get('mass', 1)
        self.friction = kwargs.get('friction', 0.5)
        self.life = kwargs.get('life', 1e99)
        self.growthRate = kwargs.get('growthRate', 1)
        self.kwargs = kwargs
        
        
    def update_velocity(self, particlesVectors: np.array, _force: callable):
        """update the particles velocity 

        Args:
            particlesVectors (np.array): an array of position vectors, vectors are taken between where THIS particle is and
            where ALL OTHER particles exterting force are.
            _force (force.calc_force_attract | force.calc_force_repulse): force to be applied to particle based on distance from neighboing particles
        """
        
        start = time.perf_counter()
        log.debug("calculating forces on particle...")
        log.debug(f"particles vectors {particlesVectors}")

        # calculate unit vectors
        particlesMag = np.sqrt(np.sum(particlesVectors**2,1))
        particlesUnitVectors = particlesVectors/particlesMag[:,np.newaxis]
        log.debug(f"repulsed unit vectors are {particlesUnitVectors}")
        
        # calculate force vectors
        particlesForces = particlesUnitVectors * np.array(list(map(_force,particlesMag)))[:,np.newaxis]
        log.debug(f"forces are {particlesForces}")
        
        # calculate resultant force
        particlesTotalForce = sum(particlesForces)
        log.debug(f"total force is {particlesTotalForce}")
        
        # convert force to velocity vector
        self.velocityVec = self.velocityVec + (particlesTotalForce/self.mass)
        log.debug(f"velocity is {self.velocityVec}")
        log.debug(f"done calculating in {time.perf_counter()-start:.4f}")
        
    def update_position(self):
        """updates the particle position based on its velocity for one iteration
        """
        
        if self.life > 0:
            self.positionVec = self.positionVec + self.velocityVec*self.growthRate
            self.velocityVec = self.velocityVec*self.friction
            self.life = self.life - 1
        # to stop ressurected partices from surging forwards, you cannot grow while dead
        if self.life <= 0:
            self.velocityVec = np.array([0,0,0])
 
         
class clump:
    """a clump of particles is what actually creates a differential growth sim.
    A clump is an ORDERED list of particles
    """
    
    def __init__(self, size: int = 10, radius = 40, insertDistance = 15, _force: force = force(), _perturbation: perturbation = perturbation(), particleKwargs = {}):
        """initialise the clump

        Args:
            size (int, optional): how many particles are in starting clump. Defaults to 10.
            radius (int, optional): radius of initial clump (circle). Defaults to 40.
            insertDistance (int, optional): distance between particles after which a new one is placed between. Defaults to 15.
            _force (force, optional): force for calculating attraction and repulsion. Defaults to force().
            _perturbation (perturbation, optional): perturbation function for initial conditions. Defaults to perturbation().
            particleKwargs (dict, optional): kwargs dictionary to be passed to particles. Defaults to {}.
        """
        self.particles = []
        self.insertDistance = insertDistance
        self._force = _force
        self._perturbation = _perturbation
        self.iterations = 0
        self.totalTime = 0
        self.particleKwargs = particleKwargs
        self.stop = False
        self.place_particles(size,radius)
        
            
    def place_particles(self, size: int, radius: float):
        """place initial clump

        Args:
            size (int): how many particles in clump
            radius (float): radius of initial clump circle
        """
        
        locations = np.linspace(0,2*3.14159,size+1)
    
        for i in range(size):
            
            self.particles.append(
                particle(
                    [
                    radius*np.cos(locations[i]),
                    radius*np.sin(locations[i]),
                    0
                    ],
                    **self.particleKwargs
                )
            )
            log.debug(f"placing particle @{self.particles[i].positionVec}")
        # perturb
        for part in self.particles:
            part.positionVec = self._perturbation.perturb(part.positionVec)
      
    def create_growing_zone(self, percentCoverage: float, life: int, growthRate: float):
        """The aim of this function is to create geometrically contiguous zones of reinvigorated growth by 
        increasing their life. The effects of this function are only visible if using normal particles with low
        lifespans.

        Args:
            percentCoverage (float): what percentage of the total perimeter of the clump is 
            available to be selected 
            life (int): how much extra life should be given to affected particles
            growthRate (float): velcoity multiplier (unused)
        """
      
        # randomly choose a starting cell
        start = int(np.floor(random.random()*len(self.particles)))
        coverage = int(np.floor(percentCoverage*len(self.particles)))
        # edit all cells from starting cell to final cell based on percentCoverage
        for i in range(start, start+coverage):
            self.particles[i%len(self.particles)].life = self.particles[i%len(self.particles)].life + life

        log.info(f"inserted growing zone from {start} to {start+coverage}")
        
    def update_particle(self, particle_index: int) -> particle:
        """update a single particle based on all others in clump

        Args:
            particle_index (int): the index of the particle to be updated

        Returns:
            particle: updated particle
        """
          
        log.debug(f"updating velocity on particle {particle}")
        
        i = particle_index
        # create an array of particle position vectors
        particlesPositions = []
        particlesPositions.append([particle.positionVec for particle in self.particles])
        particlesPositions = np.array(particlesPositions[0])
       
        # if we calculate the norm, we can check boundaries and apply the correct forces
        particlesVectors = particlesPositions - self.particles[i].positionVec
        particlesNorm = np.sum(particlesVectors**2,1)
        
        # these are the particles to be repulsed
        repulsedParticles = (particlesNorm < (self._force.repulsionDistance**2)) 
        selfParticle = particlesNorm != 0
        
        allParticles = repulsedParticles & selfParticle
        log.debug(f"all repulsed particles {repulsedParticles}")
        log.debug(f"particles that arent self {selfParticle}")
        log.debug(f"combined set: {allParticles}")
        # this will calculate repuslive force
        self.particles[i].update_velocity(particlesVectors[allParticles],self._force.calc_force_repulse)
        
        # only attract particles on left and right of current
        # you would think that since all particles repel, all particles should be factored into attraction..
        # apparently no? it makes the sim ugly so 
        attractedParticles = np.array([self.particles[i-1].positionVec,self.particles[(i+1)%len(self.particles)].positionVec])
        n = self._force.neighbors
        attractedParticlesPositions = attractedParticles[i-n:(i+n)%len(self.particles)]
        attractedParticlesVectors = attractedParticlesPositions - self.particles[i].positionVec
        
        attractedParticlesVectors = attractedParticles - self.particles[i].positionVec
        particlesNorm = np.sum(attractedParticlesVectors**2,1)
        log.debug(f"attracted particles {attractedParticles}")
        self.particles[i].update_velocity(attractedParticlesVectors[particlesNorm > self._force.repulsionDistance],self._force.calc_force_attract)
    
        return self.particles[i]
    
    def update(self):
        """update all particles in clump then insert
        """
        
        start  = time.perf_counter()
        
        # multiprocessed calculate all particles velocity
        # since each calculation is independent
        particlesIndex = list(range(0,len(self.particles),1))
        with multiprocessing.Pool() as pool:
            results = pool.map(self.update_particle,particlesIndex)
       
        self.particles = results   
        log.info("updated particles...")
        
        
        for particle in self.particles:
            particle.update_position()
        log.info("updated particle positions...")
        
        
        log.info("inserting new nodes...")
        self.insert()
        log.info(f'node count after insertion {len(self.particles)}')
        log.debug(f'updated in {time.perf_counter()-start:.4f}')    
    
    def insert(self):
        """when two particles are too far apart, they need new particles inserted between them
        this does that...
        """
        
        i = 0
        # while instead of for allows to compensate for inserting new particles and checking against newly added particles
        while i < len(self.particles):
            
            # check distance between particle and its next neighbor
            diff = self.particles[(i+1)%len(self.particles)].positionVec - self.particles[i].positionVec
            mag = np.linalg.norm(diff)
            
            # if norm of difference between vectors higher than threshold add particle (insert)
            if mag > self.insertDistance:
                log.debug("inserted")
                newParticlePos = self.particles[i].positionVec + 0.5*diff
                newParticle = particle(newParticlePos, **self.particleKwargs)
                self.particles.insert(i+1,newParticle)
            else: 
                # we only move to the next particle if a new one WASN'T added
                # this covers the case that the distance between the old particle and new 
                # particle is too big, even after one insert
                log.debug("no insertion required")
               
                i = i+1 
             
    def run_sim(self, iterations: int):
        """grows the clump for a certain amount of iterations

            Args:
                iterations (int): number of iterations to run for      
        """
        
        start = time.perf_counter()
        for i in range(iterations):
            log.info(f'running iteration {i+1}')
            self.update()
            
            if self.stop:
                break
            
        log.info(f"total runtime {time.perf_counter()-start}")
        log.info(f"avg time per iteration {(time.perf_counter()-start)/iterations}")
    
    def run_sim_iteration(self):
        """run a single iteration of the clump
        """
        start = time.perf_counter()
        self.iterations = self.iterations + 1
        log.info(f'running iteration {self.iterations}')
        self.update() 
        currentItTime = (time.perf_counter()-start)
        self.totalTime = self.totalTime + currentItTime
        log.info(f"total runtime {self.totalTime}")
        log.info(f"avg time per iteration {self.totalTime/self.iterations}")
        log.info(f"curent it time {currentItTime}\n\n\n")
        

    


