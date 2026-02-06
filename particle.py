import numpy as np
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time
from typing import Self
import logging
import multiprocessing
import random
logging.basicConfig(
    # filename='app.log', 
    level=logging.INFO, 
    format='%(asctime)s - %(levelname)s - %(message)s',
    # filemode='a' # 'a' for append (default), 'w' to overwrite each time
)

log = logging.getLogger("logger")

class force:
    
    def __init__(self, rDistance = 75, aDistance = 100, k = 20, clamp = 100):
        self.repulsionDistance = rDistance
        self.attractionDistance = aDistance
        self.springConstant = k
        self.clamp = clamp
        

    def calc_force(self, norm):
        
        if norm == 0:
            return 0
        force = (self.springConstant/(norm*norm)) 
        
        if force > self.clamp:
            force = self.clamp
            
        return -force
    
    def calc_force_repulse(self, norm):
        
        if norm == 0:
            return 0
        force = (self.springConstant/(norm*norm)) 
        
        if force > self.clamp:
            force = self.clamp
            
        return -force
    
    def calc_force_attract(self, norm):
        
        if norm == 0:
            return 0
        force = (self.springConstant/(norm*norm)) 
        
        if force > self.clamp:
            force = self.clamp
            
        return force

class perturbation():
    
    def __init__(self, prob = 0.75, magnitude=1):
        self.magnitude = magnitude
        self.prob = prob
        
    def perturb(self, vector: np.array) -> np.array:
        
        # move the whole vector by a bit, with a chance
        
        if random.random() < self.prob:
            amount = random.random()
            return vector + np.array([1,1,0])*amount*self.magnitude
        else:
            return vector
        
             
class particle:
    
    def __init__(self, initialPos: list[float] = [0,0,0], mass = 1, friction = 0.5, life = 1e99):
        
        self.positionVec = np.array(initialPos)
        self.velocityVec = np.array([0,0,0])
        self.mass = mass
        self.friction = friction
        self.life = life
        
    def update_velocity(self, particlesVectors: np.array, _force: force.calc_force):
        
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
        
    def update_position(self, incrementZ = False):
        
        if self.life > 0:
            self.positionVec = self.positionVec + self.velocityVec
            self.velocityVec = self.velocityVec*self.friction
            self.life = self.life - 1
        
        print(incrementZ)
        if incrementZ:
            self.positionVec = self.positionVec + np.array([0,0,1])
            

class clump:
    
    def __init__(self, size: int = 10, radius = 40, insertDistance = 15, _force: force = force(), _perturbation: force = None):
        
        self.particles = []
        self.insertDistance = insertDistance
        self.force_func = force()
        self.p = _perturbation
        self.force_func = _force
        self.iterations = 0
        self.totalTime = 0
        
  
        self.place_particles(size,radius)
        
            
    def place_particles(self, size: int, radius: float):
        # basic version places particles in a circle only
        
        locations = np.linspace(0,2*3.14159,size+1)
        
        for i in range(size):
            
            self.particles.append(
                particle(
                    [
                    radius*np.cos(locations[i]),
                    radius*np.sin(locations[i]),
                    0
                    ]
                )
            )
            log.debug(f"placing particle @{self.particles[i].positionVec}")

        for thing in self.particles:
            
            thing.positionVec = self.p.perturb(thing.positionVec)
            
    def update_particle(self, particle_index):
        
        log.debug(f"updating velocity on particle {particle}")
        
        i = particle_index
        # create an array of particle position vectors
        particlesPositions = []
        particlesPositions.append([particle.positionVec for particle in self.particles])
        particlesPositions = np.array(particlesPositions[0])
        # print(particlesPositions)
        
        # if we calculate the norm, we can check boundaries and apply the correct forces
        particlesVectors = particlesPositions - self.particles[i].positionVec
        particlesNorm = np.sum(particlesVectors**2,1)
        
        # these are the particles to be repulsed
        repulsedParticles = (particlesNorm < (self.force_func.repulsionDistance**2)) 
        selfParticle = particlesNorm != 0
        
        allParticles = repulsedParticles & selfParticle
        log.debug(f"all repulsed particles {repulsedParticles}")
        log.debug(f"particles that arent self {selfParticle}")
        log.debug(f"combined set: {allParticles}")
        # this will calculate repuslive force
        self.particles[i].update_velocity(particlesVectors[allParticles],self.force_func.calc_force_repulse)
        
        # these are the particles to be attracted
        # notRepulsedParticles = particlesNorm > self.force_func.repulsionDistance**2
        # attractedParticles = particlesNorm < self.force_func.attractionDistance**2
        # selfParticle = particlesNorm != 0
        # allParticles = notRepulsedParticles & attractedParticles & selfParticle
        # log.debug(f"all not repulsed particles {notRepulsedParticles}")
        # log.debug(f"particles that are in attraction range {attractedParticles}")
        # log.debug(f"particles that arent self {selfParticle}")
        # log.debug(f"combined set: {allParticles}")
        # # this will calculate attractive force
        # self.particles[i].update_velocity(particlesVectors[allParticles],self.force_func.calc_force_attract)
        
        # only attract particles on left and right of current
        attractedParticles = np.array([self.particles[i-1].positionVec,self.particles[(i+1)%len(self.particles)].positionVec])
        attractedParticlesVectors = attractedParticles - self.particles[i].positionVec
        particlesNorm = np.sum(attractedParticlesVectors**2,1)
        log.debug(f"attracted particles {attractedParticles}")
        self.particles[i].update_velocity(attractedParticlesVectors[particlesNorm > self.force_func.repulsionDistance],self.force_func.calc_force_attract)
        # this will calculate peturbation force
        
    
    
        return self.particles[i]
    
    def update(self, incrementZ = False):
        
        start  = time.perf_counter()
        
        # multiprocessed calculate all particles velocity
        particlesIndex = list(range(0,len(self.particles),1))
        with multiprocessing.Pool() as pool:
            results = pool.map(self.update_particle,particlesIndex)
       
        self.particles = results   
        log.info("updated particles...")
        
        
        for particle in self.particles:
            particle.update_position(incrementZ = incrementZ)
        log.info("updated particle positions...")
        
        log.info("inserting new nodes...")
        self.insert()
        log.info(f'node count after insertion {len(self.particles)}')
        log.debug(f'updated in {time.perf_counter()-start:.4f}')
    
    
    def insert(self):
        
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
                newParticle = particle(newParticlePos)
                self.particles.insert(i+1,newParticle)
            else: 
                # we only move to the next particle if a new one WASN'T added
                # this covers the case that the distance between the old particle and new 
                # particle is too big, even after one insert
                log.debug("no insertion required")
               
                i = i+1 

    def initial_peturbation(self,particle):
        pass
    
    def peturbation_function(self):
        return 0
        
        
    def run_sim(self, iterations):
        
        start = time.perf_counter()
        for i in range(iterations):
            log.info(f'running iteration {i+1}')
            self.update()
        log.info(f"total runtime {time.perf_counter()-start}")
        log.info(f"avg time per iteration {(time.perf_counter()-start)/iterations}")
    
    def run_sim_iteration(self, incrementZ = False):
        start = time.perf_counter()
        self.iterations = self.iterations + 1
        log.info(f'running iteration {self.iterations}')
        self.update(incrementZ = incrementZ) 
        currentItTime = (time.perf_counter()-start)
        self.totalTime = self.totalTime + currentItTime
        log.info(f"total runtime {self.totalTime}")
        log.info(f"avg time per iteration {self.totalTime/self.iterations}")
        log.info(f"curent it time {currentItTime}\n\n\n")
        
        

if __name__ == "__main__":
    log.setLevel(logging.INFO)
    
    plt.ion() # Turn on interactive mode
    fig, ax = plt.subplots()
    
    
    cl = clump()
    
    for i in range(1200):
        cl.run_sim_iteration()
    
        # fig.clf()
        ax.cla()
        
        # draw as particles
        # for particle in cl.particles:
        #     ax.plot(particle.positionVec[0],particle.positionVec[1], marker='o', color = 'red')
        
        # draw as lines
        positions =  [particle.positionVec for particle in cl.particles]
        # connect last point
        positions.append(positions[0])
        xPos = [position[0] for position in positions]
        yPos = [position[1] for position in positions]
        ax.plot(xPos,yPos, color = 'red')
        ax.set_xlim(-500,500)
        ax.set_ylim(-500,500)
        
        fig.canvas.draw()
        # Process events and pause briefly to allow the GUI to update
        fig.canvas.flush_events()
    plt.ioff()
    plt.show()
    


