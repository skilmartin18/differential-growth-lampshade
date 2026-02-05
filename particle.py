import numpy as np
import matplotlib.pyplot as plt
import time
from typing import Self
import logging


logging.basicConfig(
    # filename='app.log', 
    level=logging.INFO, 
    format='%(asctime)s - %(levelname)s - %(message)s',
    # filemode='a' # 'a' for append (default), 'w' to overwrite each time
)

log = logging.getLogger("logger")

class force:
    
    def __init__(self):
        self.repulsionDistance = 25
        self.attractionDistance = 60
        self.springConstant = 1
        
    def calc_force_vector(self, difference: np.array) -> np.array:
        
        norm = np.linalg.norm(difference,ord=2)
        
        if norm == 0:
            return np.array([0,0,0])
        # are the particles close enough to repluse?
        if norm < self.repulsionDistance:
                
            # to calculate force normalise vector to a unit vector
            direction = difference/norm
            
            # equation for force (similar to many forces like gravity/electromagnetic, has fall-off with r^2 and a constant multiplier)
            forceVector = (-self.springConstant/(norm*norm)) * direction
            log.debug(forceVector)
        else:
            forceVector = np.array([0,0,0])
            
        return forceVector

    def calc_force_scalar(self, norm):
        
        if norm == 0:
            return 0
        force = (self.springConstant/(norm*norm)) 
        
        return -force
        # if norm == 0:
        #     return 0
        # # are the particles close enough to repluse?
        # if norm < self.repulsionDistance:
                
        #     # equation for force (similar to many forces like gravity/electromagnetic, has fall-off with r^2 and a constant multiplier)
        #     force = (-self.springConstant/(norm*norm)) 
        #     return force
            
        # elif norm > self.attractionDistance:
        #     return 0
        # else:
        #     return 0
            
             
class particle:
    
    def __init__(self, initialPos: list[float] = [0,0,0], mass = 1, friction = 0.5):
        
        self.positionVec = np.array(initialPos)
        self.velocityVec = np.array([0,0,0])
        self.mass = mass
        self.friction = friction
        
    
    def update_velocity(self, particles: np.array, _force: force = None):
        
        start = time.perf_counter()
        log.debug("calculating forces on particle...")
        log.debug(f"particles compared against {particles}")
        
        
        totalForceVector = np.array([0,0,0])
        
      
        particlesVectors = particles - self.positionVec
        log.debug(f"particle vector {particlesVectors}")
        particlesNorm = np.sum(particlesVectors**2,1)
        log.debug(particlesNorm)
        
        repulsedParticles = particlesNorm < _force.repulsionDistance**2
        repulsedParticles = particlesNorm != 0

        repulsedParticlesVectors = particlesVectors[repulsedParticles]
        
        log.debug(f"repulsed particles vectors {repulsedParticlesVectors}")
   
        repulsedParticlesMag = np.sqrt(np.sum(repulsedParticlesVectors**2,1))
        repulsedParticlesUnitVectors = repulsedParticlesVectors/repulsedParticlesMag[:,np.newaxis]
        log.debug(f"repulsed unit vectors are {repulsedParticlesUnitVectors}")
        repulsedParticlesForces = repulsedParticlesUnitVectors * np.array(list(map(_force.calc_force_scalar,repulsedParticlesMag)))[:,np.newaxis]
        log.debug(f"forces are {repulsedParticlesForces}")
        repulsedParticlesTotalForce = sum(repulsedParticlesForces)
        log.debug(f"total force is {repulsedParticlesTotalForce}")
        self.velocityVec = self.velocityVec + (repulsedParticlesTotalForce/self.mass)
        log.debug(f"velocity is {self.velocityVec}")
        log.debug(f"done calculating in {time.perf_counter()-start:.4f}")
        
    def update_position(self):
        self.positionVec = self.positionVec + self.velocityVec*self.friction
        

class clump:
    
    def __init__(self, size: int = 8, radius = 20, insertDistance = 5, force_func: force = None, initial_peturbation: force = None, peturbation: force = None):
        
        self.particles = []
        self.insertDistance = insertDistance
        self.place_particles(size,radius)
        self.force_func = force()
        
        plt.ion() # Turn on interactive mode
        self.fig, self.ax = plt.subplots()
            
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
    
    def update(self):
        
        start  = time.perf_counter()
        
        # calculate all force vectors
        log.info("updating particles...")
        log.info("updating particle velocities...")
        for i in range(len(self.particles)):
            
            log.debug(f"updating velocity on particle {i}")
            
            # create an array of particle position vectors
            particlesPositions = []
            particlesPositions.append([particle.positionVec for particle in self.particles])
            particlesPositions = np.array(particlesPositions[0])
            # print(particlesPositions)
            
            # if we calculate the norm, we can check boundaries and apply the correct forces
            particlesVectors = particlesPositions - self.particles[i].positionVec
            particlesNorm = np.sum(particlesVectors**2,1)
            
            # create a vector with only repulsed particles and update
            repulsedParticles = particlesNorm < self.force_func.repulsionDistance**2
            repulsedParticles = particlesNorm != 0
            
            # print(repulsedParticles)
            self.particles[i].update_velocity(particlesPositions[repulsedParticles],self.force_func)
            
            # create a vector with only attacted particles
            # this will calculate repuslive force
           
            
            # this will calculate attractive force
            
            # this will calculate peturbation force
            
        
        log.info("updating particle positions...")
        for particle in self.particles:
            particle.update_position()
        
        log.info("inserting new nodes...")
        # insert
        self.insert()
        log.info(f'node count after insertion {len(self.particles)}')
        
        log.info(f'updated in {time.perf_counter()-start:.4f}')
    
    
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

    def initial_peturbation(self):
        pass
    
    def peturbation_function(self):
        pass
    
    def draw(self):
        
        # self.fig.clf()
        self.ax.cla()
        
        # draw as particles
        # for particle in self.particles:
        #     self.ax.plot(particle.positionVec[0],particle.positionVec[1], marker='o', color = 'red')
        
        # draw as lines
        positions =  [particle.positionVec for particle in self.particles]
        # connect last point
        positions.append(positions[0])
        xPos = [position[0] for position in positions]
        yPos = [position[1] for position in positions]
        self.ax.plot(xPos,yPos, color = 'red')
        self.ax.set_xlim(-50,50)
        self.ax.set_ylim(-50,50)
        
        self.fig.canvas.draw()
        # Process events and pause briefly to allow the GUI to update
        self.fig.canvas.flush_events()
        
        
    def run_sim(self, iterations):
        
        start = time.perf_counter()
        for i in range(iterations):
            log.info(f'running iteration {i+1}')
            cl.update()
            cl.draw()
            log.info(f"\n\n\n")
            
        plt.ioff()
        log.info(f"total runtime {time.perf_counter()-start}")
        plt.show()
        
        

if __name__ == "__main__":
    
    log.setLevel(logging.WARNING)
    
    cl = clump()
    cl.run_sim(iterations=30)


