import numpy as np
import matplotlib.pyplot as plt



class particle:
    
    def __init__(self, initialX: float = 0, initialY: float = 0, initialZ: float = 0):
        
        self.positionVec = [initialX, initialY, initialZ]
        self.velocityVec = [0,0,0]
        self.accelVec = [0,0,0]
        
    def update_velocity():
        pass
    
    def update_position():
        pass

class force:
    
    def calc_force(particle: particle, particles: list[particle]) -> list[float]:
        pass

class clump:
    
    def __init__(self, size: int = 10, shape: str = 'circle', radius = 50, insertDistance = 10, initial_peturbation: force = None, peturbation: force = None):
        
        self.particles = []
        self.insertDistance = insertDistance
        self.place_particles(size, shape, radius)
        pass
    
    def place_particles(self, size: int, shape: str, radius: float):
        
        # basic version places particles in a circle only
        locations = np.linspace(0,2*3.14159,size)
        
        for i in range(size):
            self.particles.append(
                particle(
                    initialX= radius*np.cos(locations[i]),
                    initialY= radius*np.sin(locations[i])
                )
            )
    
    def update(self):
        
        # calculate all force vectors
        
        for i in range(len(self.particles)):
            pass
            # calculate respulsive force
            # repulsive force will be calculated against ALL particles
            
            
            # calculate attractive force
            # attractive force will only be calculated on two nearest list neighbors
            
            
            # calculate peturbation
            
        
        # update positions based on force
        
        
        # insert
        
        
        
        
        
        pass
    
    def insert(self):
        
        i = 0
        
        while i < len(self.particles):
            
            # if evaluating i = 1, and requires insert, inserted at i = 2. Next loop iterates with i = 2 against i = 3 which was
            # previous i =2 particle. No particles are skipped.
            
        
            diff = np.array(self.particles[i+1]) - np.array(self.particles[i])
            mag = np.linalg.norm(diff)
            
            if mag > self.insertDistance:
                
                newParticle = self.particles[i]
        
            
            # if l2 norm of difference between vectors higher than threshold add particle
            pass
            
            
    def initial_peturbation(self):
        pass
    
    def peturbation_function(self):
        pass
    
    def draw(self):
        
        for particle in self.particles:
            plt.plot(particle.positionVec[0],particle.positionVec[1], marker='o')
        
        plt.show()
        

if __name__ == "__main__":
    
    cl = clump()
    cl.draw()


