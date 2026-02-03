import np as np



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
    
    def __init__(self, size: int = 10, shape: str = 'circle', radius = 50, initial_peturbation: force = None, peturbation: force = None):
        
        self.particles = []
        self.place_particles(size, shape, radius)
        pass
    
    def place_particles(self, size: int, shape: str, radius: float):
        
        # basic version places particles in a circle only
        locations = np.linspace(0,2*3.14159,size)
        
        for i in range(size):
            self.particles.add(
                particle(
                    initialX= radius*np.cos(locations[i]),
                    initialY= radius*np.sin(locations[i])
                )
            )
            
            
    
    def initial_peturbation(self):
        pass
    
    def peturbation_function(self):
        pass


