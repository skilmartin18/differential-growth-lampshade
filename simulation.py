import matplotlib.pyplot as plt
import logging

import particle


log = logging.getLogger("logger")
log.setLevel(logging.INFO)
    
plt.ion() # Turn on interactive mode
# fig, ax = plt.subplots()


fig = plt.figure(figsize=(10, 8))
# Add a 3D subplot
ax = fig.add_subplot(projection='3d')



### SIMULATION PARAMS ###
ITERATIONS = 1200
_perturbation_func = particle.perturbation(prob=0.75, magnitude=2)
_general_force = particle.force(rDistance=75, aDistance=100, k=30, clamp = 1e999)

SIZE = 9
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
cl = particle.clump(size=SIZE, radius=RADIUS, insertDistance=INSERT_DIST, _perturbation=_perturbation_func, _force=_general_force, particleKwargs= particleParams)


## grow clump
for i in range(ITERATIONS):
    
    cl.run_sim_iteration(incrementZ=False)
    
    
    ## draw clump
    
    # if i % 100 == 0:
    #     ax.cla()
    for p in self.clump.particles:
        ax.plot(particle.positionVec[0],particle.positionVec[1], marker='o', color = 'red')
        
    ax.cla()
    positions =  [particle.positionVec for particle in cl.particles]
    positions.append(positions[0])
    
    xPos = [position[0] for position in positions]
    yPos = [position[1] for position in positions]
    zPos = [position[2] for position in positions]
    ax.plot(xPos,yPos, zPos, color = 'green')
    ax.set_xlim(-400,400)
    ax.set_ylim(-400,400)
    
    fig.canvas.draw()
    fig.canvas.flush_events()
    
plt.ioff()
plt.show()