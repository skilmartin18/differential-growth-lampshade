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
_perturbation_func = particle.perturbation(prob=0.75, magnitude=0.5)
_general_force = particle.force(rDistance=75, aDistance=100, k=30, clamp = 1e999)
SIZE = 8
RADIUS = 40
### RUN SIMULATION ###

## create clump
cl = particle.clump(size=SIZE, radius=RADIUS, _perturbation=_perturbation_func, _force=_general_force)



## grow clump
for i in range(ITERATIONS):
    
    cl.run_sim_iteration(incrementZ=False)
    
    
    ## draw clump
    
    # if i % 100 == 0:
    #     ax.cla()
    ax.cla()
    positions =  [particle.positionVec for particle in cl.particles]
    positions.append(positions[0])
    
    xPos = [position[0] for position in positions]
    yPos = [position[1] for position in positions]
    zPos = [position[2] for position in positions]
    ax.plot(xPos,yPos, zPos, color = 'red')
    # ax.set_xlim(-400,400)
    # ax.set_ylim(-400,400)
    
    fig.canvas.draw()
    fig.canvas.flush_events()
    
plt.ioff()
plt.show()