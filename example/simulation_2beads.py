import sys
sys.path.insert(0, '..')
import numpy as np
import pyMD


#===================================================================#
# Inputs

box = np.array([
    [26, 0, 0],
    [0, 26, 0],
    [0, 0, 26],
])

positions = np.array([
    [0, 0, 0],
    [1, 0, 0],
], dtype=np.float64)

velocities = np.array([
    [ 5e-5, 0, 0],
    [-5e-5, 0, 0]
], dtype=np.float64)

masses = 18.0*np.ones(len(positions))



#===================================================================#
# Define the Simulation

simulation = pyMD.Simulation(
    positions=positions,
    masses=masses,
    box=box,
    velocities=velocities
)
simulation.set_integrator( pyMD.VelocityVerletIntegrator(1.0) )
simulation.add_pair_style( pyMD.LennardJones(epsilon=2, sigma=3.5, cutoff=11) )
simulation.add_reporter( pyMD.PositionMonitor(n_dump=1000) )
simulation.add_reporter( pyMD.ThermodynamicsReporter(n_dump=1000) )

#===================================================================#
# run the simulation
simulation.run(int(1e6))
