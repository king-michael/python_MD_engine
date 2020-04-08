import sys
sys.path.insert(0, '..')
import numpy as np
import pyMD


#===================================================================#
# Inputs

box = np.array([
    [30, 0, 0],
    [0, 30, 0],
    [0, 0, 30],
])

positions = np.array([
    [7, 0, 0],
    [23, 0, 0],
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
simulation.set_integrator( pyMD.VelocityVerletIntegrator(1.0, units='real') )
simulation.add_pair_style( pyMD.LennardJones(epsilon=1e-3, sigma=3.5, cutoff=11) )
simulation.add_reporter( pyMD.InMemoryTrajectoryReporter(n_dump=1000) )
simulation.add_reporter( pyMD.InMemoryThermodynamicsReporter(n_dump=1000) )

# visualzie with :  vmd -lammpstrj test.lammpstrj
simulation.add_reporter( pyMD.LammpsTrajectoryReporter('test.lammpstrj', n_dump=1000))
#===================================================================#

# run the simulation
simulation.run(int(1e6))
