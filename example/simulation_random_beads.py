import sys
sys.path.insert(0, '..')
import numpy as np
import pyMD


#===================================================================#
# Inputs

box = np.array([
    [10, 0, 0],
    [0, 10, 0],
    [0, 0, 10],
])

n_atoms = 200

positions = np.random.rand(n_atoms,3)
positions[:] = np.dot(positions, box.T)

masses = 1*np.ones(n_atoms)

temperature = 0.48 # just a scaling atm
velocities = temperature*np.random.rand(n_atoms,3)

simulation = pyMD.Simulation(
    positions=positions,
    masses=masses,
    box=box,
    velocities=velocities
)
simulation.set_integrator( pyMD.VelocityVerletIntegrator(1.0) )
simulation.add_pair_style( pyMD.LennardJones(epsilon=2, sigma=1, cutoff=2.5) )
#simulation.add_reporter( pyMD.InMemoryReporter(n_dump=1000) )
#simulation.add_reporter( pyMD.ThermodynamicsReporter(n_dump=1000) )

simulation.add_reporter( pyMD.LammpsTrajectoryReporter('test.lammpstrj', n_dump=1000))

#===================================================================#
# run the simulation
simulation.run(int(1e5))