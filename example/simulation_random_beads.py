import sys
sys.path.insert(0, '..')
import numpy as np
import pyMD


#===================================================================#
# Inputs

n_atoms = 200

system = pyMD.builder.RandomBeadsBox(n_atoms=n_atoms,
                                     sidelength=20,
                                     mode='mc',
                                     mc_r_min=3.0)
box = system.box
positions = system.positions

masses = 1*np.ones(n_atoms)

temperature = 0.02 # just a scaling atm
velocities = pyMD.builder.create_velocities(n_atoms, masses, temperature)


simulation = pyMD.Simulation(
    positions=positions,
    masses=masses,
    box=box,
    velocities=velocities
)

simulation.set_integrator( pyMD.VelocityVerletIntegrator(0.005, units='LJ') )
simulation.add_pair_style( pyMD.LennardJones(epsilon=0.25, sigma=1, cutoff=2.5) )
#simulation.add_reporter( pyMD.InMemoryReporter(n_dump=1000) )
#simulation.add_reporter( pyMD.ThermodynamicsReporter(n_dump=1000) )

simulation.add_reporter( pyMD.LammpsTrajectoryReporter('test.lammpstrj', n_dump=100))

simulation.minimize(10000)
#===================================================================#
# run the simulation
simulation.run(int(1e5))
