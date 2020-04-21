import sys
sys.path.insert(0, '..')
import numpy as np
import pyMD

distances = np.linspace(1.0, 100, num=1000)
pe = np.zeros(distances.size)

box = [
    [26, 0, 0],
    [0, 26, 0],
    [0, 0, 26],
]

positions = np.array([
    [0, 0, 0],
    [1, 0, 0],
], dtype=np.float64)

velocities = np.array([
    [5e-5, 0, 0],
    [-5e-5, 0, 0]
], dtype=np.float64)

masses = 18.0 * np.ones(len(positions))




def theory_lj(distances, epsilon, sigma, cutoff, shift=True):
    pe = np.zeros(distances.size)
    mask = distances<cutoff
    pe[mask] = 4*epsilon * (np.power(sigma / distances[mask], 12) - np.power(sigma / distances[mask], 6))
    pe[mask] -=  4*epsilon * (np.power(sigma / cutoff, 12) - np.power(sigma / cutoff, 6))
    return pe


def test_LennardJones():

    pe_theory = theory_lj(distances, epsilon=2, sigma=3.5, cutoff=11)

    simulation = pyMD.Simulation(
        positions=positions,
        masses=masses,
        box=box,
        velocities=velocities
    )
    simulation.set_integrator(pyMD.VelocityVerletIntegrator(1.0))
    simulation.add_pair_style(pyMD.LennardJones(epsilon=2, sigma=3.5, cutoff=11))
    simulation._init_step()

    simulation.positions[0] = (0, 0, 0)

    for i, d in enumerate(distances):
        simulation.positions[1] = (d, 0, 0)
        # simulation.neigh.apply_pbc()
        # simulation.neigh.build_neighborlist()

        pe[i] = simulation.pair_styles[0].get_current_state(return_forces=False)


    np.testing.assert_almost_equal(pe, pe_theory)

if __name__ == '__main__':
    pe_theory = theory_lj(distances, epsilon=2, sigma=3.5, cutoff=11)

    simulation = pyMD.Simulation(
        positions=positions,
        masses=masses,
        box=box,
        velocities=velocities
    )
    simulation.set_integrator(pyMD.VelocityVerletIntegrator(1.0))
    simulation.add_pair_style(pyMD.LennardJones(epsilon=2, sigma=3.5, cutoff=11))
    simulation._init_step()

    simulation.positions[0] = (0, 0, 0)

    for i, d in enumerate(distances):
        simulation.positions[1] = (d, 0, 0)
        # simulation.neigh.apply_pbc()
        # simulation.neigh.build_neighborlist()

        pe[i] = simulation.pair_styles[0].get_current_state(return_forces=False)

    import matplotlib.pyplot as plt

    plt.plot(distances, pe, label='pyMD')

    plt.plot(distances, pe_theory, ls='--', label='theory')
    plt.ylim(pe.min()*1.1, -pe.min()*1.1)
    plt.xlabel('distance')
    plt.ylabel('potential energy')
    plt.legend()
    plt.show()