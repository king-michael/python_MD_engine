import sys
import os
sys.path.insert(0, '..')
import numpy as np
import pyMD
from pyMD.misc.lammps_dump_custom import read_dump




fields = ('id', 'type', 'x', 'y', 'z', 'fx', 'fy', 'fz', 'c_peatom')


def get_lammps_data(dump_file, fields=('id', 'type', 'x', 'y', 'z', 'fx', 'fy', 'fz', 'c_peatom')):
    lammps_dump = np.array([data for data in read_dump(dump_file)])

    fields_xyz = [fields.index(fs) for fs in 'xyz']
    fields_fxyz = [fields.index(fs) for fs in ('fx', 'fy', 'fz')]

    force_lammps = lammps_dump[:, :, fields_fxyz]
    pe_lammps = lammps_dump[:, :, fields.index('c_peatom')].sum(axis=1)
    xyz_lammps = lammps_dump[:, :, fields_xyz]

    return xyz_lammps, pe_lammps, force_lammps

def test_simulation_distance_2LJ_beads():
    path_reference = os.path.join('resources', 'distance_2LJ_beads')
    dump_file = os.path.join(path_reference, 'dump.lammpstrj.gz')

    xyz_lammps, pe_lammps, force_lammps = get_lammps_data(dump_file)

    box = [
        [6.0, 0, 0],
        [0, 6.0, 0],
        [0, 0, 6.0],
    ]

    positions = np.array([
        [0, 0, 0],
        [0.4, 0, 0],
    ], dtype=np.float64)

    velocities = np.array([
        [0, 0, 0],
        [0, 0, 0]
    ], dtype=np.float64)

    masses = 1.0 * np.ones(len(positions))

    simulation = pyMD.Simulation(
        positions=positions,
        masses=masses,
        box=box,
        velocities=velocities
    )
    simulation.nn_settings['skin'] = 0.3
    simulation.set_integrator( pyMD.VelocityVerletIntegrator(0.005))
    simulation.add_pair_style( pyMD.LennardJones(epsilon=1, sigma=1.0, cutoff=2.5, shift=False))
    simulation._init_step()

    pe_pyMD = np.zeros(pe_lammps.shape)
    force_pyMD = np.zeros(force_lammps.shape)

    for i, xyz in enumerate(xyz_lammps):
        simulation.positions[:] = xyz

        simulation.neigh.build_neighborlist()

        force, pe = simulation.pair_styles[0].get_current_state(return_forces=True)
        pe_pyMD[i] = pe
        force_pyMD[i] = force

    np.testing.assert_array_almost_equal(pe_lammps, pe_pyMD, err_msg='Potential Energy mismatch')
    np.testing.assert_array_almost_equal(force_lammps, force_pyMD, err_msg='Forces mismatch')


if __name__ == '__main__':
    test_simulation_distance_2LJ_beads()
