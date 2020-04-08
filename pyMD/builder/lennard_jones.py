import numpy as np
from .boxes import RectengularBox
from .velocities import create_velocities


class RandomBeadsBox(RectengularBox):

    def __init__(self, n_atoms, sidelength, mode='mc', mc_r_min=1.0):
        """

        Parameters
        ----------
        n_atoms
        sidelength
        mode : str, optional
            mode can be ('mc', 'plain')
            'mc' uses Monte-Carlo moves to place the atoms.
            'plain' just place them randomly
        """
        super().__init__(n_atoms,
                         a=sidelength,
                         b=sidelength,
                         c=sidelength)


        if mode == 'mc':
            self.place_atoms_via_montecarlo(self.positions, self.box, r_min=mc_r_min)
        else:
            self.place_atoms_random(self.positions, self.box)

    def place_atoms_random(self, positions, box):
        print("place atoms randomly in the box")
        # Handle box
        box_vectors = np.linalg.norm(box, axis=1)
        positions[:] = np.random.rand(*positions.shape)*box_vectors
        return positions

    def place_atoms_via_montecarlo(self, positions, box, r_min=1.0):
        print('place atoms via MonteCarlo moves')
        n_atoms = len(positions)

        # Handle box
        (lx, _, _), (xy, ly, _), (xz, yz, lz) = box
        xlo, xhi = 0.0, lx
        ylo, yhi = 0.0, ly
        zlo, zhi = 0.0, lz
        box_vectors = np.linalg.norm(box, axis=1)

        cutoff_sq = r_min**2
        for i in range(n_atoms):
            while True:
                test_position = np.random.random(3)*box_vectors-np.array([xlo, ylo, zlo])

                delx = test_position - positions[:i]
                self.fix_pbc(delx, box)

                rsq = np.sum(np.power(delx, 2), axis=1)

                if not np.any(rsq < cutoff_sq):
                    positions[i] = test_position
                    break

        return positions