import numpy as np

from abc import abstractmethod


class Minimizer:
    constants = {
        'LJ': {'ftm2v': 1.0},
        'real': {'ftm2v': 0.00041839999997254776},
        'metal': {'ftm2v': 9648.533823273016},
    }

    def __init__(self, n_dump: int=10, units: str='LJ'):
        self.n_dump = n_dump
        
        self.set_unit_system(units)

    def set_unit_system(self, units: str):
        try:
            self.ftm2v = self.constants[units]['ftm2v']
        except KeyError as e:
            raise KeyError('"ftm2v" for Unit("") not found'.format(units))

    def connect(self, engine):
        #print('connect {} with {}'.format(self, engine))
        self.engine = engine

    def setup(self):
        pass

    def get_potential_energy(self):
        pe = sum(p.pe for p in self.engine.pair_styles)
        return pe

    @abstractmethod
    def minimize(self, n_steps : int):
        pass


class SimpleMinimizer(Minimizer):
    def __init__(self, max_ds: float=0.1, n_dump: int=10, units: str='LJ'):
        super().__init__(n_dump, units)
        self.max_ds = max_ds

    def minimize(self,
                 max_iter: int=100,
                 etol: float=1e-4,
                 ftol: float=1e-6):
        n_atoms = len(self.engine.positions)
        positions = self.engine.positions.T
        forces = self.engine.forces

        force_magnitude = np.zeros(n_atoms, dtype=np.float64)

        if not self.engine.initialized:
            self.engine._init_step()

        self.engine.neigh.apply_pbc()
        self.engine.neigh.build_neighborlist()

        self.update_forces()
        pe = old_pe = self.get_potential_energy()

        print('Run minimzation (max steps : {})'.format(max_iter))
        print("{:6d} | pe : {:g}".format(0, pe), flush=True)

        force_magnitude[:] = np.sum(np.power(forces, 2), axis=1)
        force_max = force_magnitude.max()
        if force_max < ftol:
            return

        step=0
        for step in range(1, max_iter+1):


            positions += self.max_ds / force_magnitude.max() * force_magnitude
            self.update_forces()

            pe = self.get_potential_energy()
            if not step % self.n_dump:
                print("{:6d} | pe : {:g}".format(step, pe), flush=True)


            if abs(pe-old_pe) < etol:

                break
            old_pe = pe

            force_magnitude[:] = np.sum(np.power(forces, 2), axis=1)
            force_max = force_magnitude.max()
            if force_max < ftol:
                break
        print("finished minization ({} steps) | pe : {:g}".format(step + 1, pe), flush=True)


        self.engine.neigh.apply_pbc()
        self.engine.neigh.build_neighborlist()


    def update_forces(self):
        forces = self.engine.forces
        forces[:] = 0
        for pair_style in self.engine.pair_styles:
            pair_style.compute(forces)