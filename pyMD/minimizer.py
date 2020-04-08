import numpy as np

from abc import abstractmethod


class Minimizer:
    constants = {
        'LJ': {'ftm2v': 1.0, 'dt' : 0.005},
        'real': {'ftm2v': 0.00041839999997254776, 'dt': 1.0},
        'metal': {'ftm2v': 9648.533823273016, 'dt': 1.0},
    }

    def __init__(self, n_dump: int=10, units: str='LJ'):
        self.n_dump = n_dump
        self.units = units
        self.set_unit_system(units)

    def set_unit_system(self, units: str):
        self.ftm2v = self.constants[units]['ftm2v']
        self.dt = self.constants[units]['dt']


    def connect(self, engine):
        #print('connect {} with {}'.format(self, engine))
        self.engine = engine
        if hasattr(engine, 'dt'):
            self.dt = engine.dt

    def setup(self):
        pass

    def get_potential_energy(self):
        pe = sum(p.pe for p in self.engine.pair_styles)
        return pe

    @abstractmethod
    def minimize(self, n_steps : int):
        return "Criterial"

    def update_forces(self):
        forces = self.engine.forces
        forces[:] = 0
        for pair_style in self.engine.pair_styles:
            pair_style.compute(forces)


class MonteCarloMinimizer(Minimizer):
    def __init__(self, max_ds: float=0.1, n_dump: int=10, units: str='LJ'):
        super().__init__(n_dump, units)
        self.max_ds = max_ds

    def minimize(self,
                 max_iter: int=100,
                 etol: float=1e-4,
                 ftol: float=1e-6,
                 temp=0.1):
        n_atoms = len(self.engine.positions)
        positions = self.engine.positions.T
        position_old = positions.copy()
        forces = self.engine.forces
        masses = self.engine.masses

        ds = self.max_ds

        force_magnitude = np.zeros(n_atoms, dtype=np.float64)

        if not self.engine.initialized:
            self.engine._init_step()

        self.engine.neigh.apply_pbc()
        self.engine.neigh.build_neighborlist()

        self.update_forces()
        pe = pe_old = self.get_potential_energy()

        max_f_mag = ds / self.dt*2 * self.ftm2v / masses
        print('Run minimzation (max steps : {})'.format(max_iter))
        print("{:6d} | pe : {:g} | ds : {:g}".format(0, pe, ds), flush=True)

        reason = ""
        force_magnitude[:] = np.linalg.norm(forces, axis=1)
        force_max = force_magnitude.max()
        if force_max < ftol:
            forces[:] = 0
            reason = 'ftol'
            return

        step = 0
        for step in range(1, max_iter + 1):
            if not step % self.n_dump:
                print("{:6d} | pe : {:g} | ds : {:g}".format(step, pe, ds), flush=True)

            position_old[:] = positions
            pe_old = pe

            positions += np.clip(forces.T / max_f_mag * ds, -ds, ds)
            self.update_forces()

            pe = self.get_potential_energy()
            pe_diff = pe - pe_old
            if abs(pe_diff) < etol:
                reason = "etol < {}".format(etol)
                break

            force_magnitude[:] = np.linalg.norm(forces, axis=1)
            force_max = force_magnitude.max()
            if force_max < ftol:
                forces[:] = 0
                reason = 'ftol < {}'.format(ftol)
                return

            if pe_diff < 0: # try : np.random.rand()
                continue
            else:
                print('decrease ds')
                ds *= 0.5
                max_f_mag = ds / self.dt * 2 * self.ftm2v / masses
                positions[:] = position_old
                pe = pe_old
                if ds < 1e-5:
                    reason = 'ds < 1e-5'
                    break
        else:
            reason = 'maxit'

        print("finished minization ({} steps) | pe : {:g}".format(step + 1, pe), flush=True)

        forces[:] = 0
        self.engine.neigh.apply_pbc()
        self.engine.neigh.build_neighborlist()
        return reason

