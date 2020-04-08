import numpy as np
from functools import wraps

from .neighborlist import NeighborList
from .minimizer import Minimizer, MonteCarloMinimizer
from .reporters.base import Reporter
from .integrators import Integrator
from .pair_styles import PairStyle

def debug(*args, **kwargs):
    print(*args, **kwargs)


class Simulation():
    nn_settings = {
        'pagesize': 2000,  # max neighborlist size
        'skin': 2.0,  # neighbor skin
        'every': 10,
    }

    def __init__(self, positions, masses, box, velocities=None, timestep=1.0):
        self.n_atoms = len(positions)

        self.positions = positions.astype(np.float64).copy()
        self.masses = masses.astype(np.float64).copy()
        self.box = box

        if velocities is None:
            self.velocities = np.zeros(positions.shape, dtype=np.float64)
        else:
            self.velocities = velocities.astype(np.float64).copy()

        self.forces = np.zeros(self.positions.shape, dtype=np.float64)

        self.timestep = timestep

        self.set_minimizer(MonteCarloMinimizer(n_dump=100, max_ds=0.02))
        # CUSTOM

        self.pair_styles = []
        self.reporters = []

        self.initialized = False

    def _init_forces(self):
        debug("Init Forces")
        for pair_style in self.pair_styles:
            pair_style.connect(self)

    def _init_neighborlist(self):
        debug("Init Neighborlist")
        self.neigh = NeighborList(self.positions, self.box,
                                  pair_cutoffs=[f.cutoff for f in self.pair_styles])

    def add_pair_style(self, pair_style: PairStyle):
        debug("add pair_style")
        self.pair_styles.append(pair_style)
        return self

    def set_integrator(self, integrator: Integrator):
        "Sets the integrator"
        debug("set integrator")
        self.integrator = integrator
        self.integrator.connect(self)
        self.integrator.setup()
        return self

    def minimize(self,
                 max_iter: int=100,
                 etol: float=1e-4,
                 ftol: float=1e-6):
        raise NotImplementedError("No Minimizer set")
        #print(self.minimizer.minimize)
        #self.minimizer.minimize(max_iter, etol, ftol)

    def set_minimizer(self, minimzer: Minimizer):
        debug("set minimizer")
        self.minimizer = minimzer
        self.minimizer.connect(self)
        self.minimizer.setup()

        def wrapped_minimze(f):
            @wraps(f)
            def wrapper(*args, **kwargs):
                reason = f(*args, **kwargs)
                print('Criteria : {} reached'.format(reason))
                return self
            return wrapper
        self.minimize = wrapped_minimze(self.minimizer.minimize)

        return self

    def add_reporter(self, reporter: Reporter):
        debug("add reporter")
        reporter.connect(self)
        reporter.setup()
        self.reporters.append(reporter)
        return self

    def _init_step(self):
        self._init_neighborlist()
        self.neigh.build_neighborlist()

        self._init_forces()

        self.initialized = True

    def _step(self):
        # for step routine see : https://lammps.sandia.gov/doc/Developer.pdf
        self.integrator.initial_integrate()

        if self.nflag:
            self.neigh.apply_pbc()
            self.neigh.build_neighborlist()
            self.nflag = False

        self.forces[:] = 0
        for pair_style in self.pair_styles:
            pair_style.compute(self.forces)

        # CONSTRAINTS
        self.integrator.final_integrate()

    def get_potential_energy(self):
        pe = sum(p.pe for p in self.pair_styles)
        return pe

    def run(self, n_steps):
        if not self.initialized:
            self._init_step()

        for reporter in self.reporters:
            reporter.report(0)

        self.nflag = True
        nn_every = self.neigh.every
        reporter_every_min = min(r.n_dump for r in self.reporters)
        print("Run {} steps".format(n_steps))
        print("Step : {:10d} / {} ({:6.2f}%)".format(0, n_steps, 0.0), end='', flush=True)
        for step in range(1, n_steps + 1):

            if not step % 1000:
                pe = self.get_potential_energy()
                print("\rStep : {:10d} / {} ({:6.2f}%) | pe : {}".format(step, n_steps, 100.0 * step / n_steps, pe), end='',
                      flush=True)

            # activate neighborlist build
            if not step % nn_every:
                self.nflag = True

            self._step()

            if not step % reporter_every_min:
                for reporter in self.reporters:
                    if not step % reporter.n_dump:
                        reporter.report(step)

        return self