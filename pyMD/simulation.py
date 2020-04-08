import numpy as np
from .neighborlist import NeighborList


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

        # CUSTOM

        self.pair_styles = []
        self.reporters = []

        self.initialized = False

    def init_forces(self):
        debug("Init Forces")
        for pair_style in self.pair_styles:
            pair_style.connect(self)

    def init_neighborlist(self):
        debug("Init Neighborlist")
        self.neigh = NeighborList(self.positions, self.box,
                                  pair_cutoffs=[f.cutoff for f in self.pair_styles])

    def add_pair_style(self, pair_style):
        debug("add pair_style")
        self.pair_styles.append(pair_style)

    def set_integrator(self, integrator):
        "Sets the integrator"
        debug("set integrator")
        self.integrator = integrator
        self.integrator.connect(self.positions, self.velocities, self.masses)
        self.integrator.setup()

    def add_reporter(self, reporter):
        debug("add reporter")
        reporter.connect(self)
        reporter.setup()
        self.reporters.append(reporter)

    def _init_step(self):
        self.init_neighborlist()
        self.neigh.build_neighborlist()

        self.init_forces()

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

        for step in range(1, n_steps + 1):

            if not step % 1000:
                print("\rStep : {:10d} / {} ({:6.2f}%)".format(step, n_steps, 100.0 * step / n_steps), end='',
                      flush=True)

            # activate neighborlist build
            if not step % nn_every:
                self.nflag = True

            self._step()

            if not step % reporter_every_min:
                for reporter in self.reporters:
                    if not step % reporter.n_dump:
                        reporter.report(step)