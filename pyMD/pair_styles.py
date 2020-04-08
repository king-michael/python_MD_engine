import numpy as np
from abc import abstractmethod

class PairStyle:
    def __init__(self, cutoff):
        self.cutoff = cutoff
        self.cutoff = cutoff
        self.cutsq = cutoff ** 2

    def connect(self, engine):
        # connects integrator with engine
        self.positions = engine.positions
        self.box = engine.box
        self.pe = 0.0

        self.neigh = engine.neigh
        self.neighborlist = engine.neigh.neighborlist
        self.ghost_images = engine.neigh.ghost_images
        self.n_locals = engine.neigh.n_locals

    def fix_pbc(self, delx, ghost_images):
        img = self.neigh.convert_ghostimage_to_3d(ghost_images)
        delx -= np.dot(img, self.box)


    def setup(self):
        # setup variables, reimplement if needed
        pass

    @abstractmethod
    def compute(self, forces):
        # computes forces, needs to be reimplemented
        pass

    def get_current_state(self, return_forces=True):
        "gets the potential energy and forces of the current snapshot"
        forces = np.zeros(self.positions.shape, dtype=np.float64)
        old_pe = self.pe
        self.compute(forces)
        new_pe, self.pe = self.pe, old_pe
        if return_forces:
            return forces, new_pe
        else:
            return new_pe


class LennardJones(PairStyle):
    def __init__(self, epsilon, sigma, cutoff, shift=True):
        super().__init__(cutoff)

        self.epsilon = epsilon
        self.sigma = sigma

        self.c12 = 4 * epsilon * sigma ** 12
        self.c6 = 4 * epsilon * sigma ** 6

        self.f_c12 = 12 * self.c12
        self.f_c6 = 6 * self.c6

        if shift:
            r6inv = 1.0 / cutoff ** 6
            self.offset = r6inv * (self.c12 * r6inv - self.c6)
        else:
            self.offset = 0.0

    def compute(self, forces):
        for i, xtmp in enumerate(self.positions):
            n_local = self.n_locals[i]
            if n_local > 0:
                neighbors = self.neighborlist[i, :n_local]
                delx = xtmp - self.positions[neighbors]

                self.fix_pbc(delx, self.ghost_images[i, :n_local])

                rsq = np.sum(np.power(delx, 2), axis=1)
                idx_in_cutoff = np.where(rsq < self.cutsq)[0]
                if len(idx_in_cutoff) > 0:
                    rsq = np.clip(rsq[idx_in_cutoff],1e-6,None)
                    r2inv = 1.0 / rsq#[idx_in_cutoff]
                    r6inv = r2inv * r2inv * r2inv

                    forcelj = r6inv * (self.f_c12 * r6inv - self.f_c6)
                    forcelj *= r2inv

                    forces[i] += np.sum( np.dot(forcelj,delx[idx_in_cutoff]) )

                    forces[neighbors[idx_in_cutoff]] -= np.dot(forcelj,delx[idx_in_cutoff])

                    self.pe += np.sum(r6inv * (self.c12 * r6inv - self.c6) - self.offset)
