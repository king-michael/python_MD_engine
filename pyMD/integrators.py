import numpy as np
from abc import abstractmethod


class Integrator:
    def __init__(self, dt):
        self.dt = dt

    def connect(self, enigne):
        # connects integrator with engine
        self.positions = enigne.positions
        self.velocities = enigne.velocities
        self.masses = enigne.masses
        self.forces = enigne.forces
        enigne.dt =self.dt

    @abstractmethod
    def setup(self):
        # setup variables
        pass

    @abstractmethod
    def initial_integrate(self):
        # before force evaluations
        pass

    @abstractmethod
    def final_integrate(self):
        # after force evaulations
        pass


class VelocityVerletIntegrator(Integrator):
    constants = {
        'LJ': {'ftm2v': 1.0},
        'real': {'ftm2v': 0.00041839999997254776},
        'metal': {'ftm2v': 9648.533823273016},
    }

    def __init__(self, dt, units='real'):
        super().__init__(dt)
        self.set_unit_system(units)

    def set_unit_system(self, units: str):
        try:
            self.ftm2v = self.constants[units]['ftm2v']
        except KeyError as e:
            raise KeyError('"ftm2v" for Unit("") not found'.format(units))

    def setup(self):
        # setup variables
        self.dtv = self.dt
        self.dtf = 0.5 * self.dt * self.ftm2v
        self.dtfm = np.atleast_2d(self.dtf / self.masses).T

    def initial_integrate(self):
        self.velocities += self.dtfm * self.forces
        self.positions += self.dtv * self.velocities

    def final_integrate(self):
        self.velocities += self.dtfm * self.forces
