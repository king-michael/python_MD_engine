import numpy as np
from abc import abstractmethod

class Reporter:
    def __init__(self, n_dump):
        self.n_dump = n_dump

    def connect(self, engine):
        self.engine = engine

    @abstractmethod
    def report(self, step):
        pass

class PositionMonitor(Reporter):
    def __init__(self, n_dump=1000):
        super().__init__(n_dump)

        self.trajectory = []

    def report(self, step):
        self.trajectory.append( self.engine.positions.copy() )

class ThermodynamicsReporter(Reporter):
    def __init__(self, n_dump=1000):
        super().__init__(n_dump)

        self.trajectory = []

    def report(self, step):
        pe = self.engine.get_potential_energy()
        self.trajectory.append( (step, pe) )