from .base import Reporter

class InMemoryTrajectoryReporter(Reporter):
    def __init__(self, n_dump=1000):
        super().__init__(n_dump)

        self.trajectory = []

    def report(self, step):
        self.trajectory.append( self.engine.positions.copy() )


class InMemoryThermodynamicsReporter(Reporter):
    def __init__(self, n_dump=1000):
        super().__init__(n_dump)

        self.trajectory = []

    def report(self, step):
        pe = self.engine.get_potential_energy()
        self.trajectory.append( (step, pe) )