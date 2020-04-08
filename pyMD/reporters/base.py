from abc import abstractmethod

class Reporter:
    def __init__(self, n_dump):
        """
        Order
        -----
        * creation : `__init__`
        * after adding to engine : -> `connect` -> `setup`
        * evoke every `n_dump` -> `report`

        Parameters
        ----------
        n_dump
        """
        self.n_dump = n_dump

    def connect(self, engine):
        self.engine = engine

    def setup(self):
        "called after connect and befor use"
        pass

    @abstractmethod
    def report(self, step):
        pass