import numpy as np



class Velocities:
    constants = {
        'LJ': {'mvv2e': 1.0, 'boltz' : 1.0},
        'real': {'mvv2e': 2390.0573615334906, 'boltz' : 0.0019872067},
        'metal': {'mvv2e': 0.00010364269, 'boltz' : 8.617343e-05},
    }
    def __init__(self, n_atoms, masses, temperature, units='LJ'):
        self.n_atoms = n_atoms
        self.masses = masses
        self.temperature = temperature

        self.velocities = self.create_velocities(n_atoms, masses, temperature)

    @classmethod
    def calculate_temperature(self, velocities, masses, extra_dof=0, units='LJ'):
        dof = velocities.size - extra_dof
        mvv2e = self.constants[units]['mvv2e']
        boltz = self.constants[units]['boltz']

        if dof > 0:
            tfactor = mvv2e / (dof * boltz)
        else:
            tfactor = 0.0

        temperature = np.sum(np.power(velocities, 2).T * masses)
        temperature*= tfactor
        return temperature

    @classmethod
    def create_velocities(self,
                          n_atoms: int,
                          masses: np.ndarray,
                          temperature: float,
                          mode : str = 'gaussian',
                          extra_dof=0,
                          units='LJ'):
        """

        Parameters
        ----------
        n_atoms : int

        masses : np.ndarray
            Mass array of shape (n_atoms)
        temperature : float
            Temperature
        mode : str, optional
            mode can be 'normal' or 'gaussian'

        Returns
        -------

        """
        assert len(masses) == n_atoms, "Number of masses does not match the number of atoms"
        velocities = np.zeros((3, n_atoms), dtype=np.float64)
        if mode == 'gaussian':
            velocities[:] = np.random.normal(size=velocities.size).reshape((3, -1, ))
        elif mode == 'normal':
            velocities[:] = np.random.random(velocities.size).reshape(3, -1) -0.5
        else:
            raise UserWarning('unknown "mode". Use "normal" or "gaussian"')

        velocities *= 1.0 / np.sqrt(masses)

        t = self.calculate_temperature(velocities.T, masses, extra_dof=extra_dof, units=units)
        Velocities.rescale(velocities.T,t, temperature)
        return velocities.T

    @staticmethod
    def rescale(velocities: np.ndarray, t_old: float, t_new: float):
        factor = np.sqrt(t_new / t_old)
        velocities *= factor
        return velocities


create_velocities = Velocities.create_velocities