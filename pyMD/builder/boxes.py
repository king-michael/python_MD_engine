from abc import abstractmethod
import numpy as np


class RectengularBox:
    def __init__(self, n_atoms, a, b, c, periodic=(True, True, True)):
        self.box = np.array([
            [a, 0, 0],
            [0, b, 0],
            [0, 0, c],
        ])
        self.box_sides = np.linalg.norm(self.box, axis=1)
        self.box_sides_half = self.box_sides / 2.0

        self.n_atoms = n_atoms

        self.positions = np.zeros((n_atoms, 3), dtype=np.float64)
        self.velocities = np.zeros((n_atoms, 3), dtype=np.float64)

        self.periodic_x, self.periodic_y, self.periodic_z  = periodic


    def apply_pbc(self):
        positions = self.positions
        if self.periodic_x:
            positions[positions[:, 0] < 0.0] += self.box[0]
            positions[positions[:, 0] > self.box_sides[0]] -= self.box[0]
        if self.periodic_y:
            positions[positions[:, 1] < 0.0] += self.box[1]
            positions[positions[:, 1] > self.box_sides[1]] -= self.box[1]
        if self.periodic_z:
            positions[positions[:, 2] < 0.0] += self.box[2]
            positions[positions[:, 2] > self.box_sides[2]] -= self.box[2]

    @staticmethod
    def fix_pbc( delx, box, periodic=(True, True, True)):
        periodic_x, periodic_y, periodic_z = periodic

        box_sides = np.linalg.norm(box, axis=1)
        box_sides_half = box_sides / 2.0

        if periodic_x:
            m_x_lo = delx[:, 0] < -box_sides_half[0]
            delx[m_x_lo] += box[0]

            m_x_hi = delx[:, 0] > box_sides_half[0]
            delx[m_x_hi] -= box[0]

        if periodic_y:
            m_y_lo = delx[:, 1] < -box_sides_half[1]
            delx[m_y_lo] += box[1]

            m_y_hi = delx[:, 1] > box_sides_half[1]
            delx[m_y_hi] -= box[1]

        if periodic_z:
            m_z_lo = delx[:, 2] < -box_sides_half[2]
            delx[m_z_lo] += box[2]

            m_z_hi = delx[:, 2] > box_sides_half[2]
            delx[m_z_hi] -= box[2]