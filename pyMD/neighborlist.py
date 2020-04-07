import numpy as np


class NeighborList:
    # HARD CODED
    periodic_x = True
    periodic_y = True
    periodic_z = True
    triclinic = False

    IMG_x_lo = 1 << 0
    IMG_x_hi = 1 << 1
    IMG_y_lo = 1 << 2
    IMG_y_hi = 1 << 3
    IMG_z_lo = 1 << 4
    IMG_z_hi = 1 << 5

    def __init__(self,
                 positions,
                 box,
                 pair_cutoffs=[],
                 pagesize=2000,
                 skin=2.0,
                 every=10,
                 newton=True,
                 ):
        self.pagesize = pagesize
        self.skin = skin
        self.every = every

        self.n_atoms = len(positions)
        self.positions = positions
        self.box = box
        self.box_sides = np.linalg.norm(box, axis=1)
        self.box_sides_half = self.box_sides / 2.0

        self.pair_cutoffs = pair_cutoffs

        self.init_neighborlist(self.pair_cutoffs)

        if newton:
            self.build_neighborlist = self.build_neighborlist_newton
        else:
            raise NotImplementedError("Only neighborlist with newton on are implemented")

    def init_neighborlist(self, pair_cutoffs):
        self.nn_cutoff = max(pair_cutoffs) + self.skin
        self.nn_cutsq = self.nn_cutoff ** 2

        if self.triclinic:
            raise NotImplementedError()
        else:
            np.testing.assert_almost_equal(self.box_sides, np.sum(self.box, axis=1), err_msg='box is triclinic')
        if np.any(self.box_sides < 2 * self.nn_cutoff):
            raise UserWarning("Box is to small")

        self.n_locals = np.zeros(self.n_atoms, dtype=np.int32)
        self.neighborlist = np.zeros((self.n_atoms, self.pagesize), dtype=np.int32)
        self.ghost_images = np.zeros((self.n_atoms, self.pagesize), dtype=np.byte)

    def build_neighborlist_newton(self):
        positions = self.positions
        self.ghost_images[:] = 0

        for i, x in enumerate(positions[:-1]):

            delx = x - positions[i + 1:]

            if self.triclinic:
                raise NotImplementedError("...")
            else:
                if self.periodic_x:
                    m_x_lo = delx[:, 0] < -self.box_sides_half[0]
                    delx[m_x_lo] += self.box[0]

                    m_x_hi = delx[:, 0] > self.box_sides_half[0]
                    delx[m_x_hi] -= self.box[0]

                if self.periodic_y:
                    m_y_lo = delx[:, 1] < -self.box_sides_half[1]
                    delx[m_y_lo] += self.box[1]

                    m_y_hi = delx[:, 1] > self.box_sides_half[1]
                    delx[m_y_hi] -= self.box[1]

                if self.periodic_z:
                    m_z_lo = delx[:, 2] < -self.box_sides_half[2]
                    delx[m_z_lo] += self.box[2]

                    m_z_hi = delx[:, 2] > self.box_sides_half[2]
                    delx[m_z_hi] -= self.box[2]

            rsq = np.sum(np.power(delx, 2), axis=1)
            mask = rsq < self.nn_cutsq
            neighbors = np.where(mask)[0]
            n_local = len(neighbors)
            if n_local > 0:
                self.neighborlist[i, :n_local] = neighbors + i + 1

                ghost_images = self.ghost_images[i, :n_local]
                if self.periodic_x:
                    np.bitwise_or(ghost_images, self.IMG_x_lo, out=ghost_images, where=m_x_lo[neighbors])
                    np.bitwise_or(ghost_images, self.IMG_x_hi, out=ghost_images, where=m_x_hi[neighbors])
                if self.periodic_y:
                    np.bitwise_or(ghost_images, self.IMG_y_lo, out=ghost_images, where=m_y_lo[neighbors])
                    np.bitwise_or(ghost_images, self.IMG_y_hi, out=ghost_images, where=m_y_hi[neighbors])
                if self.periodic_z:
                    np.bitwise_or(ghost_images, self.IMG_z_lo, out=ghost_images, where=m_z_lo[neighbors])
                    np.bitwise_or(ghost_images, self.IMG_z_hi, out=ghost_images, where=m_z_hi[neighbors])

            self.n_locals[i] = n_local

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

    @classmethod
    def convert_ghostimage_to_1d(self, image):
        ghost_images = np.zeros(len(image), dtype=np.byte)

        np.bitwise_or(ghost_images, self.IMG_x_lo, out=ghost_images, where=image[:, 0] < 0)
        np.bitwise_or(ghost_images, self.IMG_x_hi, out=ghost_images, where=image[:, 0] > 0)

        np.bitwise_or(ghost_images, self.IMG_y_lo, out=ghost_images, where=image[:, 1] < 0)
        np.bitwise_or(ghost_images, self.IMG_y_hi, out=ghost_images, where=image[:, 1] > 0)

        np.bitwise_or(ghost_images, self.IMG_z_lo, out=ghost_images, where=image[:, 2] < 0)
        np.bitwise_or(ghost_images, self.IMG_z_hi, out=ghost_images, where=image[:, 2] > 0)

        return ghost_images

    @classmethod
    def convert_ghostimage_to_3d(self, ghost_images):
        images_out = np.zeros((ghost_images.size, 3), dtype=np.int8)
        images_out[(ghost_images & self.IMG_x_lo).astype(bool), 0] = -1
        images_out[(ghost_images & self.IMG_x_hi).astype(bool), 0] = 1

        images_out[(ghost_images & self.IMG_y_lo).astype(bool), 1] = -1
        images_out[(ghost_images & self.IMG_y_hi).astype(bool), 1] = 1

        images_out[(ghost_images & self.IMG_z_lo).astype(bool), 2] = -1
        images_out[(ghost_images & self.IMG_z_hi).astype(bool), 2] = 1

        return images_out