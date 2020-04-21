import sys

sys.path.insert(0, '..')
import numpy as np
from pyMD.neighborlist import NeighborList


def test_convert_ghostimage_to_3d():
    position = np.array([[0, 0, 0],
                         [8, 0, 0]])
    box = np.array([
        [10, 0, 0],
        [0, 10, 0],
        [0, 0, 10]
    ])
    neigh = NeighborList(
        position, box, [2.0]
    )

    neigh.build_neighborlist()

    img = neigh.convert_ghostimage_to_3d(neigh.ghost_images[0, :1])

    np.testing.assert_array_equal(np.array([[-1, 0, 0]], dtype=np.int8), img)


def test_convert_ghostimage_to_1d():
    image = np.array([
        [1, 0, 0],
        [1, 1, 0],
        [1, -1, -1],
        [0, 0, 1],
    ])

    img = NeighborList.convert_ghostimage_to_3d(
        NeighborList.convert_ghostimage_to_1d(image)
    )

    np.testing.assert_array_equal(image, img)


def test_build_neighborlist_periodic_boundaries(verbose=False):

    cutoff = 2.0
    box_size = 10.0
    min_dist = 0.1

    if verbose:
        print("run : test_build_neighborlist")
        print(" - cutoff   : {}".format(cutoff))
        print(" - box_size : {}".format(box_size))
        print(" - min_dist : {}".format(min_dist))

    position = np.array([[0, 0, 0],
                         [0.1, 0, 0]])
    box = np.array([
        [box_size, 0, 0],
        [0, box_size, 0],
        [0, 0, box_size]
    ])
    neigh = NeighborList(
        position, box, [cutoff], skin=0.0
    )

    ref_positions = np.linspace(min_dist, box_size-min_dist, num=20)

    ref_nlocals = np.zeros(ref_positions.size, dtype=np.int32)
    ref_nlocals[ref_positions <= cutoff] = 1
    ref_nlocals[ref_positions >= box_size-cutoff] = 1

    ref_ghostimages_x = np.zeros(ref_positions.size, dtype=np.byte)
    ref_ghostimages_x[ref_positions >= box_size-cutoff] = 1

    test_nlocals = np.zeros(ref_positions.size, dtype=np.int32)
    test_ghost_images_x = np.zeros(ref_positions.size, dtype=np.byte)

    for i in range(ref_positions.size):
        neigh.positions[1,0] = ref_positions[i]
        neigh.build_neighborlist()
        test_nlocals[i] = neigh.n_locals[0]
        if neigh.n_locals[0]:
            test_ghost_images_x[i] = neigh.ghost_images[0, :neigh.n_locals[0]]

    np.testing.assert_equal(ref_nlocals, test_nlocals, err_msg='n_locals differ.')
    np.testing.assert_equal(ref_ghostimages_x, test_ghost_images_x, err_msg='ghost_images are wrong')

    if verbose:
        print('n_locals (ref, test)')
        print(ref_nlocals)
        print(test_nlocals)

        print('ghost_images (ref, test)')
        print(ref_ghostimages_x)
        print(test_ghost_images_x)


if __name__ == '__main__':
    test_convert_ghostimage_to_3d()
    test_convert_ghostimage_to_1d()

    test_build_neighborlist_periodic_boundaries(verbose=True)
