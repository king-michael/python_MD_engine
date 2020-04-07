import sys
sys.path.insert(0, '..')
import numpy as np
from pyMD.neighborlist import NeighborList


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


position=np.array([[0,0,0],
                   [8,0,0]])
box = np.array([
    [10,0,0],
    [0,10,0],
    [0,0,10]
])
neigh = NeighborList(
    position, box, [2.0]
    )

neigh.build_neighborlist()

img = neigh.convert_ghostimage_to_3d(neigh.ghost_images[0,:1])

np.testing.assert_array_equal(np.array([[-1,  0,  0]], dtype=np.int8), img)