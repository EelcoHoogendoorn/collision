import numpy as np


import Collision
import time

def test_performance():
    # warmup run for mem allocation
    lengthscale = 0.03
    points = np.random.rand(3000000, 3).astype(np.float32)
    grid = Collision.Grid3d(points, lengthscale)

    # run on unsorted data
    points = np.random.rand(3000000, 3).astype(np.float32)
    start = time.clock()
    grid = Collision.Grid3d(points, lengthscale)
    unsorted = time.clock() - start

    # sort data
    points = points[grid.permutation]

    # run on sorted data
    start = time.clock()
    grid = Collision.Grid3d(points, lengthscale)
    sorted = time.clock() - start

    print(unsorted, sorted)


def test_basic():
    lengthscale = 0.4
    points = np.random.rand(30, 2).astype(np.float32)
    grid = Collision.Grid2d(points, lengthscale)
    print(grid.permutation)
    print(grid.pivots)#[:grid.n_buckets+1])
    print(grid.n_buckets)
    print(grid.cells[grid.permutation])


def test_mesh():
    pass

test_basic()
test_performance()