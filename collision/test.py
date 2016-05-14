import numpy as np


import Collision


points = np.random.rand(30, 3).astype(np.float32)
lengthscale = 0.5
grid = Collision.Grid3d(points, lengthscale)
print(grid.permutation)
print(grid.pivots)#[:grid.n_buckets+1])
print(grid.n_buckets)
for i,c in zip(grid.permutation, grid.cells):
    print(i, c)