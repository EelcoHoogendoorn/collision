import numpy as np


import Collision


points = np.random.rand(30, 2).astype(np.float32)
lengthscale = 0.5
grid = Collision.Grid2d(points, lengthscale)
print(grid.indices)
print(grid.pivots[:grid.n_buckets+1])
print(grid.n_buckets)
# for i,id in zip(grid.indices, grid.cell_ids):
#     print(i, id)