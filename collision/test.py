import numpy as np


import Collision


points = np.random.rand(30, 3).astype(np.float32)
lengthscale = 0.5
grid = Collision.VertexGridHash(points, lengthscale)
print(grid.indices)
print(grid.pivots)
print(grid.n_buckets)
# for i,id in zip(grid.indices, grid.cell_ids):
#     print(i, id)