import numpy as np


import Collision


points = np.random.rand(48, 3).astype(np.float32)
lengthscale = 0.5
grid = Collision.VertexGridHash(points, lengthscale)
print(grid.indices)
print(grid.pivots)
print(points.min(), points.max())
