import numpy as np


import Collision


points = np.random.rand(100, 3).astype(np.float32)
lengthscale = 0.1
grid = Collision.VertexGridHash(points, lengthscale)
print(grid.pivots)