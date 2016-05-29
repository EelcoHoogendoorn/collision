
import time

import numpy as np

import numpy_indexed as npi

import  collision.Collision as spatial
from collision.mesh import Mesh, icosphere



def test_basic():
    lengthscale = 0.5
    points = np.random.rand(30, 2).astype(np.float32)
    grid = spatial.Grid2d(points, lengthscale)
    print(grid.permutation)
    print(grid.pivots)#[:grid.n_buckets+1])
    print(grid.n_buckets)
    print(grid.cells[grid.permutation])


def test_performance():
    # warmup run for mem allocation
    lengthscale = 0.03
    n = 300000
    points = np.random.rand(n, 3).astype(np.float32)
    spec = spatial.Spec3d(points, lengthscale)

    start = time.clock()
    grid = spatial.Grid3d(spec, points)
    warmup = time.clock() - start

    # run on unsorted data
    points = np.random.rand(n, 3).astype(np.float32)
    start = time.clock()
    grid = spatial.Grid3d(spec, points)
    unsorted = time.clock() - start

    # update from existing
    start = time.clock()
    grid = grid.update(points)
    update = time.clock() - start

    # sort data
    points = points[grid.permutation]

    # run on sorted data
    start = time.clock()
    grid = spatial.Grid3d(spec, points)
    sorted = time.clock() - start

    print(warmup, unsorted, update, sorted)


def test_mesh():
    mesh = icosphere(0.5, refinement=3)
    edge_length = float(mesh.edge_lengths().mean())
    normals = mesh.vertex_normals()
    cmesh = spatial.Mesh(
                    mesh.vertices,
                    mesh.vertex_normals(),
                    mesh.faces,
                    edge_length, edge_length)
    boxes = cmesh.boxes.reshape(-1, 3, 2).transpose(0, 2, 1)
    for box in boxes[:10]:
        print(box)


def test_collision():
    meshes = [icosphere(0.5, refinement=3), icosphere(0.5, [0.95, 0, 0], refinement=3)]
    lengthscale = 0.2
    grids = [spatial.Grid3d(m.vertices, lengthscale) for m in meshes]
    ctmeshes = [spatial.Mesh(
        m.vertices,
        m.vertex_normals(),
        m.faces,
        0.1, 0)
            for m in meshes
    ]

    for i, mi in enumerate(meshes):
        for j, mj in enumerate(meshes):
            print (i,j)
            info = spatial.Info(grids[i], ctmeshes[j], i == j)
            mask = info.triangle != -1
            print(info.triangle[mask])
            print(info.depth[mask])


def GridSpec(points, scale):
    spectypemap = {2: spatial.Spec2d, 3: spatial.Spec3d}
    n_points, n_dim = points.shape
    return spectypemap[n_dim](points, scale)

def PointGrid(spec, points, offset):
    gridtypemap = {2: spatial.Grid2d, 3: spatial.Grid3d}
    n_points, n_dim = points.shape
    return gridtypemap[n_dim](spec, points, offset)

def BoxGrid(spec, boxes):
    gridtypemap = {3: spatial.BoxGrid3d}
    n_boxes, _, n_dim = boxes.shape
    return gridtypemap[n_dim](spec, boxes.reshape(n_boxes, -1))


def test_point_performance():
    import itertools
    stencil = [-1, 0, 1]
    ndim = 3
    stencil = itertools.product(*[stencil]*ndim)
    stencil = np.array(list(stencil)).astype(np.int32)
    print (stencil.strides)
    points = np.random.rand(4000, ndim).astype(np.float32)
    scale = 0.1

    spec = GridSpec(points, float(scale))
    offsets = spec.stencil(stencil).astype(np.int32)
    print(offsets)

    start = time.clock()
    grid = PointGrid(spec, points, offsets)
    construction = time.clock() - start

    start = time.clock()
    grid = PointGrid(spec, points[grid.permutation()], offsets)
    construction = time.clock() - start

    start = time.clock()
    pairs = grid.pairs()
    pairing = time.clock() - start

    from scipy.spatial import cKDTree
    start = time.clock()
    tree = cKDTree(points)
    tree_construction = time.clock() - start
    start = time.clock()
    tree_pairs = tree.query_pairs(scale, output_type='ndarray')
    tree_pairing = time.clock() - start

    print(points.shape, pairs.shape, len(pairs) / len(points))
    print(construction, pairing, pairing / construction)
    print(tree_construction, tree_pairing)


def test_point_correctness():
    import itertools
    stencil = [-1, 0, 1]
    ndim = 3
    n = 2000
    stencil = itertools.product(*[stencil]*ndim)
    stencil = np.array(list(stencil)).astype(np.int32)

    points = (np.random.rand(n, ndim) * [1, 2, 3]).astype(np.float32)
    scale = 0.1

    spec = GridSpec(points, float(scale))
    offsets = spec.stencil(stencil).astype(np.int32)
    grid = PointGrid(spec, points, offsets)

    pairs = grid.pairs()

    from scipy.spatial import cKDTree
    tree = cKDTree(points)
    tree_pairs = tree.query_pairs(scale, output_type='ndarray')
    print(tree_pairs)
    print(pairs)

    assert np.alltrue(npi.unique(tree_pairs) == npi.unique(np.sort(pairs, axis=1)))


def test_box_grid():

    n = 40
    ndim = 3
    points = (np.random.rand(n, ndim) * [1, 2, 3]).astype(np.float32)
    scale = 0.1
    bounds = np.array([[0, 0, 0], [1, 2, 3]])

    boxes = points[:10, None, :]
    boxes = np.concatenate([boxes, boxes + 0.2], axis=1)


    spec = GridSpec(bounds.astype(np.float32), float(scale))
    grid = BoxGrid(spec, boxes)
    print (grid.permutation)



test_box_grid()

