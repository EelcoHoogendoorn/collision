import numpy as np

import  collision.Collision as spatial
from collision.mesh import Mesh, icosphere

import time

def test_performance():
    # warmup run for mem allocation
    lengthscale = 0.03
    points = np.random.rand(300000, 3).astype(np.float32)
    grid = spatial.Grid3d(points, (lengthscale))

    # run on unsorted data
    points = np.random.rand(300000, 3).astype(np.float32)
    start = time.clock()
    grid = spatial.Grid3d(points, (lengthscale))
    unsorted = time.clock() - start

    # sort data
    points = points[grid.permutation]

    # run on sorted data
    start = time.clock()
    grid = spatial.Grid3d(points, lengthscale)
    sorted = time.clock() - start

    print(unsorted, sorted)


def test_basic():
    lengthscale = 0.4
    points = np.random.rand(30, 2).astype(np.float32)
    grid = spatial.Grid2d(points, lengthscale)
    print(grid.permutation)
    print(grid.pivots)#[:grid.n_buckets+1])
    print(grid.n_buckets)
    print(grid.cells[grid.permutation])


def test_mesh():
    meshes = [icosphere(0.5, refinement=3), icosphere(0.5, [1, 0, 0], refinement=3)]
    normals = [mesh.vertex_normals() for mesh in meshes]
    cmeshes = [spatial.Mesh(
                    mesh.vertices,
                    mesh.vertex_normals(),
                    mesh.faces,
                    float(mesh.edge_lengths().mean()))
               for mesh in meshes]
    # boxes = cmesh.boxes.reshape(-1, 3, 2).transpose(0, 2, 1)
    # for box in boxes[:10]:
    #     print(box)
    # mesh = icosphere(0.5, refinement=3)
    vispy_plot(meshes)


def test_collision():
    meshes = [icosphere(0.5, refinement=3), icosphere(0.55, [1, 0, 0], refinement=3)]
    lengthscale = 0.5
    grids = [spatial.Grid3d(m.vertices, lengthscale) for m in meshes]
    ctmeshes = [spatial.Mesh(
        m.vertices,
        m.vertex_normals(),
        m.faces,
        float(m.edge_lengths().mean()))
                for m in meshes]

    for i, mi in enumerate(meshes):
        for j, mj in enumerate(meshes):
            if i == j: continue
            print (i,j)
            info = spatial.Info(grids[i], ctmeshes[j], False)
            mask = info.triangle != -1
            print(info.triangle[mask])
            print(info.depth[mask])


def vispy_plot(meshes):

    from vispy import app, scene, io

    # Prepare canvas
    canvas = scene.SceneCanvas(keys='interactive', size=(800, 600), show=True)
    canvas.measure_fps()

    # Set up a viewbox to display the image with interactive pan/zoom
    view = canvas.central_widget.add_view()

    vis_meshes = [scene.visuals.Mesh(
        mesh.vertices * 100,
        mesh.faces[:, ::-1],
        shading='flat',
        parent=view.scene) for mesh in meshes]

    fov = 60.
    cam1 = scene.cameras.FlyCamera(parent=view.scene, fov=fov, name='Fly')
    view.camera = cam1

    app.run()


test_collision()
