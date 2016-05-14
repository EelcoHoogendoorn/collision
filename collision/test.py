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
    from collision.mesh import Mesh
    mesh = Mesh.load_stl(r'part0.stl')
    cmesh = Collision.Mesh(mesh.vertices, mesh.vertex_normals(), mesh.faces, float(mesh.edge_lengths().mean()))
    boxes = cmesh.boxes.reshape(-1, 3, 2).transpose(0, 2, 1)
    for box in boxes[:10]:
        print(box)
    vispy_plot(mesh)



def vispy_plot(mesh):

    from vispy import app, scene, io

    # Prepare canvas
    canvas = scene.SceneCanvas(keys='interactive', size=(800, 600), show=True)
    canvas.measure_fps()

    # Set up a viewbox to display the image with interactive pan/zoom
    view = canvas.central_widget.add_view()

    meshvis = scene.visuals.Mesh(
        mesh.vertices * 100,
        mesh.faces[:, ::-1],
        shading='flat',
        parent=view.scene)

    # Create three cameras (Fly, Turntable and Arcball)
    fov = 60.
    cam1 = scene.cameras.FlyCamera(parent=view.scene, fov=fov, name='Fly')
    view.camera = cam1

    app.run()


test_performance()
