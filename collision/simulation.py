import numpy as np

import  collision.Collision as spatial
from collision.mesh import Mesh, icosphere
from collision import math


class Actor(object):
    def __init__(self, mesh):
        self.mesh = mesh

        self.position = self.mesh.vertices      # alias mesh vertices
        self.force = np.zeros_like(self.position)
        self.length_scale = float(mesh.edge_lengths().mean())

        self.gravity = np.array([0, 0, -1])

    def integrate(self, dt):
        pass

    def get_spatial_grid(self):
        return spatial.Grid3d(
            self.position,
            self.length_scale
        )
    def get_spatial_mesh(self):
        return spatial.Mesh(
            self.position,
            self.mesh.vertex_normals(),
            self.mesh.faces,
            self.length_scale
        )
    def spatial_collide(self, other):
        spatial.Info(self.get_spatial_grid(), other.get_spatial_mesh())
    def permute(self):
        """apply permutation of previously computed grid to all state variables"""

    def compute_forces(self):
        self.forces[:] = self.gravity

class StaticActor(Actor):
    pass

class RigidActor(Actor):
    """integrates all external forces into angular momentum quaternion"""
    pass

class Nonridid(Actor):

    def __init__(self, mesh, elasticity):
        super(Nonridid, self).__init__(mesh)
        self.vertex_incidence = inc = self.mesh.compute_vertex_incidence()  # ExV
        self.elasticity = elasticity

        self.velocity = np.zeros_like(self.position)
        self.mass = self.mesh.vertex_areas()

        self.compute_rest_length()

    def compute_rest_length(self):
        edges = self.vertex_incidence * self.mesh.vertices
        self.rest_length = np.linalg.norm(edges)

    def compute_edge_forces(self):
        edges = self.vertex_incidence * self.mesh.vertices
        diff = np.linalg.norm(edges) - self.rest_length
        dir = math.normalized(edges)
        force = diff * self.elasticity
        return self.vertex_incidence.T * (dir * force)

    def integrate(self, dt):
        self.velocity += self.force * dt / self.mass
        self.position += self.velocity * dt

    def collect_forces(self):
        super(Nonridid, self).collect_forces()
        self.force += self.compute_edge_forces()


class Balloon(Nonridid):

    def __init__(self, mesh, elasticity, compressibility):
        super(Balloon, self).__init__(mesh, elasticity)
        self.compressibility = compressibility
        self.compute_rest_length()

    def compute_rest_volume(self):
        self.rest_volume = self.mesh.volume()

    def compute_volume_forces(self):
        diff = self.mesh.volume() - self.rest_volume
        normals = self.mesh.vertex_normals() * self.mass    # should actually be non-normalized normals...
        return normals * (diff * self.compressibility)

    def collect_forces(self):
        super(Balloon, self).collect_forces()
        self.force += self.compute_volume_forces()


class Scene(object):
    def __init__(self, actors):
        self.actors = actors

    def integrate(self, dt):
        for a in self.actors:
            a.collect_forces()
        for a in self.actors:
            a.integrate(dt)

    def plot(self):

        from vispy import app, scene, io

        # class Canvas(app.Canvas):
        #     def __init__(self):
        #         app.Canvas.__init__(self, keys='interactive', size=(800, 600))
        #
        #     def on_draw(self, ev):
        #         gloo.set_viewport(0, 0, *self.physical_size)
        #         gloo.clear(color='black', depth=True)
        #
        #         for mesh in self.meshes:
        #             mesh.draw()

        # Prepare canvas
        canvas = scene.SceneCanvas(keys='interactive', size=(800, 600), show=True)
        canvas.measure_fps()

        # Set up a viewbox to display the image with interactive pan/zoom
        view = canvas.central_widget.add_view()

        vis_meshes = [scene.visuals.Mesh(
            actor.mesh.vertices * 100,
            actor.mesh.faces[:, ::-1],
            shading='flat',
            parent=view.scene) for actor in self.actors]

        fov = 60.
        cam1 = scene.cameras.FlyCamera(parent=view.scene, fov=fov, name='Fly')
        view.camera = cam1

        app.run()


if __name__=='__main__':
    turtle = Mesh.load_stl('part0.stl')
    # normalize orientation
    u, s, v = np.linalg.svd(turtle.vertices, full_matrices=0)
    turtle.vertices = turtle.vertices.dot(v)
    turtle.faces = turtle.faces[:, ::-1]

    ico = icosphere(0.1, refinement=3)
    ball = lambda p : icosphere(0.1, p, 3)
    e = 1
    c = 1
    actors = [StaticActor(turtle),
              Balloon(ball([0,0,0]), e, c),
              Balloon(ball([0,0,1]), e, c)]

    scene = Scene(actors)
    scene.plot()