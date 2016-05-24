import numpy as np

import  collision.Collision as spatial
from collision.mesh import Mesh, icosphere
from collision import math
import collision.visualize.spheres

import itertools


class Actor(object):

    def collect_forces(self):
        self.force[:] = 0


class ParticleActor(Actor):
    def __init__(self, position, scale):
        self.position = np.ascontiguousarray(position.astype(np.float32))
        self.velocity = np.zeros_like(position)
        self.scale = scale
        self.mass = np.ones_like(self.position[:, 0]) * 1e-4

        self.force = np.zeros_like(self.position)

        self.gravity = np.array([0, 0, -5])

        self.stencil_prepare()
        self.collision_prepare()

    def stencil_prepare(self):
        stencil = [-1, 0, 1]
        stencil = itertools.product(*[stencil] * 3)
        self.stencil = np.array(list(stencil)).astype(np.int32)


    def collision_prepare(self):
        # self.spatial_grid = self.spatial_grid.update(self.position)
        # bounds = np.array([[-10]*3, [10]*3]).astype(np.float32)
        self.spec = spatial.Spec3d(self.position, self.scale)
        self.offsets = self.spec.stencil(self.stencil).astype(np.int32)
        self.spatial_grid = spatial.Grid3d(self.spec, self.position, self.offsets)

    def integrate(self, dt):
        self.velocity += self.force * dt / self.mass[:, None]
        self.position += self.velocity * dt

    def collect_forces(self):
        super(ParticleActor, self).collect_forces()
        self.force += self.gravity * self.mass[:, None]

    def init_visual(self, scene, parent):
        self.visual = collision.visualize.spheres.Spheres(self.scale/2, parent=parent)
        self.visual.set_data(self.position)


class HardParticleActor(ParticleActor):
    def collision_self(self):
        pairs = self.spatial_grid.get_pairs()
        s, e = pairs.T
        diff = self.position[s] - self.position[e]
        dist = np.linalg.norm(diff, axis=1)
        normal = diff / dist[:, None]
        depth = self.scale - dist
        stiffness = 2e0
        force = normal * depth[:, None] * stiffness
        np.add.at(self.force, s, force)
        np.add.at(self.force, e, -force)


class FluidParticleActor(ParticleActor):
    """
    http://www.ligum.umontreal.ca/Clavet-2005-PVFS/pvfs.pdf
    """

    def __init__(self, position, scale):
        self.scale = scale
        self.k_far = 0.004
        self.k_near = 0.01
        self.neutral_density = 10

        self.position = np.ascontiguousarray(position.astype(np.float32))
        self.previous_position = self.position.copy()
        self.velocity = np.zeros_like(self.position)
        self.mass = np.ones_like(self.position[:, 0]) * 1e-4

        self.force = np.zeros_like(self.position)

        self.far = np.zeros_like(self.mass)
        self.near = np.zeros_like(self.mass)

        self.gravity = np.array([0, 0, -5])

        self.stencil_prepare()
        self.collision_prepare()

    def collision_self(self):
        pass

    def integrate(self, dt):
        dt=dt*100
        self.velocity = (self.position - self.previous_position) / dt
        self.velocity += self.force * dt
        self.previous_position = self.position.copy()
        self.position += self.velocity * dt
        self.relaxation(dt)

    def relaxation(self, dt):
        pairs = self.spatial_grid.get_pairs()
        s, e = pairs.T
        diff = self.position[s] - self.position[e]
        dist = np.linalg.norm(diff, axis=1)

        weight = 1 - dist / self.scale

        self.far[:] = 0
        far  = weight * weight
        np.add.at(self.far, s, far)
        np.add.at(self.far, e, far)

        self.near[:] = 0
        near = far * weight
        np.add.at(self.near, s, near)
        np.add.at(self.near, e, near)

        p_far  = (self.far - self.neutral_density) * self.k_far
        p_near = self.near * self.k_near

        displacement = (p_far[pairs].mean(axis=1) * weight + p_near[pairs].mean(axis=1) * far) * (dt**2 / 2)
        displacement = (displacement / dist)[:, None] * diff
        np.add.at(self.position, s, displacement)
        np.add.at(self.position, e, -displacement)


class MeshActor(Actor):
    def __init__(self, mesh):
        assert mesh.volume() > 0
        assert mesh.is_orientated()

        self.mesh = mesh

        self.position = self.mesh.vertices      # alias mesh vertices
        self.velocity = np.zeros_like(self.position)
        self.force = np.zeros_like(self.position)
        self.length_scale = float(mesh.edge_lengths().mean())

        self.gravity = np.array([0, 0, -5])

        self.spatial_grid = self.get_spatial_grid()

    def integrate(self, dt):
        pass

    def collision_prepare(self):
        # self.spatial_grid = self.spatial_grid.update(self.position)
        self.spatial_grid = self.get_spatial_grid()
        self.spatial_mesh = self.get_spatial_mesh()

    def get_spatial_grid(self):
        stencil = np.zeros((1), np.int32)
        return spatial.Grid3d(
            spatial.Spec3d(self.position, self.length_scale * 2),
            self.position,
            stencil
        )
    def get_spatial_mesh(self):
        return spatial.Mesh(
            self.position,
            self.mesh.vertex_normals(),
            np.ascontiguousarray(self.mesh.faces[:, ::+1]),
            self.length_scale, 0
        )

    def permute(self):
        """apply permutation of previously computed grid to all state variables"""

    def init_visual(self, scene, parent):
        self.visual = scene.visuals.Mesh(
            self.mesh.vertices,
            self.mesh.faces[:, ::+1],
            shading='flat',
            parent=parent)


class StaticActor(MeshActor):
    def __init__(self, mesh):
        super(StaticActor, self).__init__(mesh)

        self.spatial_grid = self.get_spatial_grid()
        self.spatial_mesh = self.get_spatial_mesh()
        # self.get_spatial_grid = lambda : grid
        # self.get_spatial_mesh = lambda : mesh

    def collision_prepare(self):
        pass

class RigidActor(MeshActor):
    """integrates all external forces into angular momentum quaternion"""
    pass

class Nonridid(MeshActor):

    def __init__(self, mesh, elasticity, damping):
        super(Nonridid, self).__init__(mesh)
        self.vertex_incidence = self.mesh.compute_vertex_incidence()  # ExV
        self.elasticity = elasticity
        self.damping = damping

        self.velocity = np.zeros_like(self.position)
        self.mass = self.mesh.vertex_areas()

        self.compute_rest_length()

    def compute_rest_length(self):
        edges = self.vertex_incidence * self.mesh.vertices
        self.rest_length = np.linalg.norm(edges, axis=1)

    def compute_edge_forces(self):
        edges = self.vertex_incidence * self.position
        compression = (self.rest_length - np.linalg.norm(edges, axis=1)) / self.rest_length
        dir = math.normalize(edges)
        relative_velocity = self.vertex_incidence * self.velocity
        relative_velocity = math.dot(relative_velocity, dir)
        force = compression * self.elasticity - relative_velocity * self.damping
        return self.vertex_incidence.T * (dir * force[:, None])

    def integrate(self, dt):
        self.velocity += self.force * dt / self.mass[:, None]
        self.position += self.velocity * dt

    def collect_forces(self):
        super(Nonridid, self).collect_forces()
        self.force += self.gravity * self.mass[:, None]
        self.force += self.compute_edge_forces()


class Balloon(Nonridid):

    def __init__(self, mesh, elasticity, damping, compressibility):
        super(Balloon, self).__init__(mesh, elasticity, damping)
        self.compressibility = compressibility
        self.compute_rest_volume()

    def compute_rest_volume(self):
        self.rest_volume = self.mesh.volume() * 3

    def compute_volume_forces(self):
        compression = (self.rest_volume - self.mesh.volume()) / self.rest_volume
        volume_gradient = self.mesh.vertex_volume_gradient()
        return volume_gradient * (compression * self.compressibility)

    def collect_forces(self):
        super(Balloon, self).collect_forces()
        self.force += self.compute_volume_forces()


class Scene(object):
    def __init__(self, actors):
        self.actors = actors

    def integrate(self, dt):
        """integrate the state of the scene with a timestep dt"""
        for a in self.actors:
            a.collect_forces()
        self.collide()
        for a in self.actors:
            a.integrate(dt)

    def collide(self):
        for a in self.actors:
            a.collision_prepare()

        for i, ai in enumerate(self.actors):
            if isinstance(ai, ParticleActor):
                ai.collision_self()

            for j, aj in enumerate(self.actors):
                if isinstance(aj, ParticleActor):
                    continue
                if i == j: continue
                info = spatial.Info(ai.spatial_grid, aj.spatial_mesh, i==j)
                mask = info.triangle != -1
                active = np.flatnonzero(mask)   # active vertex idx
                if np.any(mask):# and not isinstance(ai, StaticActor):
                    triangle = info.triangle[active]
                    bary = info.bary[active]

                    velocity = aj.velocity[aj.mesh.faces[triangle]]
                    velocity = np.einsum('vtc,vt->vc', velocity, bary)
                    relative_velocity = ai.velocity[active] - velocity

                    friction = 1e-2
                    stiffness = 1e1
                    force = info.depth[active][:, None] * info.normal[active] * stiffness - relative_velocity * friction
                    assert not np.any(np.isnan(force))

                    np.add.at(ai.force, active, force)
                    corners = aj.mesh.faces[triangle]
                    for i in range(3):
                        np.add.at(aj.force, corners[:, i], -bary[:, [i]] * force)
                        # aj.force[corners[:,i]] -= info.bary[:, i] * force

                    # aj.force
                    # print(info.triangle[mask])
                    # print(info.depth[mask])

    def plot(self):

        from vispy import app, scene, io

        # Prepare canvas
        canvas = scene.SceneCanvas(keys='interactive', size=(800, 600), show=True)
        canvas.measure_fps()

        # Set up a viewbox to display the image with interactive pan/zoom
        view = canvas.central_widget.add_view()

        for a in self.actors:
            a.init_visual(scene, view.scene)

        fov = 60.
        cam1 = scene.cameras.FlyCamera(parent=view.scene, fov=fov, name='Fly')
        view.camera = cam1

        # def on_resize(self, event):
        #     176  # setup the new viewport
        #     gloo.set_viewport(0, 0, *event.physical_size)
        #     w, h = event.size
        #     self.projection = perspective(45.0, w / float(h), 1.0, 1000.0)
        #     self.program['u_projection'] = self.projection

        dt = 0.002
        def update(event):
            # update the simulation
            for i in range(1):
                self.integrate(dt)
            # upload new state to gpu
            for i in range(len(self.actors)):
                actor = self.actors[i]
                if isinstance(actor, MeshActor):
                    if not isinstance(actor, StaticActor):
                        actor.visual.set_data(vertices=actor.position[actor.mesh.faces])
                if isinstance(actor, ParticleActor):
                    actor.visual.set_data(actor.position)

        timer = app.Timer(interval=dt, connect=update)

        timer.start()
        app.run()


if __name__=='__main__':
    turtle = Mesh.load_stl('part0.stl')
    # normalize orientation
    u, s, v = np.linalg.svd(turtle.vertices, full_matrices=0)
    # v[:, [1, 2, 0]] on other machine
    turtle.vertices = turtle.vertices.dot(v) * 3 + np.array([0,0,0], np.float32)
    # turtle.faces = turtle.faces[:, ::-1]

    ico = icosphere(0.1, refinement=3)
    ball = lambda p : icosphere(0.1, p, 3)
    e = 1.5e-1
    c = 1e2
    d = .02

    grid_points = np.mgrid[0:6, 0:6, 0:40].reshape(3, -1).T * 0.05
    actors = [StaticActor(turtle),
              Balloon(ball([0,0,-0.8]), e, d, c),
              # Balloon(ball([0,0.2,-0.6]), e, d, c),
              # Balloon(ball([0,0,-0.2]), e, d, c),
              # Balloon(ball([0,0,0]), e, d, c),
              HardParticleActor(np.random.rand(1000, 3) * [2,2,5] + [[0,-0.5,-2.8]], scale=0.1)
              # FluidParticleActor(grid_points+ [[0, -0.5, -1.8]], scale=0.1)
              ]

    scene = Scene(actors)
    scene.plot()