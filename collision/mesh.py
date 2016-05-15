
import numpy as np
import numpy_indexed as npi
import scipy.sparse
import scipy.sparse.linalg
import scipy.spatial

import collision.math

class PolyData(object):

    def merge(self, other):
        vertices = np.concatenate([self.vertices, other.vertices], axis=0)
        faces = np.concatenate([self.faces, other.faces + len(self.vertices)], axis=0)
        _, _idx, _inv = npi.unique(vertices, return_index=True, return_inverse=True)
        return type(self)(vertices[_idx], _inv[faces])

    def squeeze(self):
        """compact geometry description, removing unused vertices, and adjusting faces accordingly"""
        active, inv = np.unique(self.faces, return_inverse=True)
        return type(self)(self.vertices[active], inv.reshape(self.faces.shape))

    @staticmethod
    def order_edges(edges):
        return np.where((edges[:, 0] < edges[:, 1])[:, None], edges[:, ::+1], edges[:, ::-1])


class Mesh(PolyData):
    def __init__(self, vertices, faces):
        vertices = np.asarray(vertices, dtype=np.float32)
        faces = np.asarray(faces, dtype=np.int32)

        assert vertices.shape[1] == 3
        assert faces.shape[1] == 3

        self.vertices = vertices
        self.faces = faces

    def edges(self):
        """construct 3 edges for each triangle"""
        edges = self.faces[:, [[1, 2], [2, 0], [0, 1]]]
        return edges.reshape(-1, 2)

    def ordered_edges(self):
        return self.order_edges(self.edges())

    def compute_face_incidence(self):
        unsorted_edges = self.edges().reshape(-1, 2)
        sorted_edges = np.sort(unsorted_edges, axis=-1)

        unique_edges, edge_indices = npi.unique(sorted_edges, return_inverse=True)
        face_indices = np.arange(self.faces.size) // 3
        orientation = sorted_edges[:, 0] == unsorted_edges[:, 0]
        incidence = scipy.sparse.csr_matrix((orientation * 2 - 1, (edge_indices, face_indices)))
        return incidence, unique_edges

    def compute_vertex_incidence(self):
        unsorted_edges = self.edges().reshape(-1, 2)
        sorted_edges = np.sort(unsorted_edges, axis=-1)

        vertex_indices = npi.unique(sorted_edges)
        edge_indices = np.arange(vertex_indices.size) // 2
        orientation = vertex_indices == vertex_indices[:, 0:1]

        incidence = scipy.sparse.csr_matrix(
            ((orientation * 2 - 1).flatten(), (edge_indices, vertex_indices.flatten())))
        return incidence

    def is_orientated(self):
        return npi.all_unique(self.edges())

    def boundary(self):
        edges = self.edges()
        return edges[npi.multiplicity(self.order_edges(edges)) == 1]

    def is_orientated(self):
        return npi.all_unique(self.edges())

    def face_normals(self):
        edges = np.diff(self.vertices[self.faces], axis=1)
        return np.cross(edges[:, 1], edges[:, 0]) / 2

    def vertex_areas(self):
        """compute area associated with each vertex,
        whereby each face shares its area equally over all its component vertices

        Returns
        -------
        vertex_area : ndarray, [n_vertices], float
            vertex area per vertex
        """
        face_areas = np.linalg.norm(self.face_normals(), axis=-1)
        areas = np.zeros_like(self.vertices[:, 0])
        for i in range(3):
            np.add.at(areas, self.faces[:, i], face_areas)
        return areas / 3

    def vertex_normals(self):
        """
        Compute vertex normals, by averaging the weighted sum of its incident face normals.
        """
        assert self.is_orientated(), 'Vertex normals require an oriented segment'

        fnormals = self.face_normals()
        vnormals = np.zeros_like(self.vertices)
        # process the three corners of all faces one-by-one; would be nice to vectorize, but unsure if we can
        for i in range(3):
            np.add.at(vnormals, self.faces[:, i], fnormals)
        return collision.math.normalize(vnormals)

    def volume(self):
        return (self.face_normals() * self.face_centroids()).sum() / 3

    def compute_angles(self):
        """compute angles for each triangle-vertex"""
        edges = self.edges().reshape(-1, 3, 2)
        vecs = np.diff(self.vertices[edges], axis=2)[:, :, 0]
        vecs = collision.math.normalize(vecs)
        angles = np.arccos(-collision.math.dot(vecs[:, [1, 2, 0]], vecs[:, [2, 0, 1]]))
        assert np.allclose(angles.sum(axis=1), np.pi, rtol=1e-3)
        return angles

    def remap_edges(self, field):
        """given a quantity computed on each triangle-edge, sum the contributions from each adjecent triangle"""
        edges = self.edges().reshape(-1, 3, 2)
        sorted_edges = np.sort(edges, axis=-1)
        _, field = npi.group_by(sorted_edges.reshape(-1, 2)).sum(field.flatten())
        return field

    def hodge_edge(self):
        """compute edge hodge based on cotan formula; corresponds to circumcentric calcs"""
        cotan = 1 / np.tan(self.compute_angles())
        return self.remap_edges(cotan) / 2

    def laplacian_vertex(self):
        """compute cotan/area gradient based laplacian"""
        hodge = self.hodge_edge()
        hodge = scipy.sparse.dia_matrix((hodge, 0), shape=(len(hodge),) * 2)
        incidence = self.compute_vertex_incidence()
        return incidence.T * hodge * incidence

    def compute_gradient(self, field):
        """compute gradient of scalar function on vertices on faces"""
        normals = self.face_normals()
        face_area = np.linalg.norm(normals, axis=1)
        normals = collision.math.normalize(normals)

        edges = self.edges().reshape(-1, 3, 2)
        vecs = np.diff(self.vertices[edges], axis=2)[:, :, 0, :]
        gradient = (field[self.faces][:, :, None] * np.cross(normals[:, None, :], vecs)).sum(axis=1)
        return gradient / (2 * face_area[:, None])

    def compute_divergence(self, field):
        """compute divergence of vector field at faces at vertices"""
        edges = self.edges().reshape(-1, 3, 2)
        sorted_edges = np.sort(edges, axis=-1)
        vecs = np.diff(self.vertices[sorted_edges], axis=2)[:, :, 0, :]
        inner = collision.math.dot(vecs, field[:, None, :])
        cotan = 1 / np.tan(self.compute_angles())
        vertex_incidence = self.compute_vertex_incidence()
        return vertex_incidence.T * self.remap_edges(inner * cotan) / 2

    def face_centroids(self):
        return self.vertices[self.faces].mean(axis=1)

    def edge_lengths(self):
        edges = self.vertices[self.edges()]
        lengths = np.linalg.norm(edges[:, 1]- edges[:, 0], axis=1)
        return self.remap_edges(lengths) / 2

    def geodesic(self, seed, m=1):
        """Compute geodesic distance map

        Notes
        -----
        http://www.multires.caltech.edu/pubs/GeodesicsInHeat.pdf
        """
        laplacian = self.laplacian_vertex()
        mass = self.vertex_areas()
        t = self.edge_lengths().mean() ** 2 * m
        heat = lambda x : mass * x - laplacian * (x * t)
        operator = scipy.sparse.linalg.LinearOperator(shape=laplacian.shape, matvec=heat)

        diffused = scipy.sparse.linalg.minres(operator, seed.astype(np.float64), tol=1e-5)[0]
        gradient = -collision.math.normalize(self.compute_gradient(diffused))
        # self.plot(facevec=gradient)
        rhs = self.compute_divergence(gradient)
        phi = scipy.sparse.linalg.minres(laplacian, rhs)[0]
        return phi - phi.min()

    def save_STL(self, filename):
        """save a mesh to plain stl. vertex ordering is assumed to be correct"""
        header      = np.zeros(80, '<c')
        triangles   = np.array(len(self.faces), '<u4')
        dtype       = [('normal', '<f4', 3,),('vertex', '<f4', (3,3)), ('abc', '<u2', 1,)]
        data        = np.empty(triangles, dtype)

        data['abc']    = 0     #standard stl cruft
        data['vertex'] = self.vertices[self.faces]
        data['normal'] = self.normalize(self.face_normals())

        with open(filename, 'wb') as fh:
            header.   tofile(fh)
            triangles.tofile(fh)
            data.     tofile(fh)

    def load_stl(filename):
        dtype       = [('normal', '<f4', 3,),('vertex', '<f4', (3,3)), ('abc', '<u2', 1,)]

        with open(filename, 'rb') as fh:
            header    = np.fromfile(fh, '<c', 80)
            triangles = np.fromfile(fh, '<u4', 1)[0]
            data      = np.fromfile(fh, dtype, triangles)

        vertices, triangles = npi.unique(data['vertex'].reshape(-1, 3), return_inverse=True)
        return Mesh(vertices, triangles.reshape(-1, 3))


def icosahedron():
    """
    Generate a icosahedron.

    Returns
    -------
    mesh: Mesh
    """
    t = (1 + np.sqrt(5)) / 2
    norm_f = np.sqrt(t**2 + 1)

    a = np.array([-1, 1, -1, 1]) / norm_f
    b = np.array([t, t, -t, -t]) / norm_f
    c = np.zeros(4)

    vertices_part = np.array([a, b, c]).reshape((3, 4)).T
    vertices = np.array([np.roll(vertices_part, shift, axis=1) for shift in range(3)]).reshape((12, 3))

    return triangulate_convex(vertices)


def triangulate_convex(vertices):
    """triangulate a convex set of vertices

    Parameters
    ----------
    vertices : ndarray, [n, 3], float
        convex set of vertices

    Returns
    -------
    skcg.Mesh instance, with outward-pointing normals
    """
    hull = scipy.spatial.ConvexHull(vertices)
    return Mesh(vertices, hull.simplices)#.fix_orientation().fix_curvature()


def refine_sphere(mesh):
    """given a spherical mesh, insert a new vertex on every edge

    Parameters
    ----------
    mesh : skcg.Mesh instance

    Returns
    -------
    skcg.Mesh instance
    """
    vertices = mesh.vertices
    edges = npi.unique(np.sort(mesh.edges, axis=2).reshape(-1, 2))
    new_vertices = vertices[edges].mean(axis=1)
    new_vertices /= np.linalg.norm(new_vertices, axis=1, keepdims=True)
    return triangulate_convex(np.concatenate((vertices, new_vertices)))


def icosphere(radius=1, position=(0, 0, 0), refinement=1):
    """
    generate an icosphere. Naively generating a sphere with sin and cos will cause the vertices to have an
    inhomogeneous distribution over the sphere surface. An icosphere avoids this.

    Parameters
    ----------
    radius: (optional) float
        the radius of the sphere being generated. The default value is 1
    position: (optional) tuple, (float, float)
        the position of the sphere. The default value is (0, 0, 0)
    refinement: int
        higher refinement will cause the mesh to have more vertices

    Returns
    -------
    icosphere_mesh: mesh
    """
    ico = icosahedron()

    for _ in range(refinement):
        ico = refine_sphere(ico)

    return Mesh((ico.vertices * radius) + position, ico.faces)
