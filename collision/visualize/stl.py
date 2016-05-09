"""
stl handling module
"""
import numpy as np

import numpy_indexed as npi


def save_STL(filename, mesh):
    """save a mesh to plain stl. vertex ordering is assumed to be correct"""
    header      = np.zeros(80, '<c')
    triangles   = np.array(len(mesh.faces), '<u4')
    dtype       = [('normal', '<f4', 3,),('vertex', '<f4', (3,3)), ('abc', '<u2', 1,)]
    data        = np.empty(triangles, dtype)

    data['abc']    = 0     #standard stl cruft
    data['vertex'] = mesh.vertices[mesh.faces]
    data['normal'] = util.normalize(mesh.face_normals())

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
    return vertices, triangles.reshape(-1, 3)
