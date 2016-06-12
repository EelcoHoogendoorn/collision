
import numpy as np

def dot(A, B):
    """dot arrays of vecs; contract over last indices"""
    return np.einsum('...i,...i->...', A, B)

def normalize(listofvec):
    """normalize an array of vecs"""
    norm = np.sqrt(dot(listofvec, listofvec))
    norm[norm==0] = 1
    return listofvec / np.expand_dims( norm, -1)
