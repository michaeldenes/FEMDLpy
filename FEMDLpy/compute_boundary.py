import numpy as np

def compute_boundary(triangles):
    """
    Compute the boundary of a triangulation, where the boundary is defined
    as any node of an edge that appears only once (that is, not shared
    between any two triangle faces).
    """
    edges = np.concatenate((triangles[:,[0,1]],triangles[:,[0,2]],triangles[:,[1,2]]))
    edges = np.sort(edges, axis=1)

    uniqueedges, uniquecounts = np.unique(edges,axis=0,return_counts=True)
    uniqueedges = uniqueedges[np.where(uniquecounts == 1)]

    boundary_indices = np.unique(uniqueedges.flatten())
    return boundary_indices
