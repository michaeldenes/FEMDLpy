from numpy.linalg import norm

def normed(x):
    """
    Scale by inf norm of an array.
    """
    y = x/norm(x,inf)
    return y
