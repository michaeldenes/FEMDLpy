import numpy as np

def softthresh(x,lam):
    """
    Compute the softthreshhold for SEBA.
    """
    return np.sign(x)*np.maximum((np.abs(x)-lam),0)
