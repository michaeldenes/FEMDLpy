import numpy as np

def moving_average(a, n=2) :
    """
    Compute an n window moving average.
    """
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n
