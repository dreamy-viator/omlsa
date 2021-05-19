import numpy as np
from scipy.signal import hanning

def lnshift(x, t):
    szX = x.shape
    if szX[0] > 1:
        n = szX[0]
        y = np.hstack((x[t:n], x[:t]))
    else:
        n = szX[1]
        # TODO
    return y


def mat_hanning(n):
    y = hanning(n+1, False)
    y = y[1:]
    return y