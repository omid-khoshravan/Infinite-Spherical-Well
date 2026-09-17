import mpmath as mp
import numpy as np
from scipy.special import spherical_jn

def j(l, x):
    '''Spherical Bessel function of the first kind.'''
    return spherical_jn(l, x)

def beta(N, l):
    '''Nth zero of the spherical Bessel function.'''
    spherical_jn_zeros = mp.besseljzero(l + 0.5, N)
    return float(spherical_jn_zeros)