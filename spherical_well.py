import numpy as np
from spherical_bessel import j, beta

def A(N, l, a):
    '''
    Returns the normalization constant for the radial eigenfunction of an infinite spherical well.

    Parameters
    ----------
    N : int
        Radial quantum number.
    l : int
        Angular momentum quantum number.
    a : float
        The radius of the infinite spherical well.

    Returns
    -------
    float
        The normalization constant for the infinite spherical well eigenfunction.
    '''
    return np.sqrt(2) / ( a ** 1.5 * abs(j(l + 1, beta(N, l))))

def R(r, N, l, a):
    '''
    Returns the radial eigenfunction of an infinite spherical well.

    Parameters
    ----------
    r : float or numpy.ndarray
        Radial distance from the center of the well.
    N : int
        Radial quantum number.
    l : int
        Angular momentum quantum number.
    a : float
        The radius of the well.

    Returns
    -------
    float or numpy.ndarray
        The radial eigenfunction of the well.
    '''
    return A(N, l, a) * j(l, beta(N, l) * r / a)

def E(N, l, a, hbar = 1.0, m = 1.0):
    '''
    Returns the eigenenergy of an infinite spherical well.

    Parameters
    ----------
    N : int
        Radial quantum number.
    l : int
        Angular momentum quantum number.
    a : float
        The radius of the well.
    hbar : float, optinal
        Planck's reduced constant.
    m : float, optional
        Mass of the particle.

    Returns
    -------
    float
        The eigenenergy of the particle inside the well.
    '''
    return hbar ** 2 * beta(N, l) ** 2 / (2 * m * a ** 2)