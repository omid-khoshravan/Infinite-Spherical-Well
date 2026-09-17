import numpy as np
from scipy.integrate import quad
from spherical_bessel import j, beta
from spherical_well import A, R, E



def alpha_N_prime(N_prime, N, l, a):
    '''
    Returns the expansion coefficient alpha_N' after the doubling of the well radius from a to 2a.

    Parametes
    ---------
    N_prime : int
        Radial quantum number of the eigenstate in the expanded well.
    N : int
        Initial radial quantum number.
    l : int
        Angular momentum quantum number.
    a : float
        Initial radius of the well.

    Returns
    -------
    float
        Expansion coefficient alpha_N'.
    '''

    Normalization = A(N, l, a) * A(N_prime, l, 2 * a)

    b = beta(N, l)
    b_prime = beta(N_prime, l)

    denominator = 4 * b ** 2 - b_prime ** 2

    if np.isclose(denominator, 0):
        integral, error = quad(lambda r: r ** 2 * j(l, b * r / a) * j(l, b_prime * r / (2 * a)), 0, a)
        return float(Normalization * integral)
    else:
        numerator = 4 * a ** 3 * b * j(l, b_prime / 2) * j(l + 1, b)
        return float(Normalization * numerator / denominator)


    
def R_t(r, t, N, l, a, N_max, hbar = 1, m =1):
    '''
    Returns the time-dependent radial wavefunction after the well's radius doubles.

    Parameters
    ----------
    r : float
        Radial distance from the center of the well.
    t : float
        Time.
    N : int
        Initial radial quantum number.
    l : int
        Angular momentum quantum number.
    a : float
        Initial radius of the well.
    N_max : int
        Maximum radial quantum number included in the expansion.
    hbar : float, optional
        Planck's reduced constant.
    m : float, optional
        Mass of the particle.

    Return
    ------
    complex
        Time-dependant radial wavefunction R(r, t).
    '''
    sum = np.zeros_like(r, dtype = complex)

    for N_prime in range(1, N_max + 1):
        alpha = alpha_N_prime(N_prime, N, l, a)

        sum += alpha * R(r, N_prime, l, 2 * a) * np.exp(-1j * E(N_prime, l, 2 * a, hbar, m) * t / hbar)

    return sum



def P_r(r, t, N, l, a, N_max):
    '''
    Returns the radial probability density after the doubling of the radius.

    Prameters
    ---------
    r : float
        Radial distance from the center of the well.
    t : float
        Time.
    N : int
        Initial radial quantum number.
    l : int
        Angular momentum quantum number.
    a : float
        Initial radius of the well.
    N_max : int
        Maximum radial quantum number included in the expansion.

    Returns
    -------
    float
        Radial probability density r^2 |R(r, t)|^2
    '''

    return r ** 2 * abs(R_t(r, t, N, l, a, N_max)) ** 2