# Infinite-Spherical-Well

A numerical solution and visualization of a particle confined to an infinite spherical well that suddenly doubles in radius.

## Overview

This project is to numerically solve for the radial wavefunction of a particle that was initially at an eigenstate of an infinite spherical well, that suddenly doubles in radius from $a$ to $2a$ at time $t = 0$.

The project was inspired by a QM class assignment. The original problem was not a coding project; we were asked to analytically or numerically derive the time-dependent wavefunction after the expansion, for the case of a particle that started out in the ground state. This repository develops a numerical implementation of that problem and allows the initial radial quantum number and angular momentum quantum number to be specified arbitrarily.

The calculation is performed using spherical Bessel functions and an expansion of the initial state in the eigenstates of the expanded well.

## Physical Setup

For an infinite spherical well of radius $a$, the potential is 

$$
V(r) = 
\begin{cases}
0, & r < a, \\
\infinity, & r\geq a.
\end{cases}
$$

The stationary states can be written as 
$$\psi_{Nlm}(r, \theta, \phi) = R_{Nl}(r)Y_l^m(\theta, \phi)$$
where $N$ is the radial quantum number, $l$ is the angular momentum quantum number, and $m$ is the magnetic quantum number.

The radial part of the normalized eigenfunction is
$$R_{Nl}^{(a)}(r) = A_{Nl}^{(a)}j_l(\frac{\beta_{Nl}r}{a})$$
where $j_l$ is a spherical Bessel function and $\beta_{Nl}$ is its N-th positive zero.

The normalization constant is
$$A_{Nl}^{(a)} = \frac{\sqrt{2}}{a^{3/2}|j_{l+1}(\beta_{Nl})|}$$

The corresponding energy eigenvalue is
$$E_{Nl}^{(a)} = \frac{\hbar^2\beta_{Nl}^2}{a^{3/2}|j_{l + 1}(\beta_{Nl})|$$

## Sudden Expansion

Initially, the particle is in an eigenstate of the well of radius $a$,
$$\psi(r, 0^-) = \psi_{Nlm}^{(a)}(r, \theta, \phi)$$

When the radius suddenly changes to $2a$, the wavefunction does not have time to change instantaneously. Therefor,
$$\psi(r, 0^+) = \psi(r, 0^-)$$
However, this state is no longer an eigenstate of the expanded well. It must instead be expressed as a superposition of the eigenstates of the new well:
$$\psi(\mathbf{r}, t) = \sum_{N' = 1}^{\infinity} \alpha_{N'}\psi_{N'lm}^{(2a)}(\mathbf{r})e^{-iE_{N'l}^{(2a)}t/\hbar}$$
Because the expansion preserves spherical symmetry, $l$ and $m$ remain unchanged. Only the radial quantum number is summed over.

The expansion coefficients are given by
$$\alpha_{N'} = \langleN'lm;2a|Nlm;a\rangle$$
After carrying out the angular integration, this inner product reduces to:
$$\alpha_{N'} = \int_0^a r^2R_{Nl}^{(2a)}(r)dr$$
The coefficients are evaluated numerically in the code, with a direct numerical integration used in the special case where the closed-form expression becomes numerically singular.

## Time-Dependent Radial Wavefunction

The radial part of the time-dependent wavefunction is therefor
$$R(r, t) = \sum_{N' = 1}^{N_max} alpha_{N'}R_{N'l}^{(2a)}(r)e^{-iE_{N'l}^{(2a)}t/\hbar},$$
where $N_{\max}$ determines the number of expanded-well eigenstates included in the numerical approximation.

The quantity visualized in the animation is the radial probability density,
$$P(r, t) = r^2 |R(r, t)|^2$$
It satisfies
$$\int_0^{2a}P(r, t)dr = 1$$
Thus, the animation shows how the radial probability distribution evolves after the sudden expansion of the well.

## Numerical Implementation

The calculation is divided into several Python modules:
- `spherical_bessel.py` Provides spherical Bessel functions and calculates their positive zeros.
- `spherical_well.py` Defines the normalized radial eigenfunctions and eigenenergies of the spherical well.
- `sudden_expansion.py` Calculates the expansion coefficients and constructs the time-dependent radial wavefunction and radial probability density.
- `expanding_well.ipynb` Sets the physical parameters, evaluates the radial probability density, and generates the time-evolution animation.

## Visualization

The radial probability density is animated over a chosen time interval following the sudden expansion. The resulting animation is saved as a GIF:
![Radial probability density animation](infinite_spherical_well.gif)

## Requirements
The project uses:
- Python
- NumPy
- SciPy
- mpmath
- Matplotlib
- Jupyter Notebook

## Notes
The default implementation uses natural units with
$$
\hbar = 1,
m = 1
$$
unless different values are supplied.

The number of terms in the expansion is finite in the numerical calculation. Increasing $N_{\max}$ improves the representation of the initial state in the expanded-well basis, at the cost of longer computation time.
