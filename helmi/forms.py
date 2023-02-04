import skfem
from skfem.helpers import dot
import numpy as np


@skfem.BilinearForm
def helmholtz_laplace(u, v, w):
    # Laplace part in LHS of scalar Helmholtz equation for z-component of a field `phi` (Jin, p. 77, eq. 4.1):
    # (alpha_x * d²/dx² + alpha_y d²/dy² - beta) * phi_z = -f
    # simplified with alpha_x = alpha_y = alpha to
    # [alpha * (d²/dx² + d²/dy²) - beta] * phi_z = -f
    #
    # --> F = integrate [alpha * (d²/dx² + d²/dy²) * phi_z] dx dy
    #
    # elemental bilinear coefficients (Jin, p. 85, eq. 4.30):
    # K_ij(u, v) = alpha_x * du/dx * dv/dx + alpha_y * du/dy * dv/dy
    return w['alpha'] * dot(u.grad, v.grad)

@skfem.BilinearForm
def helmholtz_mass(u, v, w):
    # Mass part in LHS of scalar Helmholtz equation for z-component of a field `phi` (Jin, p. 77, eq. 4.1):
    # (alpha_x * d²/dx² + alpha_y d²/dy² - beta) * phi_z = -f
    # simplified with alpha_x = alpha_y = alpha to
    # [alpha * (d²/dx² + d²/dy²) - beta] * phi_z = -f
    #
    # --> F = integrate [beta * phi_z] dx dy
    #
    # elemental bilinear coefficients (Jin, p. 85, eq. 4.30):
    # K_ij(u, v) = beta * u * v
    return w['beta'] * u * v


@skfem.LinearForm
def helmholtz_excitation(u, w):
    # RHS of scalar Helmholtz equation for z-component of a field `phi` (Jin, p. 85, eq. 4.31):
    # (alpha_x * d²/dx² + alpha_y d²/dy² - beta) * phi_z = -f
    # simplified with alpha_x = alpha_y = alpha to
    # [alpha * (d²/dx² + d²/dy²) - beta] * phi_z = -f
    #
    # F = integrate [f * phi_z] dx dy
    # elemental linear coefficients (Jin, p. 79, eq. 4.8):
    # b_i(u) = f * u
    return w['f'] * u


@skfem.Functional
def helmholtz_near2far(w):
    r = w['r']
    e_r = r / np.sqrt(r[0] ** 2 + r[1] ** 2)
    r_s = w['x']
    n_s = w['n']
    k = w['k']
    phi = w['phi']

    return (dot(e_r, n_s) * phi + 1j / k * dot(n_s, phi.grad)) * np.exp(1j * k * dot(e_r, r_s))

