import numpy as np


def get_curvature(x: np.array, y: np.array) -> np.array:
    r = np.sqrt(x ** 2 + y ** 2)
    phi = np.arctan2(y, x)
    # sorts phi from -pi to pi (negative rotation; consider during differentiation)
    i_sort = np.argsort(phi)
    r = r[i_sort]
    phi = np.unwrap(phi[i_sort])

    # first and second derivative of r w.r.t. phi
    dr_dphi = np.gradient(r, phi)
    d2r_dphi2 = np.gradient(dr_dphi, phi)

    # Bronstein p. 255, eq. (3.502)
    return (r ** 2 + 2 * dr_dphi ** 2 - r * d2r_dphi2) / ((r ** 2 + dr_dphi ** 2) ** 1.5)
