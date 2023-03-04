import numpy as np


def get_curvature(x: np.array, y: np.array, wrapped=True) -> np.array:
    # magnitude and phase
    r = np.sqrt(x ** 2 + y ** 2)
    phi = np.arctan2(y, x)

    # sorts und unsorts phi from -pi to pi (flipped for positive rotation)
    i_sort = np.flip(np.argsort(phi))
    i_unsort = np.argsort(np.flip(i_sort))
    r_sorted = r[i_sort]
    phi_sorted = np.unwrap(phi[i_sort])

    if wrapped:
        # extend edges by wrapping it around
        r_ext = np.array([r_sorted[-2], r_sorted[-1], *r_sorted, r_sorted[0], r_sorted[1]])
        phi_ext = np.array([-1 * phi_sorted[-2], -1 * phi_sorted[-1], *phi_sorted, -1 * phi_sorted[0], -1 * phi_sorted[1]])
    else:
        r_ext = r_sorted
        phi_ext = phi_sorted

    # first and second derivative of r w.r.t. phi
    dr_dphi = np.gradient(r_ext, phi_ext, edge_order=2)
    d2r_dphi2 = np.gradient(dr_dphi, phi_ext, edge_order=1)

    if wrapped:
        dr_dphi = dr_dphi[2:-2]
        d2r_dphi2 = d2r_dphi2[2:-2]

    # Bronstein p. 255, eq. (3.502)
    k_sorted = (r_sorted ** 2 + 2 * dr_dphi ** 2 - r_sorted * d2r_dphi2) / ((r_sorted ** 2 + dr_dphi ** 2) ** 1.5)
    return k_sorted[i_unsort]
