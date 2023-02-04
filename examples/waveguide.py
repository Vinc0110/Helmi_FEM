from helmi import Helmholtz
import skfem
import numpy as np

x_pts = np.linspace(0, 100, 101)
y_pts = np.linspace(-5, 5, 21)
mesh = skfem.MeshTri.init_tensor(x_pts, y_pts)
mesh = mesh.with_subdomains({'air': lambda x: x[0] < 50,
                             'plastic': lambda x: x[0] >= 50})
mesh = mesh.with_boundaries({'bound_xmin': lambda x: np.isclose(x[0], x_pts[0]),
                             'bound_xmax': lambda x: np.isclose(x[0], x_pts[-1]),
                             'bound_ymin': lambda x: np.isclose(x[1], y_pts[0]),
                             'bound_ymax': lambda x: np.isclose(x[1], y_pts[-1])})

element = skfem.ElementTriP2()
fem = Helmholtz(mesh, element)

k0 = 0.5
eps_air = 1
mu_air = 1
eps_plastic = 2 - 0.1j
mu_plastic = 1
fem.assemble_subdomains(alpha={'air': 1 / mu_air,
                               'plastic': 1 / mu_plastic},
                        beta={'air': -1 * k0 ** 2 * eps_air,
                              'plastic': -1 * k0 ** 2 * eps_plastic},
                        f={'air': 0,
                           'plastic': 0})

fem.assemble_boundaries_dirichlet(value={'bound_ymin': 0,
                                         'bound_ymax': 0})

fem.assemble_boundaries_3rd(gamma={'bound_xmin': 1 / mu_air * 1j * k0,
                                   'bound_xmax': 1 / mu_plastic * 1j * k0},
                            q={'bound_xmin': 1 / mu_air * 2j * k0,
                               'bound_xmax': 0})

fem.solve()

from skfem.visuals.matplotlib import plot
import matplotlib.pyplot as mplt

fig, ax = mplt.subplots(2, 1)
plot(fem.basis, fem.phi_re, ax=ax[0])
plot(fem.basis, fem.phi_im, ax=ax[1])
ax[0].set_aspect(1)
ax[1].set_aspect(1)
mplt.tight_layout()
mplt.show()
