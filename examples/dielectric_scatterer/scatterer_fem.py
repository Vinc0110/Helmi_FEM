import numpy as np
import skfem
from helmi import Helmholtz
from helmi.helper import get_curvature
from skfem.visuals.matplotlib import draw, show, plot
import matplotlib.pyplot as mplt
from timeit import default_timer as timer

# wavelength (in layout units) and normalized wave number
lambda_laser = 20
k0 = 2 * np.pi / lambda_laser

# angle (direction of travel) of incident wave; in radiant [0, 2*pi]
theta_laser = 0

# material parameters
eps_air = 1.0
mu_air = 1.0
eps_core = 1.4475 ** 2
#eps_core = 1
mu_core = 1.0
eps_cladding = 1.444 ** 2
#eps_cladding = 1
mu_cladding = 1.0

print('Loading mesh')
t1 = timer()
mesh = skfem.Mesh.load('./scatterer.msh')
t2 = timer()
print(f'Loading took {t2 - t1:.3f} s\n')

print('Init FEM')
element = skfem.ElementTriP3()
fem = Helmholtz(mesh, element)

print('\nAssembly')
t1 = timer()

# simulation parameters for TM polarization
alpha_air = 1 / mu_air
alpha_core = 1 / mu_core
alpha_cladding = 1 / mu_cladding
beta_air = -1 * k0 ** 2 * eps_air
beta_core = -1 * k0 ** 2 * eps_core
beta_cladding = -1 * k0 ** 2 * eps_cladding

# assemble all three subdomains
fem.assemble_subdomains(alpha={'air': alpha_air,
                               'core': alpha_core,
                               'cladding': alpha_cladding},
                        beta={'air': beta_air,
                              'core': beta_core,
                              'cladding': beta_cladding})

# assemble boundary conditions

# get boundary locations (coordinates) and curvature
x, y = fem.basis.doflocs[:, fem.basis.get_dofs()]
r = np.sqrt(x ** 2 + y ** 2)
phi = np.arctan2(y, x)
kappa = get_curvature(x, y, wrapped=True)

# calculate incident field at boundary
kr = k0 * (x * np.cos(theta_laser) + y * np.sin(theta_laser))
phi0 = 1.0 * np.exp(-1j * kr)
dphi0_dn = -1j * kr / r * phi0
d2phi0_ds2 = -1j * k0 / r * (y * np.cos(theta_laser) + x * np.sin(theta_laser) ** 2) * phi0

# first order absorber
g1 = 1j * k0 + kappa / 2
g2 = 0

# second order absorber
#g1 = 1j * k0 + kappa / 2 - 1j * kappa ** 2 / (8 * (1j * kappa - k0))
#g2 = -1j / (2 * (1j * kappa - k0))

q = alpha_air * (dphi0_dn + g1 * phi0 + g2 * d2phi0_ds2)
fem.assemble_boundaries_3rd(gamma={'bound': alpha_air * g1},
                            gamma2={'bound': alpha_air * g2},
                            q={'bound': q})
t2 = timer()
print(f'Assembly took {t2 - t1:.3f} s\n')

print(f'Solving for {2 * fem.basis.N} unknowns')
t1 = timer()
fem.solve(direct=True, cuda=False)
t2 = timer()
print(f'Solving took {t2 - t1:.3f} s\n')

#plot(fem.basis, fem.phi.real, shading='gouraud', colorbar=True)
#show()

print('near2far()')
t1 = timer()
x_farfield = 10000
y_farfield = np.linspace(-10000, 10000, 1001)
phi_farfield = np.zeros_like(y_farfield, dtype=complex)
for i in range(len(y_farfield)):
    phi_farfield[i] = fem.near2far(r=(x_farfield, y_farfield[i]), k=k0, field=fem.phi, boundaries=['bound'])
t2 = timer()
print(f'near2far() took {t2 - t1:.3f} s\n')

fig, ax = mplt.subplots(1, 2, figsize=(13, 7), gridspec_kw={'width_ratios': [3, 1]})
plot(fem.basis, fem.phi.real, shading='gouraud', colorbar=True, ax=ax[0])
draw(mesh, boundaries_only=True, ax=ax[0])
ax[0].set_aspect(1)
ax[0].set_title('Real part of transverse field')
ax[1].plot(np.real(phi_farfield), y_farfield, label='Real')
ax[1].plot(np.imag(phi_farfield), y_farfield, label='Imaginary')
ax[1].plot(np.abs(phi_farfield), y_farfield, label='Magnitude')
ax[1].legend()
ax[1].invert_xaxis()
ax[1].yaxis.set_label_position('right')
ax[1].yaxis.tick_right()
ax[1].set_ylabel('y farfield')
ax[1].set_title(f'Farfield at $x={x_farfield}$')
mplt.tight_layout()
#mplt.show()
mplt.savefig('./scatterer_fiber.png', dpi=300)
mplt.close()
