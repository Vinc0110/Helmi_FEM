import numpy as np
import skfem
from helmi import Helmholtz
from helmi.helper import get_curvature
from skfem.visuals.matplotlib import draw, show, plot
import matplotlib.pyplot as mplt
from timeit import default_timer as timer
import slit_mesh

# wavelength (in layout units) and normalized wave number
wavelength = 2
k0 = 2 * np.pi / wavelength

# slit parameters
n_slits = 3
w_slit = 10
pitch = 20

# angle (direction of travel) of incident wave; in radiant [0, 2*pi]
theta_laser = 0

# material parameters
eps_air = 1.0
mu_air = 1.0

print('Loading mesh')
t1 = timer()
slit_mesh.mesh(n_slits, w_slit, pitch)
mesh = skfem.Mesh.load('./slit.msh')
t2 = timer()
print(f'Loading took {t2 - t1:.3f} s\n')

print('Init FEM')
element = skfem.ElementTriP2()
fem = Helmholtz(mesh, element)

print('\nAssembly')
t1 = timer()

# calculate incident field at left-side boundaries
x, y = fem.basis.doflocs[:, fem.basis.get_dofs('bounds_left')]
r = np.sqrt(x ** 2 + y ** 2)
kr = k0 * (x * np.cos(theta_laser) + y * np.sin(theta_laser))
phi0 = 1.0 * np.exp(-1j * kr)
dphi0_dn = -1j * kr / r * phi0

# simulation parameters for TM polarization
alpha_air = 1 / mu_air
beta_air = -1 * k0 ** 2 * eps_air
fem.assemble_subdomains(alpha={'air': alpha_air}, beta={'air': beta_air})
fem.assemble_boundaries_dirichlet(value={'bounds_aperture': 0})

# first order absorber
g_left = 1j * k0
x_right, y_right = fem.basis.doflocs[:, fem.basis.get_dofs('bounds_right')]
g_right = 1j * k0 + get_curvature(x_right, y_right) / 2

fem.assemble_boundaries_3rd(gamma={'bounds_left': alpha_air * g_left,
                                   'bounds_right': alpha_air * g_right},
                            q={'bounds_left': alpha_air * (dphi0_dn + g_left * phi0),
                               'bounds_right': 0})
t2 = timer()
print(f'Assembly took {t2 - t1:.3f} s\n')

print(f'Solving for {2 * fem.basis.N} unknowns')
t1 = timer()
fem.solve(direct=True, cuda=False)
t2 = timer()
print(f'Solving took {t2 - t1:.3f} s\n')

print('near2far()')
t1 = timer()
x_farfield = 1000
y_farfield = np.linspace(-1000, 1000, 401)
phi_farfield = np.zeros_like(y_farfield, dtype=complex)
for i in range(len(y_farfield)):
    phi_farfield[i] = fem.near2far(r=(x_farfield, y_farfield[i]), k=k0, field=fem.phi, boundaries=['bounds_right'])
t2 = timer()
print(f'near2far() took {t2 - t1:.3f} s\n')

# theoretical intensity distribution in diffraction gratings
# https://sites.ualberta.ca/~pogosyan/teaching/PHYS_130/FALL_2010/lectures/lect36/lecture36.html
sin_theta = y_farfield / x_farfield
slits = np.array(range(n_slits))
delta_theta = (slits[:, None] - 0.5 * (n_slits - 1)) * k0 * pitch * sin_theta
r = np.sqrt(x_farfield ** 2 + y_farfield ** 2)
phi_farfield_ideal = np.amax(np.abs(phi_farfield)) * np.sinc(k0 / 2 / np.pi * w_slit * sin_theta) * np.mean(np.exp(-1j * (k0 * r + delta_theta)), axis=0)

fig, ax = mplt.subplots(1, 2, figsize=(13, 7), gridspec_kw={'width_ratios': [3, 1]})
plot(fem.basis, fem.phi.real, shading='gouraud', colorbar=True, ax=ax[0])
ax[0].set_aspect(1)
ax[1].plot(np.abs(phi_farfield), y_farfield, label='FEM')
ax[1].plot(np.abs(phi_farfield_ideal), y_farfield, label='Theory')
ax[1].invert_xaxis()
ax[1].yaxis.set_label_position('right')
ax[1].yaxis.tick_right()
ax[1].legend()
mplt.tight_layout()
#mplt.show()
mplt.savefig(f'./slit_{n_slits}_{w_slit}_{pitch}.png', dpi=600)
mplt.close()
