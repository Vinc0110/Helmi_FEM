import numpy as np
import skfem
from helmi import Helmholtz
from skfem.visuals.matplotlib import plot
from scipy.constants import epsilon_0, mu_0
import matplotlib.pyplot as mplt
from timeit import default_timer as timer
import horn_mesh

# parameters
unit = 1e-3
plane = 'h'
f = 140e9
n = 1
eps_r = 1.0
mu_r = 1.0
k0 = 2 * np.pi * f * np.sqrt(epsilon_0 * mu_0) * unit
z0 = np.sqrt(mu_0 / epsilon_0) * unit

print('Meshing...')
t1 = timer()
horn_mesh.mesh(plane)
mesh = skfem.Mesh.load(f'./horn_{plane}-plane.msh')
t2 = timer()
print(f'Meshing took {t2 - t1:.3f} s\n')

print('Init FEM...')
element = skfem.ElementTriP2()
fem = Helmholtz(mesh, element)

print('\nAssembly...')
t1 = timer()
if plane == 'e':
    # TE; phi=Hz
    alpha = 1 / eps_r
    beta = -1 * k0 ** 2 * mu_r
    coeff_complementary = -1 * z0 / k0 * alpha
else:
    # TM; phi=Ez
    alpha = 1 / mu_r
    beta = -1 * k0 ** 2 * eps_r
    coeff_complementary = 1 / z0 / k0 * alpha
    fem.assemble_boundaries_dirichlet(value={'bound_ymin': 0, 'bound_ymax': 0})

fem.assemble_subdomains(alpha={'air': alpha}, beta={'air': beta})
fem.assemble_boundaries_3rd(gamma={'bound_feed': alpha * 1j * k0, 'bound_freespace': alpha * 1j * k0},
                            q={'bound_feed': alpha * 2j * k0, 'bound_freespace': 0})
t2 = timer()
print(f'Assembly took {t2 - t1:.3f} s')

print('\nSolving...')
t1 = timer()
fem.solve(direct=True)
t2 = timer()
print(f'Solving took {t2 - t1:.3f} s\n')

print('near2far()...')
t1 = timer()
r_farfield = 100
theta_farfield = np.linspace(-0.5 * np.pi, 0.5 * np.pi, 181)
phi_farfield = np.zeros_like(theta_farfield, dtype=complex)
for i, theta_i in enumerate(theta_farfield):
    r_x = r_farfield * np.cos(theta_i)
    r_y = r_farfield * np.sin(theta_i)
    phi_farfield[i] = fem.near2far(r=(r_x, r_y), k=k0, field=fem.phi, boundaries=['bound_freespace'])
t2 = timer()
print(f'near2far() took {t2 - t1:.3f} s\n')

print('Post processing...')
phi_mag_farfield = np.abs(phi_farfield)
phi_db_farfield = 20 * np.log10(phi_mag_farfield / np.amax(phi_mag_farfield))

# dof locations (x, y) at x=xmin and at x=xmax
x_feed, y_feed = fem.basis.doflocs[:, fem.basis.get_dofs('bound_feed')]

# indices of center locations (y=0) at x=xmin and at x=xmax
i_feed_y0 = np.argmin(np.abs(y_feed - 0))

# get fields at (x=xmin, y=0) and at (x=xmax, y=0):
phi_feed = fem.phi[fem.basis.get_dofs('bound_feed')][i_feed_y0]

# calculate reflection and transmission coefficients
r = phi_feed / 1.0 - 1
print(f'r = {np.abs(r) ** 2} = {20 * np.log10(np.abs(r)):.1f} dB')

# complementary fields (in x-y plane);
phi_comp_x_re = -1 * coeff_complementary * fem.basis.project(fem.basis.interpolate(fem.phi_im).grad[1])
#phi_comp_x_im = coeff_complementary * fem.basis.project(fem.basis.interpolate(fem.phi_re).grad[1])
#phi_comp_x = phi_comp_x_re + 1j * phi_comp_x_im

phi_comp_y_re = -1 * coeff_complementary * fem.basis.project(fem.basis.interpolate(fem.phi_im).grad[0])
#phi_comp_y_im = coeff_complementary * fem.basis.project(fem.basis.interpolate(fem.phi_re).grad[0])
#phi_comp_y = phi_comp_y_re + 1j * phi_comp_y_im

fig, ax = mplt.subplots(2, 1, figsize=(8, 6))
fig.suptitle(f'Real parts of fields at f = {f * 1e-9} GHz')
#draw(mesh, ax=ax[0])
plot(fem.basis, fem.phi_re, colorbar=True, ax=ax[0])
ax[0].set_aspect('equal')
ax[0].set_title('Transverse field (Ez or Hz)')
#draw(mesh, ax=ax[1])
plot(fem.basis, phi_comp_x_re + phi_comp_y_re, colorbar=True, ax=ax[1])
ax[1].set_aspect('equal')
ax[1].set_title('Complementary field (Hx+Hy or Ex+Ey)')
fig.tight_layout()
mplt.savefig(f'./horn_{plane}-plane_fields.png')
mplt.close()

mplt.figure()
mplt.polar(theta_farfield, phi_db_farfield)
mplt.title(f'Radiation pattern ({plane.upper()}-Plane, {f * 1e-9} GHz)')
mplt.tight_layout()
mplt.savefig(f'./horn_{plane}-plane_pattern_polar.png')
mplt.close()

mplt.figure()
mplt.plot(np.rad2deg(theta_farfield), phi_db_farfield)
mplt.title(f'Radiation pattern ({plane.upper()}-Plane, {f * 1e-9} GHz)')
mplt.tight_layout()
mplt.savefig(f'./horn_{plane}-plane_pattern_rect.png')
mplt.close()
