# Workflow
## 1. Create or load a mesh of the simulation domain
Use the mesh constructors of `skfem.mesh` for simple shapes such as rectangles, circles, or L-shapes. 
Use *gmsh* or other software for complex shapes.

https://scikit-fem.readthedocs.io/en/latest/api.html#module-skfem.mesh

http://gmsh.info/

```python
import skfem

# create a rectangular mesh with skfem:

# load a mesh from file:
mesh = skfem.Mesh.load('mymesh.msh')
```

## 2. Choose a finite element from `skfem.element`. 
https://scikit-fem.readthedocs.io/en/latest/api.html#module-skfem.element

``` python
element = skfem.ElementTriP2()
fem = Helmholtz(mesh, element)
```

## 3. Assemble the domains and the boundaries
The Helmholtz equation is implemented as a general second-order partial differential equation:

$- \alpha (\frac{\partial^2 \Phi}{\partial x^2} + \frac{\partial^2 \Phi}{\partial y^2}) + \beta \Phi = f$.

The electromagnetic wave propagation can be formulated in two different ways, depending on the polarization of interest:
1) TM polarization: $\Phi = E_Z$, $\alpha = 1 / \mu_r$, $\beta = -k_0^2 \epsilon_r$, and $f = -j k_0 Z_0 J_Z$
2) TE polarization: $\Phi = H_Z$, $\alpha = 1 / \epsilon_r$, $\beta = -k_0^2 \mu_r$, and $f = \frac{1}{\epsilon_r} (\frac{\partial J_Y}{\partial x} - \frac{\partial J_X}{\partial y})$

with the normalized free-space vacuum wave number $k_0 = \frac{2 \pi}{\lambda_0} = \frac{2 \pi f}{c_0}$. 
Here, $c_0$ is propagation velocity (speed of light) in vacuum, expressed in *mesh units per second*.

For $f = 10 GHz$ and a mesh unit of $a = 1 mm$, one gets $k_0 \approx \frac{2 \pi 10 GHz}{3e8 m/s 1 / 1mm} \approx 0.209$.

For a mesh with two different subdomains labeled 'air' and 'water', the domains are assembled with their respective parameters:
```python
k0 = 0.2
eps_air = 1
mu_air = 1
eps_water = 81 - 10j
mu_water = 1
fem.assemble_subdomains(alpha={'air': 1 / mu_air, 
                               'water': 1 / mu_water}, 
                        beta={'air': -1 * k0 ** 2 * eps_air, 
                              'water': -1 * k0 ** 2 * eps_water}, 
                        f={'air': 0, 
                           'water': 0})
```

Similarly, the boundary conditions are defined. In this example, the upper and lower boundaries, labeled as 'bound_ymax' and 'bound_ymin', are supposed to be perfectly conducting (metallic) waveguide walls. 
The corresponding essential (Dirichlet) boundary condition is $E_Z = 0$:
```python
fem.assemble_boundaries_dirichlet(value={'bound_ymin': 0, 
                                         'bound_ymax': 0})
```

The boundaries on the left and right, labeled 'bound_xmin' and 'bound_xmax', are supposed to be waveguide interfaces that are artificially extended to behave like infinite waveguides for the fundamental propagation mode.
This can be formulated with the third-order boundary condition:
$\alpha (\frac{partial \Phi}{\partial x} \hat{x} + \frac{partial \Phi}{\partial y} \hat{y}) \hat{n} + \gamma \Phi = q$.
Here, $\hat{x}, $\hat{y}$, and $\hat{n}$ are unit vectors in x and y direction, and normal to the boundary, respectively.

```python
fem.assemble_boundaries_3rd(gamma={'bound_xmin': 1 / mu_air * 1j * k0, 
                                   'bound_xmax': 1 / mu_water * 1j * k0}, 
                            q={'bound_xmin': 1 / mu_air * 2j * k0, 
                               'bound_xmax': 0})
```

## 4. Solve the linear system for the field solution
```python
fem.solve()
```

## 5. Process the solution
After solving, the field solution is stored in `fem.phi` as a complex vector. The individual real and imaginary parts are also stored seperately in `fem.phi_re` and `fem.phi_im`.

The corresponding locations of the `N` elements in the solution vector on the mesh are stored in `fem.basis.doflocs`, which has the shape `(2, N)`.

scikit-fem offers helper functions to find certain elements, for example by means of labeled subdomains or boundaries:
```python
x_bound_xmin, y_bound_xmin = fem.basis.doflocs[:, fem.basis.get_dofs('bound_xmin')]
```

Plotting is simple with the functions in `skfem.visuals`:
```python
from skfem.visuals.matplotlib import plot, show
plot(fem.basis, fem.phi_re, colorbar=True)
show()
```