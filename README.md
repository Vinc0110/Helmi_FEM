# Workflow
## 1. Create or load a mesh of the simulation domain
Use the mesh constructors of `skfem.mesh` for simple shapes such as rectangles, circles, or L-shapes. 
Use *gmsh* or other software for complex shapes.

https://scikit-fem.readthedocs.io/en/latest/api.html#module-skfem.mesh

http://gmsh.info/

```python
import skfem
mesh = skfem.Mesh.load('mymesh.msh')
```

## 2. Choose a finite element from `skfem.element`. 
https://scikit-fem.readthedocs.io/en/latest/api.html#module-skfem.element

``` python
element = skfem.ElementTriP2()
```

## 3. Assemble the domains and the boundaries
The Helmholtz equation is implemented as a general second-order partial differential equation:

$- \alpha \frac{\partial^2 \Phi}{\partial x^2} - \alpha \frac{\partial^2 \Phi}{\partial y^2}$.

Formulation for TM polarization (solve for $E_\mathrm{Z}$ field): $\alpha = 1 / \mu$ and $\beta = $

```python
fem.assemble_subdomains(alpha={'air': 1 / mu_r}, beta={'air': -1 * k0 ** 2 * eps_r})
```

## 4. Solve the linear system for the field solution
```python
fem.solve()
```
