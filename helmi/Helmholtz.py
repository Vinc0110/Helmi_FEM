import skfem
from .forms import *
import numpy as np
from scipy.sparse import coo_matrix, bmat, identity, isspmatrix_coo, block_diag
import warnings

try:
    import cupy
    from cupyx.scipy.sparse import csr_matrix as cu_csr_matrix
    from cupyx.scipy.sparse.linalg import spsolve as cu_spsolve
except ImportError:
    pass


class Helmholtz:
    def __init__(self, mesh: skfem.Mesh, element: skfem.Element) -> None:
        self.mesh = mesh
        self.element = element
        self.basis = skfem.Basis(mesh, element)

        self.phi = None
        self.phi_re = None
        self.phi_im = None
        self._boundaries_dirichlet = None
        self._phi_dirichlet_re = None
        self._phi_dirichlet_im = None
        self._A_laplace_re = None
        self._A_laplace_im = None
        self._A_mass_re = None
        self._A_mass_im = None
        self._b_re = None
        self._b_im = None

    def assemble_subdomains(self, alpha: dict = None, beta: dict = None, f: dict = None) -> None:
        alpha_re = self.basis.zeros()
        alpha_im = self.basis.zeros()
        beta_re = self.basis.zeros()
        beta_im = self.basis.zeros()
        f_re = self.basis.zeros()
        f_im = self.basis.zeros()

        # populate parameters for all subdomains
        for subdomain in self.mesh.subdomains.keys():
            dofs = self.basis.get_dofs(elements=subdomain)
            if alpha is not None and subdomain in alpha:
                alpha_re[dofs] = alpha[subdomain].real
                alpha_im[dofs] = alpha[subdomain].imag
            if beta is not None and subdomain in beta:
                beta_re[dofs] = beta[subdomain].real
                beta_im[dofs] = beta[subdomain].imag
            if f is not None and subdomain in f:
                f_re[dofs] = f[subdomain].real
                f_im[dofs] = f[subdomain].imag

        a_laplace_re = skfem.asm(helmholtz_laplace, self.basis, alpha=alpha_re)
        if self._A_laplace_re is None:
            self._A_laplace_re = a_laplace_re
        else:
            self._A_laplace_re += a_laplace_re

        a_laplace_im = skfem.asm(helmholtz_laplace, self.basis, alpha=alpha_im)
        if self._A_laplace_im is None:
            self._A_laplace_im = a_laplace_im
        else:
            self._A_laplace_im += a_laplace_im

        a_mass_re = skfem.asm(helmholtz_mass, self.basis, beta=beta_re)
        if self._A_mass_re is None:
            self._A_mass_re = a_mass_re
        else:
            self._A_mass_re += a_mass_re

        a_mass_im = skfem.asm(helmholtz_mass, self.basis, beta=beta_im)
        if self._A_mass_im is None:
            self._A_mass_im = a_mass_im
        else:
            self._A_mass_im += a_mass_im

        b_re = skfem.asm(helmholtz_excitation, self.basis, f=f_re)
        if self._b_re is None:
            self._b_re = b_re
        else:
            self._b_re += b_re

        b_im = skfem.asm(helmholtz_excitation, self.basis, f=f_im)
        if self._b_im is None:
            self._b_im = b_im
        else:
            self._b_im += b_im

    def assemble_boundaries_3rd(self, gamma: dict = None, q: dict = None) -> None:
        gamma_re = self.basis.zeros()
        gamma_im = self.basis.zeros()
        q_re = self.basis.zeros()
        q_im = self.basis.zeros()

        # populate parameters for all boundaries
        boundaries_lhs = []
        boundaries_rhs = []
        for boundary in self.mesh.boundaries.keys():
            dofs = self.basis.get_dofs(boundary)
            if gamma is not None and boundary in gamma:
                gamma_re[dofs] = gamma[boundary].real
                gamma_im[dofs] = gamma[boundary].imag
                if boundary not in boundaries_lhs:
                    boundaries_lhs.append(boundary)
            if q is not None and boundary in q:
                q_re[dofs] = q[boundary].real
                q_im[dofs] = q[boundary].imag
                if boundary not in boundaries_rhs:
                    boundaries_rhs.append(boundary)

        c = skfem.asm(helmholtz_mass, self.basis.boundary(boundaries_lhs), beta=gamma_re)
        if self._A_mass_re is None:
            self._A_mass_re = c
        else:
            self._A_mass_re += c

        c = skfem.asm(helmholtz_mass, self.basis.boundary(boundaries_lhs), beta=gamma_im)
        if self._A_mass_im is None:
            self._A_mass_im = c
        else:
            self._A_mass_im += c

        c = skfem.asm(helmholtz_excitation, self.basis.boundary(boundaries_rhs), f=q_re)
        if self._b_re is None:
            self._b_re = c
        else:
            self._b_re += c

        c = skfem.asm(helmholtz_excitation, self.basis.boundary(boundaries_rhs), f=q_im)
        if self._b_im is None:
            self._b_im = c
        else:
            self._b_im += c

    def assemble_boundaries_eigenmode(self, k0: float, n_mode: dict) -> tuple:
        gamma_re = self.basis.zeros()
        gamma_im = self.basis.zeros()
        q_re = self.basis.zeros()
        q_im = self.basis.zeros()

        # populate parameters for all boundaries
        boundaries_lhs = []
        boundaries_rhs = []
        phi_n = {}
        gamma_n = {}
        for boundary in self.mesh.boundaries.keys():
            if boundary in n_mode:
                dofs = self.basis.get_dofs(boundary)
                kt2, phi, ebasis = self.solve_eigenmodes_boundary(boundary, n_modes=n_mode[boundary])
                phi_n[boundary] = phi[:, n_mode[boundary] - 1]
                k2_eigen = k0 ** 2 - kt2[n_mode[boundary] - 1]
                if k2_eigen < 0:
                    # evanescent mode
                    gamma = -1 * np.sqrt(np.abs(k2_eigen))
                else:
                    # propagating mode
                    gamma = 1j * np.sqrt(k2_eigen)

                gamma_n[boundary] = gamma
                gamma_re[dofs] = gamma.real
                gamma_im[dofs] = gamma.imag

                q = 2 * gamma * phi_n[boundary]
                q_re[dofs] = q.real
                q_im[dofs] = q.imag

                if boundary not in boundaries_lhs:
                    boundaries_lhs.append(boundary)

                if boundary not in boundaries_rhs:
                    boundaries_rhs.append(boundary)

        c = skfem.asm(helmholtz_mass, self.basis.boundary(boundaries_lhs), beta=gamma_re)
        if self._A_mass_re is None:
            self._A_mass_re = c
        else:
            self._A_mass_re += c

        c = skfem.asm(helmholtz_mass, self.basis.boundary(boundaries_lhs), beta=gamma_im)
        if self._A_mass_im is None:
            self._A_mass_im = c
        else:
            self._A_mass_im += c

        c = skfem.asm(helmholtz_excitation, self.basis.boundary(boundaries_rhs), f=q_re)
        if self._b_re is None:
            self._b_re = c
        else:
            self._b_re += c

        c = skfem.asm(helmholtz_excitation, self.basis.boundary(boundaries_rhs), f=q_im)
        if self._b_im is None:
            self._b_im = c
        else:
            self._b_im += c

        return gamma_n, phi_n

    def assemble_boundaries_dirichlet(self, value: dict):
        self._boundaries_dirichlet = []
        self._phi_dirichlet_re = self.basis.zeros()
        self._phi_dirichlet_im = self.basis.zeros()

        for boundary in self.mesh.boundaries.keys():
            if value is not None and boundary in value:
                self._boundaries_dirichlet.append(boundary)
                if value[boundary] != 0:
                    dofs = self.basis.get_dofs(boundary)
                    self._phi_dirichlet_re[dofs] = value[boundary].real
                    self._phi_dirichlet_im[dofs] = value[boundary].imag

    @staticmethod
    def _stagger_re_im(a_re, a_im):
        if a_re.shape != a_im.shape:
            return None
        # build combined A and b for linear system (A*x=b) with staggered real and imaginary parts:

        # block staggering
        # A = bmat([[A_re, -1 * A_im], [A_im, A_re]])
        # b = np.hstack([b_re, b_im])

        # elemental staggering
        if len(a_re.shape) == 2:
            # NxN matrix
            if not isspmatrix_coo(a_re):
                a_re = a_re.tocoo()
            if not isspmatrix_coo(a_im):
                a_im = a_im.tocoo()
            a = coo_matrix(([*a_re.data, *a_re.data, *(-1 * a_im.data), *a_im.data],
                            ([*(2 * a_re.row), *(2 * a_re.row + 1), *(2 * a_im.row), *(2 * a_im.row + 1)],
                             [*(2 * a_re.col), *(2 * a_re.col + 1), *(2 * a_im.col + 1), *(2 * a_im.col)])),
                           shape=(2 * a_re.shape[0], 2 * a_re.shape[1])).tocsr()
        elif len(a_re.shape) == 1:
            # vector
            a = np.empty(2 * a_re.shape[0])
            a[0::2] = a_re
            a[1::2] = a_im
        else:
            return None

        return a

    def solve(self, direct=True, cuda=False):
        A = self._stagger_re_im(a_re=self._A_laplace_re + self._A_mass_re,
                                a_im=self._A_laplace_im + self._A_mass_im)
        b = self._stagger_re_im(a_re=self._b_re, a_im=self._b_im)

        if self._boundaries_dirichlet is not None and len(self._boundaries_dirichlet) > 0:
            dofs = self.basis.get_dofs(self._boundaries_dirichlet).flatten()
            D = np.hstack((2 * dofs, 2 * dofs + 1))
            x = np.empty(2 * self.basis.N)
            x[0::2] = self._phi_dirichlet_re
            x[1::2] = self._phi_dirichlet_im
        else:
            D = None

        if direct:
            solver = skfem.solver_direct_scipy()
        else:
            # Jacobi preconditioning
            #M = skfem.build_pc_diag(A)
            #solver = skfem.solver_iter_krylov(M=M, x0=M.dot(b))

            # iLU preconditioning
            M = skfem.build_pc_ilu(A)
            solver = skfem.solver_iter_krylov(M=M, x0=M.matvec(b))

        if D is not None:
            if cuda:
                AII, bI, xI, I = skfem.condense(A, b, x=x, D=D)

                A_cu = cu_csr_matrix(AII)
                b_cu = cupy.array(bI)
                phi = x.copy()
                phi[I] = cupy.asnumpy(cu_spsolve(A_cu, b_cu))
            else:
                phi = skfem.solve(*skfem.condense(A, b, x=x, D=D), solver=solver)
        else:
            if cuda:
                A_cu = cu_csr_matrix(A)
                b_cu = cupy.array(b)
                phi_cu = cu_spsolve(A_cu, b_cu)
                phi = cupy.asnumpy(phi_cu)
            else:
                phi = skfem.solve(A, b, solver=solver)

        self.phi_re = phi[0::2]
        self.phi_im = phi[1::2]
        self.phi = self.phi_re + 1j * self.phi_im

    def solve_eigenmodes_boundary(self, boundary, n_modes=1):
        # dofs = self.basis.get_dofs(facets=boundary)
        # x, y = self.basis.doflocs[:, dofs]
        # mesh = skfem.MeshLine(p=y)

        def mapping(p):
            # define mapping from 2d coordinates of boundary elements (p = [[x1, x2, ...], [y1, y2, ...]])
            # to 1d coordinates (positions) of line elements
            return p[1]

        mesh, facets = self.mesh.trace(facets=boundary, mtype=skfem.MeshLine, project=mapping)

        if self.element.maxdeg == 1:
            element = skfem.ElementLineP1()
        elif self.element.maxdeg == 2:
            element = skfem.ElementLineP2()
        else:
            warnings.warn('Using higher order Legendre elements; results might be wrong due to incompatible basis '
                          'functions.', stacklevel=2)
            element = skfem.ElementLinePp(self.element.maxdeg)

        basis = skfem.Basis(mesh, element)

        A = skfem.asm(helmholtz_laplace, basis, alpha=basis.ones())
        B = skfem.asm(helmholtz_mass, basis, beta=basis.ones())

        solver = skfem.solver_eigen_scipy_sym(k=n_modes, sigma=1)
        kt2, phi = skfem.solve(*skfem.condense(A, B, D=basis.get_dofs()), solver=solver)

        return kt2, phi, basis

    def solve_eigenmodes(self, psi_x: float, psi_y: float, n_modes: int = 1):
        xmin, ymin = np.amin(self.basis.doflocs, axis=1)
        xmax, ymax = np.amax(self.basis.doflocs, axis=1)

        # boundary edge dofs including corners
        dofs_xmin = self.basis.get_dofs(lambda p: p[0] == xmin).flatten()
        dofs_xmax = self.basis.get_dofs(lambda p: p[0] == xmax).flatten()
        dofs_ymin = self.basis.get_dofs(lambda p: p[1] == ymin).flatten()
        dofs_ymax = self.basis.get_dofs(lambda p: p[1] == ymax).flatten()

        # boundary corner dofs
        dof_c00 = np.intersect1d(dofs_xmin, dofs_ymin)[0]
        dof_c10 = np.intersect1d(dofs_xmax, dofs_ymin)[0]
        dof_c01 = np.intersect1d(dofs_xmin, dofs_ymax)[0]
        dof_c11 = np.intersect1d(dofs_xmax, dofs_ymax)[0]

        # boundary edge dofs without corners
        dofs_xmin = dofs_xmin[np.logical_and(dofs_xmin != dof_c00, dofs_xmin != dof_c01)]
        dofs_xmax = dofs_xmax[np.logical_and(dofs_xmax != dof_c10, dofs_xmax != dof_c11)]
        dofs_ymin = dofs_ymin[np.logical_and(dofs_ymin != dof_c00, dofs_ymin != dof_c10)]
        dofs_ymax = dofs_ymax[np.logical_and(dofs_ymax != dof_c01, dofs_ymax != dof_c11)]

        # strictly internal dofs
        dofs_int = self.basis.complement_dofs(self.basis.get_dofs()).flatten()

        # vector sizes
        n_int = len(dofs_int)
        n_xmin = len(dofs_xmin)
        n_xmax = len(dofs_xmax)
        n_ymin = len(dofs_ymin)
        n_ymax = len(dofs_ymax)

        # # real parts
        # A[2 * dofs_xmin, 2 * dofs_xmin] += A[2 * dofs_xmax, 2 * dofs_xmax] * np.cos(psi_x)
        # B[2 * dofs_xmin, 2 * dofs_xmin] += B[2 * dofs_xmax, 2 * dofs_xmax] * np.cos(psi_x)
        # A[2 * dofs_ymin, 2 * dofs_ymin] += A[2 * dofs_ymax, 2 * dofs_ymax] * np.cos(psi_y)
        # B[2 * dofs_ymin, 2 * dofs_ymin] += B[2 * dofs_ymax, 2 * dofs_ymax] * np.cos(psi_y)
        # A[2 * dof_c00, 2 * dof_c00] += A[2 * dof_c01, 2 * dof_c01] * np.cos(psi_x)
        # B[2 * dof_c00, 2 * dof_c00] += B[2 * dof_c01, 2 * dof_c01] * np.cos(psi_x)
        # A[2 * dof_c00, 2 * dof_c00] += A[2 * dof_c10, 2 * dof_c10] * np.cos(psi_y)
        # B[2 * dof_c00, 2 * dof_c00] += B[2 * dof_c10, 2 * dof_c10] * np.cos(psi_y)
        # A[2 * dof_c00, 2 * dof_c00] += A[2 * dof_c11, 2 * dof_c11] * np.cos(psi_x + psi_y)
        # B[2 * dof_c00, 2 * dof_c00] += B[2 * dof_c11, 2 * dof_c11] * np.cos(psi_x + psi_y)
        #
        # # imaginary parts
        # A[2 * dofs_xmin + 1, 2 * dofs_xmin + 1] += A[2 * dofs_xmax + 1, 2 * dofs_xmax + 1] * np.sin(-1 * psi_x)
        # B[2 * dofs_xmin + 1, 2 * dofs_xmin + 1] += B[2 * dofs_xmax + 1, 2 * dofs_xmax + 1] * np.sin(-1 * psi_x)
        # A[2 * dofs_ymin + 1, 2 * dofs_ymin + 1] += A[2 * dofs_ymax + 1, 2 * dofs_ymax + 1] * np.sin(-1 * psi_y)
        # B[2 * dofs_ymin + 1, 2 * dofs_ymin + 1] += B[2 * dofs_ymax + 1, 2 * dofs_ymax + 1] * np.sin(-1 * psi_y)
        # A[2 * dof_c00 + 1, 2 * dof_c00 + 1] += A[2 * dof_c01 + 1, 2 * dof_c01 + 1] * np.sin(-1 * psi_x)
        # B[2 * dof_c00 + 1, 2 * dof_c00 + 1] += B[2 * dof_c01 + 1, 2 * dof_c01 + 1] * np.sin(-1 * psi_x)
        # A[2 * dof_c00 + 1, 2 * dof_c00 + 1] += A[2 * dof_c10 + 1, 2 * dof_c10 + 1] * np.sin(-1 * psi_y)
        # B[2 * dof_c00 + 1, 2 * dof_c00 + 1] += B[2 * dof_c10 + 1, 2 * dof_c10 + 1] * np.sin(-1 * psi_y)
        # A[2 * dof_c00 + 1, 2 * dof_c00 + 1] += A[2 * dof_c11 + 1, 2 * dof_c11 + 1] * np.sin(-1 * (psi_x + psi_y))
        # B[2 * dof_c00 + 1, 2 * dof_c00 + 1] += B[2 * dof_c11 + 1, 2 * dof_c11 + 1] * np.sin(-1 * (psi_x + psi_y))

        ones_int = identity(n_int, format='csr')
        ones_xmin = identity(n_xmin, format='csr')
        ones_xmax = identity(n_xmax, format='csr')
        ones_ymin = identity(n_ymin, format='csr')
        ones_ymax = identity(n_ymax, format='csr')

        p_rows = [*dofs_int, *dofs_xmin, *dofs_xmax, *dofs_ymin, *dofs_ymax, dof_c00, dof_c10, dof_c01, dof_c11]
        p_cols = [*dofs_int, *dofs_xmin, *dofs_xmin, *dofs_ymin, *dofs_ymin, dof_c00, dof_c00, dof_c00, dof_c00]
        #p_rows = range(self.basis.N)
        #p_cols = [*dofs_int, *dofs_xmin, *dofs_xmin, *dofs_ymin, *dofs_ymin, dof_c00, dof_c00, dof_c00, dof_c00]
        print(len(p_rows))
        print(len(p_cols))
        p_data_re = [*ones_int.data, *ones_xmin.data, *(ones_xmax.data * np.cos(psi_x)), *ones_ymin.data,
                     *(ones_ymax.data * np.cos(psi_y)), 1, np.cos(psi_x), np.cos(psi_y), np.cos(psi_x + psi_y)]
        p_data_im = [*ones_int.data, *ones_xmin.data, *(ones_xmax.data * np.sin(-psi_x)), *ones_ymin.data,
                     *(ones_ymax.data * np.sin(-psi_y)), 1, np.sin(-psi_x), np.sin(-psi_y), np.sin(-psi_x - psi_y)]
        p_data_im_conj = [*ones_int.data, *ones_xmin.data, *(ones_xmax.data * np.sin(psi_x)), *ones_ymin.data,
                          *(ones_ymax.data * np.sin(psi_y)), 1, np.sin(psi_x), np.sin(psi_y), np.sin(psi_x + psi_y)]

        P_re = coo_matrix((p_data_re, (p_rows, p_cols))).tocsr()
        P_im = coo_matrix((p_data_im, (p_rows, p_cols))).tocsr()
        P_im_conj = coo_matrix((p_data_im_conj, (p_rows, p_cols))).tocsr()

        # P_re = bmat([[ones_int, None, None, None],
        #              [None, ones_xmin, None, None],
        #              [None, ones_xmax * np.cos(psi_x), None, None],
        #              [None, None, ones_ymin, None],
        #              [None, None, ones_ymax * np.cos(psi_y), None],
        #              [None, None, None, 1],
        #              [None, None, None, np.cos(psi_x)],
        #              [None, None, None, np.cos(psi_y)],
        #              [None, None, None, np.cos(psi_x + psi_y)]], format='csr')
        #
        # P_im = bmat([[ones_int, None, None, None],
        #              [None, ones_xmin, None, None],
        #              [None, ones_xmax * np.sin(-1 * psi_x), None, None],
        #              [None, None, ones_ymin, None],
        #              [None, None, ones_ymax * np.sin(-1 * psi_y), None],
        #              [None, None, None, 1],
        #              [None, None, None, np.sin(-1 * psi_x)],
        #              [None, None, None, np.sin(-1 * psi_y)],
        #              [None, None, None, np.sin(-1 * (psi_x + psi_y))]], format='csr')
        #
        # P_im_conj = bmat([[ones_int, None, None, None],
        #                   [None, ones_xmin, None, None],
        #                   [None, ones_xmax * np.sin(psi_x), None, None],
        #                   [None, None, ones_ymin, None],
        #                   [None, None, ones_ymax * np.sin(psi_y), None],
        #                   [None, None, None, 1],
        #                   [None, None, None, np.sin(psi_x)],
        #                   [None, None, None, np.sin(psi_y)],
        #                   [None, None, None, np.sin(psi_x + psi_y)]], format='csr')

        print(f'{P_re.transpose().shape} @ {self._A_laplace_re.shape} @ {P_re.shape}')

        A_re = P_re.transpose() @ self._A_laplace_re @ P_re
        B_re = P_re.transpose() @ self._A_mass_re @ P_re
        #A_im = P_im_conj.transpose() @ self._A_laplace_im @ P_im
        #B_im = P_im_conj.transpose() @ self._A_mass_im @ P_im

        # staggered system matrices
        #A = self._stagger_re_im(A_re, A_im)
        #B = self._stagger_re_im(B_re, B_im)

        # condense system matrices to include only independent dofs
        #dofs_dep_re = 2 * np.array([*dofs_xmax, *dofs_ymax, dof_c01, dof_c10, dof_c11])
        #dofs_dep_im = 1 + 2 * dofs_dep_re

        dofs_dep_re = np.array([*dofs_xmax, *dofs_ymax, dof_c01, dof_c10, dof_c11])
        dofs_dep_im = 1 + 2 * dofs_dep_re
        #AII, BII, _, _ = skfem.condense(A_re, B_re, D=dofs_dep_re)
        #print(len(dofs_dep_re))
        #print(AII.shape)

        solver = skfem.solver_eigen_scipy(k=n_modes, sigma=1)
        #kt2, phi = skfem.solve(A_re, B_re, solver=solver)
        kt2, phi = skfem.solve(*skfem.condense(A_re, B_re, D=dofs_dep_re), solver=solver)
        print(kt2)

    def near2far(self, r: tuple, k: float, field: skfem.DiscreteField, boundaries: list):
        r_abs = np.sqrt(r[0] ** 2 + r[1] ** 2)
        s = np.sqrt(1j * k / (8 * np.pi * r_abs)) * np.exp(-1j * k * r_abs)
        return s * skfem.asm(helmholtz_near2far, self.basis.boundary(boundaries), r=r, k=k, phi=field)

    def flux_ez(self, ez: skfem.DiscreteField, boundaries: list):
        pass

    def create_gif(self, field: skfem.DiscreteField, n_frames=36):
        import matplotlib.pyplot as mplt
        from skfem.visuals.matplotlib import plot
        import matplotlib.animation as animation
        import os
        import subprocess

        if not os.path.exists('video'):
            os.mkdir('video')

        for i in range(n_frames):
            plot(self.basis, np.real(field * np.exp(1j * i / n_frames * 2 * np.pi)), figsize=(6, 6))
            mplt.tight_layout()
            mplt.savefig(f'./video/frame_{i:02d}.png', dpi=300)
            mplt.close()

        os.chdir('video')
        #subprocess.call(['ffmpeg', '-framerate', '8', '-i', 'frame_%02d.png', '-r', '30', '-pix_fmt', 'yuv420p', 'video.mp4'])

        # convert from ImageMagick
        subprocess.call(['convert', '-loop', '0', 'frame_*.png', 'video.gif'])
