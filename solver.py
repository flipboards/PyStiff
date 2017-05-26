""" Solving a FEM problem
"""

import numpy as np
import sympy as sp
import getopt
from sympy.abc import x, y, z
import multiprocessing as mp

import polygon2d
from structure import Structure, Element
from utils import str_mat, sym2float, xyz2array, mpdot


class Solver(object):
    """ problem solver
    """
    # problem type
    ELASTIC = 1
    THREMO = 2

    # bc type
    BC_REDUCE = 1
    BC_TOONES = 2
    BC_LARGENUM = 3

    def __init__(self):
        self.problem = 0
        self.bc_type = 0

        self.stiff_mat = None
        self.force = None
        self.disp = None 
        self.bc = []

        self.structure = None

        self.elem_disp_mats = []            # matrix of general displacement
        self.elem_disp_mats_solved = []     # 
        self.elem_stress_mats = []          # matrix of general stress
        self.elem_stress_mats_solved = []   # 
        self.elem_stiff_mats = []           # matrix of stiffness
        np.set_printoptions(precision=4)

    def assemble(self, structure, *args):
        """ Assemble a structure.
        Args:
            structure: Structure instance;
            args: command-line args;
                '-p [problem_type]'
                '-b [bc_type]'

        assemble works by this procedule:
            1. Create displacement matrices for elements;
            2. Create stiffness matrices for elements;
            3. Assemble stiffness matrices;
            4. Read force and boundary conditions;
        """

        self.read_assmble_args(*args)
        self.load(structure)

        for element in structure.elements:
            print('Physical analyze of one element...')

            N = self.construct_dispmat(element.type, element.node_coord)
            print('Displacement matrix:\n', str_mat(N, symbol=True))

            B, S = self.construct_stressmat(N, structure.material)
            print('Stress matrix:\n', str_mat(S, symbol=True))

            if structure.dimension == 2:
                K = self.intergrate_2d(element.type, element.node_coord, B.T.dot(S) * structure.material.thickness)

            print('Stiffness matrix:\n', str_mat(K, True))

            self.elem_disp_mats.append(N)
            self.elem_stress_mats.append(S)
            self.elem_stiff_mats.append(K)

        self.assemble_mat(structure)
        print('Global stiffness matrix:\n', str_mat(self.stiff_mat, True))

    def solve(self):
        """ Solve the problem.

            Steps:
                1. modify boundary condition;
                2. solve matrix equation;
                3. restore boundary condition;
                4. evaluate element displacement matrices;
        """
        self.modify_bc()
        self.disp = np.linalg.solve(self.stiff_mat, self.force)

        print('Raw value for disp:\n%s' % str_mat(self.disp))

        self.restore_bc()
        print('Node displacement:\n%s' % str_mat(self.disp))

        def expand_node_id(node_id):
            return np.array([[i*2, i*2+1] for i in node_id]).flatten()
            

        for i in range(len(self.elem_disp_mats)):
            node_id_in_mat = expand_node_id(self.structure.elements[i].node_id)

            self.elem_disp_mats_solved.append(self.elem_disp_mats[i].dot(self.disp[node_id_in_mat]))
            self.elem_stress_mats_solved.append(self.elem_stress_mats[i].dot(self.disp[node_id_in_mat]))

    def read_assmble_args(self, *args):
        """ Read assemble arguments as standard CLI args.
        """
        opts, _args = getopt.getopt(args[1:], 'p:b:')
        for opt, arg in opts:
            if opt == '-p': # problem type
                if arg == 'elastic':
                    self.problem = Solver.ELASTIC
                else:
                    raise ValueError(arg)

            elif opt == '-b':
                if arg == 'reduce':
                    self.bc_type = Solver.BC_REDUCE
                elif arg == 'toones':
                    self.bc_type = Solver.BC_TOONES
                elif arg == 'largenum':
                    self.bc_type = Solver.BC_LARGENUM
                else:
                    raise ValueError(arg)

    def eval(self, *args):

        def print_eval_result(disp, stress):
            print('Displacement matrix:\n%s' % str_mat(disp))
            print('Stress matrix:\n%s' % str_mat(stress))

        opts, _args = getopt.getopt(args[1:], 'p:b:')
        for opt, arg in opts:
            if opt == '-p': # problem type
                print('Point:', arg)
                print_eval_result(*self.eval_point(xyz2array(arg)))

    def load(self, structure):
        """ Loading structure.
        """
        self.structure = structure

        self.create_force(self.structure)
        print('Force vector:\n', str_mat(self.force))

        self.create_bc(self.structure)
        print('Boundary Conditions:\n', '\n'.join(['%d: % .4e' % s for s in self.bc]))

        
    def create_force(self, structure):
        """ Create force (residual) array.
        """
        self.force = np.array(structure.loads).flatten()

    def create_bc(self, structure):
        """ create bc into tuple (pos, num)
        """
        for i in range(len(structure.bcs)):
            for j, bc in enumerate(structure.bcs[i]):
                if bc is not None:
                    self.bc.append((structure.dimension * i + j, bc))

        self.bc.sort(key=lambda x:x[0])

    def construct_dispmat(self, elem_type, node_coord):
        """ Construct displacement matrix of an element.
        Args:
            elem_type: Type of element;
            node_coord: Local coord matrix of element.
        Returns:
            disp_mat: Displacement matrix of the element.
        """
        if elem_type == Element.TRIANGLE_2D:
            return _create_dispmat_tri2d(node_coord)
        elif elem_type == Element.RECTANGLE_2D:
            return _create_dispmat_rect2d(node_coord)
        else:
            self._error_invalid_type(elem_type)

    def construct_stressmat(self, disp_mat, material):
        """ Construct stress matrix of an element.
        Args:
            disp_mat: Displacement matrix of the element;
            material: Material of the element;
        Returns:
            k: Array of symbolic function of stress;
        """
        if self.problem == Solver.ELASTIC:
            D = _create_physicalmat_2d_normal(material.poisson, material.elastic_modulus)
        else:
            raise RuntimeError('Invalid value for problem: %d' % self.problem)

        B = np.dot(_create_geomat_2d(), disp_mat.tolist())
        print('Creating stress matrix...')
        return B, D.dot(B) 

    def intergrate_2d(self, elem_type, node_coord, f):
        """ Intergrate function f over an element.
        Args:
            elem_type: type of element;
            node_coord: Local nodes coordnation of the element;
            f: symbolic function matrix;
        Returns:
            F: Numpy matrix, definite intergral of f over element.
        """

        F = np.zeros_like(f, dtype=np.float)

        if elem_type == Element.TRIANGLE_2D:
            for i in range(f.shape[0]):
                for j in range(f.shape[1]):
                    F[i, j] = sym2float(f[i, j].evalf())
            return F * polygon2d._area_2d(node_coord)

        elif elem_type == Element.RECTANGLE_2D:
            
            pool = mp.Pool(4)

            tmp = np.zeros_like(F, dtype=object)

            for i in range(f.shape[0]):
                for j in range(f.shape[1]):
                    tmp[i, j] = pool.apply_async(_integate_2d, args=(f[i, j], (x, node_coord[3,0], node_coord[0,0]), (y, node_coord[0,1], node_coord[1,1])))

            pool.close()
            pool.join()

            for i in range(f.shape[0]):
                for j in range(f.shape[1]):
                    F[i, j] = sym2float(tmp[i, j].get())
           
            return F
        else:
            self._error_invalid_type(elem_type)

    def assemble_mat(self, structure):
        """ Assemble global matrix according to node id.
        Args:
            structure: Structure instance;
        """
        d = structure.dimension
        mat_length = d * len(structure.coords)
        self.stiff_mat = np.zeros((mat_length, mat_length))

        for element, mat in zip(structure.elements, self.elem_stiff_mats):
            for i in range(len(element.node_id)): 
                for j in range(len(element.node_id)):
                    self.stiff_mat[
                        d*element.node_id[i]:d*element.node_id[i]+d, 
                        d*element.node_id[j]:d*element.node_id[j]+d
                        ] += mat[d*i:d*i+d, d*j:d*j+d]

    def modify_bc(self):
        """ Modify stiffness matrix by boundary condition according to bc_type.
        """
        
        fix_disp = np.array([x[1] for x in self.bc])
        remain_vec = np.array([True] * len(self.force)) # vector kept after reduction
        for x in self.bc:
            remain_vec[x[0]] = False

        if self.bc_type == Solver.BC_REDUCE:

            self.force = self.force[remain_vec] - self.stiff_mat[remain_vec][:,np.logical_not(remain_vec)].dot(fix_disp)
            self.stiff_mat = self.stiff_mat[remain_vec][:,remain_vec]

        elif self.bc_type == Solver.BC_TOONES:
            
            self.force -= self.stiff_mat[:,np.logical_not(remain_vec)].dot(fix_disp)
            for x in self.bc:
                self.force[x[0]] = x[1]
                self.stiff_mat[x[0]] *= 0
                self.stiff_mat[:,x[0]] *= 0
                self.stiff_mat[x[0], x[0]] = 1


        elif self.bc_type == Solver.BC_LARGENUM:
            for x in self.bc:
                self.force[x[0]] = x[1] * self.stiff_mat[x[0], x[0]] *1e20
                self.stiff_mat[x[0], x[0]] *= 1e20

        else:
            raise RuntimeError('Invalid value for bc type %d' % self.bc_type)

        print('After applying BC:\n  Stiff Mat:\n%s' % str_mat(self.stiff_mat, True))
        print('  Force:\n%s' % str_mat(self.force))


    def restore_bc(self):
        """ Restore from bounary condition.
        """

        if self.bc_type == Solver.BC_REDUCE:
            real_disp = []

            if len(self.bc) > 1:
                for i in range(len(self.bc)):
                    real_disp.append(self.bc[i][1])
                    real_disp += self.disp[self.bc[i][0]:self.bc[i+1][0]]

            real_disp += self.bc[-1]
        else:
            for idx, pos in self.bc:
                self.disp[idx] = pos


    def locate_element(self, point):
        """ Locate an element.
        """
        for i in range(len(self.structure.elements)):
            if polygon2d._within_2d(point, self.structure.elements[i].node_coord):
                return i
        return None

    def eval_point(self, point):
        """ Evaluate displacement matrix in arbitary position.
        Args:
            point: coord Array (x, y, ..)
        Returns:
            disp: general displacement Array (x, y, ..)
            stress: general strain Array (Fx, Fy, .. )
        """
        element_id = self.locate_element(point)

        if not self.elem_disp_mats_solved:
            raise RuntimeError('Problem not solved')

        disp = _eval_2d(self.elem_disp_mats_solved[element_id], point)
        stress = _eval_2d(self.elem_stress_mats_solved[element_id], point)
        return disp, stress

    def _error_invalid_type(self, type):
        raise RuntimeError('Invalid element type %d' % type)


class Partial(object):
    """ Partial operator
    """
    def __init__(self, var):
        self.var = var

    def __mul__(self, other):
        return sp.diff(other, self.var)



def _create_dispmat_tri2d(node_coord):
    """ Create displacement matrix for 2d triangle without interpolation.
    Args:
        node_coord: 3x2 Matrix ((x1, y1), (x2, y2), (x3, y3))
    Returns:
        2x6 Symbolic Matrix (N1(x,y), N2(x,y), N3(x,y))
    """
    
    N = sp.zeros(2, 6)

    a = np.roll(node_coord[:,0], -1) * np.roll(node_coord[:,1], -2) - np.roll(node_coord[:,0], -2) * np.roll(node_coord[:,1], -1)
    b = np.roll(node_coord[:,1], -1) - np.roll(node_coord[:,1], -2)
    c = np.roll(node_coord[:,0], -2) - np.roll(node_coord[:,0], -1)

    for i in range(3):
        N[0, 2*i] = N[1, 2*i+1] = a[i] + b[i]*x + c[i]*y

    return N / (2*polygon2d._area_2d(node_coord))


def _create_dispmat_rect2d(node_coord):
    """ Create displacement matrix for 2d reactangle without interpolation.
    """

    N = sp.zeros(2, 8)
    for i in range(4):
        N[0, 2*i] = N[1, 2*i+1] = (-1)**(i-1) * (x - node_coord[3-i, 0])*(y - node_coord[1-i, 1])

    return N /polygon2d._area_2d(node_coord)


def _create_geomat_2d():
    """ Geometry matrix of 2d element.
    """
    return np.array([[Partial(x), 0], [0, Partial(y)], [Partial(y), Partial(x)]])


def _create_physicalmat_2d_normal(p, E):
    """ Create Physical matrix for normal 2d element.
    Args:
        p: Poisson ratio;
        E: Elastic modulus
    Returns:
        3x3 2d strain Matrix
    """
    return E/(1 - p**2)*np.array([
        [1, p, 0], [p, 1, 0], [0, 0, (1-p)/2]
        ])


def _eval_2d(disp_expr, coord):
    """ Evaluate an 2d displacement
    Args:
        disp_expr: Array (u(x), v(x))
        coord: Array (x, y)
    Returns:
        Float array (u, v)
    """
    return np.array([sym2float(i.evalf(subs={x:coord[0], y:coord[1]})) for i in disp_expr])


def _integate_2d(f, xrange, yrange):
    """ Subroutine of 2d integration
    """
    return sp.integrate(f, xrange, yrange).evalf()

