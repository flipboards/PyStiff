""" Solving a FEM problem
"""

import numpy as np
import sympy as sp
import getopt
from structure import Structure, Element
from sympy.abc import x, y, z

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
        self.elem_origins = []              # zero point of elements
        self.elem_node_ids = []
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

        self.elem_origins.clear()
        self.elem_disp_mats.clear()
        self.elem_stiff_mats.clear()
        for element in structure.elements:
            print('Begin evaluation...')

            self.elem_node_ids.append(element.node_id)

            O = _origin_2d(element.node_coord)
            print('Element original point:\n', O)

            B = self.construct_dispmat(element.type, element.node_coord - O)
            print('Displacement matrix:\n', B)

            k = self.construct_stressmat(B, structure.material)
            print('Stress matrix:\n', k)

            K = self.intergrate_2d(element.type, element.node_coord - O, k)
            print('Stiffness matrix:\n', K)

            self.elem_origins.append(O)
            self.elem_disp_mats.append(B)
            self.elem_stress_mats.append(k)
            self.elem_stiff_mats.append(K)

        self.assemble_mat(structure)
        print('Global stiffness matrix:\n', self.stiff_mat)

        self.create_force(structure)
        print('Force vector:\n', self.force)

        self.create_bc(structure)
        print('Boundary Conditions:\n', self.bc)

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
        self.restore_bc()

        for i in range(len(self.elem_disp_mats)):
            self.elem_disp_mats_solved.append(self.elem_disp_mats[i] * self.disp[self.elem_node_ids[i]])
            self.elem_stress_mats_solved.append(self.elem_stress_mats[i] * self.disp[self.elem_node_ids[i]])

    def eval(self, point):
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

        local_point = self.localize(element_id, point)
        disp = _eval_2d(self.elem_disp_mats_solved[element_id], local_point)
        stress = _eval_2d(self.elem_stress_mats_solved[element_id], local_point)
        return disp, stress

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

        B = _create_geomat_2d() * disp_mat
        S = D.dot(B)

        return B.T.dot(S)

    def intergrate_2d(self, elem_type, node_coord, f):
        """ Intergrate function f over an element.
        Args:
            elem_type: type of element;
            node_coord: Local nodes coordnation of the element;
            f: symbolic function matrix;
        Returns:
            F: Numpy matrix, definite intergral of f over element.
        """

        F = np.zeros_like(f)

        if elem_type == Element.TRIANGLE_2D:
            for i in f.shape[0]:
                for j in f.shape[1]:
                    F[i, j] = f[i, j].evalf()
            return F * _area_2d(node_coord)

        elif elem_type == Element.RECTANGLE_2D:
            for i in f.shape[0]:
                for j in f.shape[1]:
                    F[i, j] = sp.integrate(f[i, j], (x, node_coord[3,0], node_coord[0,0]), (y, node_coord[0,1], node_coord[1,1])).evalf()
            return F
        else:
            self._error_invalid_type(elem_type)

    def assemble_mat(self, structure):
        """ Assemble global matrix according to node id.
        Args:
            structure: Structure instance;
        """
        d = structure.dimension
        mat_length = d * len(structure.elements)
        self.stiff_mat = np.zeros((mat_length, mat_length))

        for element, mat in zip(structure.elements, self.elem_stiff_mats):
            for i in range(len(element.node_id)): 
                for j in range(len(element.node_id)):
                    self.stiff_mat[
                        element.node_id[i]:element.node_id[i]+d, 
                        element.node_id[j]:element.node_id[j]+d
                        ] += mat[i:i+d, j:j+d]

    def create_force(self, structure):
        """ Create force (residual) array.
        """
        self.force = np.array(structure.loads).flatten()

    def create_bc(self, structure):
        """ create bc into tuple (pos, num)
        """
        for i in range(structure.bcs):
            for j in range(structure.dimension):
                if structure.bcs[i, j]:
                    self.bc.append((i + j, structure.bcs[i, j]))

        self.bc.sort(key=lambda x:x[0])

    def modify_bc(self):
        """ Modify stiffness matrix by boundary condition.
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
                self.stiff_mat[x[0], x[0]] *= 1e20
                self.force[x[0]] *= 1e20

        else:
            raise RuntimeError('Invalid value for bc type %d' % self.bc_type)

        print('After applying BC:Stiff Mat:\n', self.stiff_mat)
        print('Force:\n', self.force)


    def restore_bc(self):
        

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
        return 0


    def localize(self, element_id, point):
        """ Return the point coords in local coordination;
        """
        return point - self.elem_origins[element_id]

    def _error_invalid_type(self, type):
        raise RuntimeError('Invalid element type %d' % type)


class Partial(object):
    """ Partial operator
    """
    def __init__(self, var):
        self.var = var

    def __mul__(self, other):
        return sp.diff(other, self.var)




def _origin_2d(node_coords):
    """ Return the origin point of an 2d polygon. Default is left down.
    """
    return np.argmax(node_coords.dot(np.array([[-1], [-1]])), axis=0)


def _create_dispmat_tri2d(node_coord):
    """ Create displacement matrix for 2d triangle without interpolation.
    Args:
        node_coord: 3x2 Matrix ((x1, y1), (x2, y2), (x3, y3))
    Returns:
        2x6 Symbolic Matrix (N1(x,y), N2(x,y), N3(x,y))
    """
    
    N = sp.zeros(2, 6)
    for i in range(3):
        a = np.outer(node_coord[(i+1)%3], node_coord[(i+2)%3])
        b = node_coord[(i+1)%3, 1] - node_coord[(i+2)%3, 1]
        c = node_coord[(i+2)%3, 0] - node_coord[(i+1)%3, 0]
        N[0, 2*i] = N[1, 2*i+1] = a + b*x + c*y

    return N / (2*_area_2d(node_coord))


def _create_dispmat_rect2d(node_coord):
    """ Create displacement matrix for 2d reactangle without interpolation.
    """

    n = [
        -(x - node_coord[3, 0])*(y - node_coord[1, 1]),
        (x - node_coord[2, 0])*(y - node_coord[0, 1]),
        -(x - node_coord[1, 0])*(y - node_coord[3, 1]),
        (x - node_coord[0, 0])*(y - node_coord[2, 1]),
        ] 

    N = sp.zeros(2, 8)
    for i in range(4):
        N[0, 2*i] = N[1, 2*i+1] = (-1)**(i-1) * (x - node_coord[3-i, 0])*(y - node_coord[1-i, 1])

    return N / _area_2d(node_coord)


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
    """

def _area_2d(node_coords):
    """ Returns the area of an 2d polygon.
    """
    return np.sum([
        np.outer(node_coords[i], node_coords[(i+1)%len(node_coords)]) for i in range(len(node_coords))
        ]) * 0.5
