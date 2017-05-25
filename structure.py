import numpy as np
from utils import xyz2array

class Element(object):
    """ Basic element in FEM.
    """
    # type
    TRUSS_1D = 1
    TRUSS_2D = 2
    TRIANGLE_2D = 3
    RECTANGLE_2D = 4

    NODE_NUM = [None, 2, 2, 3, 4]
    COORD_NUM = [None, 1, 2, 2, 2]

    def __init__(self, elem_type, node_id, node_coord):
        """ Construct a basic element in FEM.
        Args:
            elem_type: type of element (int)
            node_id: Array, element id -> structure id
            node_coord: Matrix of node coords
        """
        self.type = elem_type
        if not len(node_id) == Element.NODE_NUM[self.type]:
            raise RuntimeError('Node number not match')

        self.node_id = node_id
        self.node_coord = node_coord


class Structure(object):

    def __init__(self):
        self.material = None
        self.elements = []
        self.coords = []
        self.node_id = {}
        self.bcs = []
        self.loads = []
        self.dimension = 0

    def apply_structure(self, model):
        for line in model.COORD:
            self.coords.append(xyz2array(line[1]))
            self.node_id[line[0]] = len(self.coords) - 1

        self.dimension = len(self.coords[0])

        print('Node positions:')
        for node in self.node_id:
            print('Node %s: %s'% (node, str(self.coords[self.node_id[node]])))

        element_ids = []

        for line in model.ELEM:
            element_ids.append(line)

        self.fill_elements(element_ids) # push elements

        print('Node id for elements:')
        for idx, element in enumerate(self.elements):
            print(idx, 'type=%d' % element.type, element.node_id)

    def fill_elements(self, element_ids):
        for nodes in element_ids:
            self.elements.append(Element(
                self.get_element_type(len(nodes)),
                np.array([self.node_id[x] for x in nodes]),
                np.array([self.coords[self.node_id[x]] for x in nodes])
                ))

    def get_element_type(self, node_num):
        if self.dimension == 1:
            if node_num == 2:
                return Element.TRUSS_1D
        if self.dimension == 2:
            if node_num == 2:
                return Element.TRUSS_2D
            elif node_num == 3:
                return Element.TRIANGLE_2D
            elif node_num == 4:
                return Element.RECTANGLE_2D

    def apply(self, model, option):
        if option == 'load':
            self.apply_load(model)
        elif option == 'bc':
            self.apply_bc(model)

    def apply_load(self, model):

        for i in range(len(self.coords)):
            self.loads.append(np.zeros(self.dimension))

        for line in model.LOAD:
            elem_num = len(line) - 1
            for e in line[:-1]:
                self.loads[self.node_id[e]] += xyz2array(line[-1])/elem_num

        print('Load of each node:')
        for node in self.node_id:
            print(node, self.loads[self.node_id[node]])

    def apply_bc(self, model):
        self.bcs = np.array([[None] * self.dimension] * len(self.coords))

        for line in model.BC:
            for e in line[:-1]:
                for c_idx, c_data in enumerate(line[-1].split(',')):
                    if c_data != 'x':
                        self.bcs[self.node_id[e]][c_idx] = float(c_data)           

        print('BC of each node:')
        for node in self.node_id:
            print(node, self.bcs[self.node_id[node]])


