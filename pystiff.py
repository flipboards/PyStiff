import sys 
import blockreader
from solver import Solver
from structure import Structure


class InputReader(object):
    def __init__(self):
        self.input_src = None
        self.input = None
        self.model = None
        self.solver = Solver()
        self.structure = Structure()

    def read(self, filename):
        self.input_src = open(filename)
        self.input = blockreader.LineReader(self.input_src)

    def run(self):
        while not self.input.eof():
            line = self.input.getline()
            print(' '.join(line.args))
            self.eval(line.header, *(line.args))

    def eval(self, command, *args):
        if command == 'read_material':
            with open(args[1]) as f:
                self.structure.material = blockreader.LineReader(f).getattrs(True)
                print('Material read:')
                print(self.structure.material)
        elif command == 'read_data':
            with open(args[1]) as f:
                self.model = blockreader.readblocks(f)
            self.structure.apply_structure(self.model)
        elif command == 'assemble':
            self.solver.assemble(self.structure, *args)
        elif command == 'apply':
            self.structure.apply(self.model, args[1])

    def close():
        self.input_src.close()


def main():
    reader = InputReader()
    reader.read(sys.argv[1])
    reader.run()
    reader.close()


if __name__ == '__main__':
    main()
