""" Read block file, config file, command lines.
"""

class CommandLine(object):
    def __init__(self, command, args):
        self.header = command
        self.args = args


class Attr(object):
    """ Indexable object 
    """
    def __init__(self, d):
        self.d = d

    def __getattr__(self, s):
        try:
            return self.d[s]
        except KeyError:
            raise RuntimeError('Attr has no attribute named %s' % s)

    def __str__(self, **kwargs):
        return '\n'.join([k + '=' + str(v) for k, v in self.d.items()])


class LineReader(object):
    """ Read linear command from file
    """

    def __init__(self, src):
        self.lines = src.readlines()
        self.line_cnt = 0

    def getline(self):
        ls = cleanline(self.lines[self.line_cnt]).split()
        self.line_cnt += 1
        return CommandLine(ls[0], ls)

    def getattrs(self, to_float=False, split_sym='='):
        attrs = {}
        for line in self.lines:
            ls = cleanline(line).split(split_sym)
            if to_float:
                attrs[ls[0]] = float(ls[1])
            else:
                attrs[ls[0]] = ls[1]

        return Attr(attrs)

    def eof(self):
        return self.line_cnt == len(self.lines)


def cleanline(line):
    return line.rstrip('\n').strip()


def readblocks(src, begin_str='BEGIN', end_str='END'):
    blocks = {}

    crt_lines = []
    title = None

    for line in src.readlines():
        line = cleanline(line)
        if line == 'END':
            append_state = False
            blocks[title] = crt_lines
            crt_lines = []
        elif line.split()[0] == 'BEGIN':
            title = line.split()[1]
        elif title:
            crt_lines.append(line.split())
           
    return Attr(blocks)