import os
import sys

ref, new = sys.argv[1], sys.argv[2]


def atoms(path):
    res = []
    for line in open(path):
        if line.startswith(('ATOM', 'HETATM')):
            res.append((line[12:16].strip(), float(line[30:38]), float(line[38:46]), float(line[46:54]),
                        float(line[70:76]), line[77:].strip()))
    return res


def numbers(path):
    res = []
    for line in open(path):
        try:
            res.append(float(line))
        except ValueError:
            pass
    return res


for name in sorted(os.listdir(ref)):
    a, b = os.path.join(ref, name), os.path.join(new, name)
    if not os.path.exists(b):
        print('compare: %s MISSING in new' % name)
        continue
    same = open(a).read() == open(b).read()
    if same:
        print('compare: %s identical' % name)
    elif name.endswith('.pdbqt'):
        x, y = atoms(a), atoms(b)
        names_eq = [t[0] for t in x] == [t[0] for t in y]
        types_eq = [t[5] for t in x] == [t[5] for t in y]
        dq = max([abs(p[4] - q[4]) for p, q in zip(x, y)] or [0])
        dxyz = max([max(abs(p[i] - q[i]) for i in (1, 2, 3)) for p, q in zip(x, y)] or [0])
        print('compare: %s atoms %d/%d names_eq=%s types_eq=%s max_dcharge=%.4f max_dxyz=%.4f'
              % (name, len(x), len(y), names_eq, types_eq, dq, dxyz))
        bx, by = dict((t[0], t) for t in x), dict((t[0], t) for t in y)
        print('compare: %s by name: same_names=%s same_atoms=%s' % (name, sorted(bx) == sorted(by),
              all(bx[k] == by.get(k) for k in bx)))
        for k in bx:
            if bx[k] != by.get(k):
                print('compare:   %s ref=%s new=%s' % (k, bx[k], by.get(k)))
        if name.startswith('ligand'):
            print('compare: --- ref ' + name)
            print(open(a).read())
            print('compare: --- new')
            print(open(b).read())
    elif name.endswith('.map'):
        x, y = numbers(a), numbers(b)
        d = max([abs(p - q) for p, q in zip(x, y)] or [0])
        print('compare: %s values %d/%d max_diff=%.5f' % (name, len(x), len(y), d))
    else:
        print('compare: %s differs' % name)
