"""Split PROTACs into warhead / linker / E3 ligand with mapped attachment points.

Architecture is linear: warhead - linker - E3 ligand.
  1. locate the E3 ligand by a reference scaffold, cut the bridge bond that yields
     the SMALLEST E3-side fragment still containing that scaffold;
  2. from the exposed atom, walk outward through linker-like atoms (acyclic, or
     small saturated/triazole rings) until a warhead ring system is reached.

Emits the linker as the core with [*:1] (warhead end) and [*:2] (E3 end).
"""
import sys
from rdkit import Chem, RDLogger

RDLogger.DisableLog('rdApp.*')

E3_SCAFFOLDS = [
    # CRBN / IMiD
    ('CRBN', 'thalidomide',  'O=C1CCC(N2C(=O)c3ccccc3C2=O)C(=O)N1'),
    ('CRBN', 'lenalidomide', 'O=C1CCC(N2Cc3ccccc3C2=O)C(=O)N1'),
    ('CRBN', 'aza-IMiD',     'O=C1CCC(N2C(=O)c3ccncc3C2=O)C(=O)N1'),
    ('CRBN', 'aza-IMiD',     'O=C1CCC(N2C(=O)c3cccnc3C2=O)C(=O)N1'),
    ('CRBN', 'aza-lenalid',  'O=C1CCC(N2Cc3ccncc3C2=O)C(=O)N1'),
    # VHL (VH032 family) - keep the tert-leucine so the exit vector is its amine
    ('VHL',  'VH032',        'CC(C)(C)C(N)C(=O)N1CC(O)CC1C(=O)NCc1ccc(-c2scnc2C)cc1'),
    ('VHL',  'VH032-Me',     'CC(C)(C)C(N)C(=O)N1CC(O)CC1C(=O)NC(C)c1ccc(-c2scnc2C)cc1'),
]
E3_QUERIES = [(fam, name, Chem.MolFromSmiles(s)) for fam, name, s in E3_SCAFFOLDS]
assert all(q is not None for _, _, q in E3_QUERIES)


def ring_info(mol):
    """Per-atom ring membership plus which rings may be walked through as linker."""
    ri = mol.GetRingInfo()
    rings = [set(r) for r in ri.AtomRings()]
    # fused rings share atoms; a fused system belongs to the warhead
    fused = set()
    for i, a in enumerate(rings):
        for b in rings[i + 1:]:
            if a & b:
                fused |= a | b
    passable = set()
    for r in rings:
        if r & fused or len(r) > 7:
            continue
        atoms = [mol.GetAtomWithIdx(i) for i in r]
        if all(not a.GetIsAromatic() for a in atoms):
            passable |= r                                    # piperazine, piperidine, cyclohexane
        elif len(r) == 5 and sum(a.GetSymbol() == 'N' for a in atoms) >= 2:
            passable |= r                                    # click triazole / tetrazole
    in_ring = {a.GetIdx() for r in rings for a in (mol.GetAtomWithIdx(i) for i in r)}
    return in_ring, passable


def bridge_sides(mol, bond):
    """Atom sets on each side of an acyclic bond."""
    a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
    seen, stack = {a}, [a]
    while stack:
        cur = stack.pop()
        for nb in mol.GetAtomWithIdx(cur).GetNeighbors():
            i = nb.GetIdx()
            if i in seen or (cur == a and i == b):
                continue
            seen.add(i)
            stack.append(i)
    return seen, set(range(mol.GetNumAtoms())) - seen


def cuttable(bond):
    return (not bond.IsInRing() and bond.GetBondType() == Chem.BondType.SINGLE
            and bond.GetBeginAtom().GetAtomicNum() > 1 and bond.GetEndAtom().GetAtomicNum() > 1)


def find_e3(mol):
    """(bond index, e3 atom set, family, scaffold name) for the tightest E3 fragment."""
    best = None
    for fam, name, q in E3_QUERIES:
        match = mol.GetSubstructMatch(q)
        if not match:
            continue
        core = set(match)
        for bond in mol.GetBonds():
            if not cuttable(bond):
                continue
            s1, s2 = bridge_sides(mol, bond)
            side = s1 if core <= s1 else (s2 if core <= s2 else None)
            if side is None or len(side) == mol.GetNumAtoms():
                continue
            if best is None or len(side) < len(best[1]):
                best = (bond.GetIdx(), side, fam, name)
    return best


def walk_linker(mol, start, in_ring, passable, forbidden):
    """Atoms reachable from `start` without entering a warhead ring system."""
    if start in in_ring and start not in passable:
        return None
    linker, stack = {start}, [start]
    while stack:
        cur = stack.pop()
        for nb in mol.GetAtomWithIdx(cur).GetNeighbors():
            i = nb.GetIdx()
            if i in linker or i in forbidden:
                continue
            if i in in_ring and i not in passable:
                continue
            linker.add(i)
            stack.append(i)
    return linker


def split(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumAtoms() > 140:
        return None
    if len(Chem.GetMolFrags(mol)) != 1:
        return None
    found = find_e3(mol)
    if not found:
        return None
    b1, e3_atoms, fam, name = found
    bond1 = mol.GetBondWithIdx(b1)
    start = bond1.GetEndAtomIdx() if bond1.GetBeginAtomIdx() in e3_atoms else bond1.GetBeginAtomIdx()

    in_ring, passable = ring_info(mol)
    linker = walk_linker(mol, start, in_ring, passable, e3_atoms)
    if linker is None or not 2 <= len(linker) <= 26:
        return None

    exits = [b for b in mol.GetBonds()
             if cuttable(b) and len({b.GetBeginAtomIdx(), b.GetEndAtomIdx()} & linker) == 1
             and b.GetIdx() != b1]
    if len(exits) != 1:
        return None                                          # branched or ambiguous tether
    bond2 = exits[0]
    b2 = bond2.GetIdx()
    warhead = set(range(mol.GetNumAtoms())) - e3_atoms - linker
    if len(warhead) < 12:
        return None

    war_mol_rings = sum(1 for r in mol.GetRingInfo().AtomRings() if set(r) <= warhead)
    if war_mol_rings < 2:
        return None

    labels = [(2, 2), (1, 1)] if b1 < b2 else [(1, 1), (2, 2)]
    frag = Chem.FragmentOnBonds(mol, sorted([b1, b2]), addDummies=True, dummyLabels=labels)
    pieces = Chem.GetMolFrags(frag, asMols=True, sanitizeFrags=True)
    if len(pieces) != 3:
        return None

    out = {}
    for p in pieces:
        iso = sorted(a.GetIsotope() for a in p.GetAtoms() if a.GetAtomicNum() == 0)
        for a in p.GetAtoms():
            if a.GetAtomicNum() == 0:
                a.SetAtomMapNum(a.GetIsotope())
                a.SetIsotope(0)
        smi = Chem.MolToSmiles(p)
        if iso == [1, 2]:
            out['linker'] = smi
        elif iso == [2]:
            out['e3'] = smi
        elif iso == [1]:
            out['warhead'] = smi
    if len(out) != 3:
        return None
    out['e3_family'] = fam
    out['e3_scaffold'] = name
    return out


if __name__ == '__main__':
    for smi in sys.argv[1:]:
        print(smi, '->', split(smi))
