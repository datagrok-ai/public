import {RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {_isSmarts, getMolSafe} from '../utils/mol-creation_rdkit';

/** Prepare coordinates for drawing only; the source reaction and its atom maps stay unchanged. */
export function prepareReactionDepiction(
  rdkit: RDModule, reaction: string, molblocks: Map<string, string>,
): string | null {
  const sides = reaction.split(/(?<!-)>/);
  if (sides.length !== 3 || reaction.includes('|') || reaction.startsWith('$RXN'))
    return null;

  const molecules: {mol: RDMol, side: number, key: string, positioned: boolean}[] = [];
  try {
    if (sides.some((s) => s && _isSmarts(s)))
      return null;
    for (let side = 0; side < sides.length; side++) {
      if (!sides[side])
        continue;
      const {mol, isQMol} = getMolSafe(sides[side], {removeHs: false}, rdkit);
      try {
        if (!mol || isQMol)
          return null;
        const {molList} = mol.get_frags();
        try {
          for (let i = 0; i < molList.size(); i++) {
            const fragment = molList.at(i);
            const entry = {mol: fragment, side, key: '', positioned: fragment.has_coords() === 2};
            molecules.push(entry);
            entry.key = fragment.get_smiles();
            const previous = molblocks.get(entry.key);
            if (previous && !entry.positioned) {
              const template = rdkit.get_mol(previous, '{"removeHs":false}');
              try {
                fragment.generate_aligned_coords(template, '{"useCoordGen":false,"acceptFailure":false}');
                entry.positioned = true;
              } finally {
                template.delete();
              }
            }
            if (!entry.positioned)
              fragment.set_new_coords(false);
          }
        } finally {
          molList.delete();
        }
      } finally {
        mol?.delete();
      }
    }

    // Prefer a substantial retained substrate over a smaller reagent or byproduct.
    const templates = molecules.filter((m) => m.side !== 1 && m.mol.get_num_atoms(true) >= 3 &&
      molecules.some((other) => other.side === 2 - m.side &&
        other.mol.get_substruct_match(m.mol) !== '{}'))
      .sort((a, b) => b.mol.get_num_atoms(true) - a.mol.get_num_atoms(true));
    for (const entry of molecules) {
      if (entry.positioned || entry.side === 1)
        continue;
      const template = templates.find((t) => t !== entry && entry.mol.get_substruct_match(t.mol) !== '{}');
      if (template)
        entry.mol.generate_aligned_coords(template.mol, '{"useCoordGen":false,"acceptFailure":false}');
    }

    const smarts: string[][] = [[], [], []];
    const coordinates: string[] = [];
    for (const entry of molecules) {
      const cx = entry.mol.get_cxsmarts();
      const match = /^(.*?) \|\(([^)]*)\)(.*)\|$/.exec(cx);
      // Maps and ordinary stereobonds are already encoded in SMARTS; keep other CX metadata on the original path.
      if (!match || !/^(?:,w[UD]:[\d.,]+|,atomProp:(?:\d+\.molAtomMapNumber\.\d+:?)+)*$/.test(match[3]))
        return null;
      smarts[entry.side].push(match[1]);
      coordinates.push(match[2]);
    }
    for (const entry of molecules)
      molblocks.set(entry.key, entry.mol.get_molblock());
    return `${smarts.map((s) => s.join('.')).join('>')} |(${coordinates.join(';')})|`;
  } catch {
    return null;
  } finally {
    for (const entry of molecules)
      entry.mol.delete();
  }
}
