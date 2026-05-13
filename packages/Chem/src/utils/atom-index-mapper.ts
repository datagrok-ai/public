/**
 * Atom-index mapper: maps 2D atom indices (RDKit canonical, from the
 * interactive atom picker) to 3D atom indices (PDB serial numbers, from
 * a docked pose or 3D coordinate file).
 *
 * The core problem: the 2D SMILES and the 3D docked pose represent the
 * same molecule but with different atom orderings. AutoDock's PDBQT
 * format reorganises atoms by torsion tree (ROOT/BRANCH), and 3D
 * structures add explicit hydrogens. A naive index+1 offset doesn't
 * work — we need a topology-based substructure match.
 *
 * Four-tier fallback chain:
 *   1.  Strict substructure match (works for SDF/molblock with bond orders)
 *   2a. SMARTS bond-agnostic match (any-bond query, preserves atom types)
 *   2b. Flattened molblock match (all single bonds, less precise)
 *   3.  Heavy-atom serial order (last resort when topology matching fails)
 */

import {RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';

// -- PDB/PDBQT → Molblock converter -----------------------------------------

export interface PdbAtom {
  serial: number;
  element: string;
  x: number; y: number; z: number;
}

/** Covalent radii (Å) for common organic elements. Sum of two radii + 0.4 Å
 *  tolerance gives the max bonding distance for that pair. */
const COVALENT_RADII: {[el: string]: number} = {
  H: 0.31, D: 0.31, C: 0.76, N: 0.71, O: 0.66, F: 0.57,
  P: 1.07, S: 1.05, Cl: 1.02, Br: 1.20, I: 1.39, Se: 1.20,
  Si: 1.11, B: 0.84,
};
const BOND_TOLERANCE = 0.45;

/** Parses ATOM/HETATM lines from PDB or PDBQT text. */
export function parsePdbAtoms(pdbText: string): PdbAtom[] {
  const atoms: PdbAtom[] = [];
  for (const line of pdbText.split('\n')) {
    const rec = line.substring(0, 6).trim();
    if (rec !== 'ATOM' && rec !== 'HETATM') continue;
    // PDB serial number: columns 7-11 (1-indexed). This is what Molstar
    // uses as the atom 'id' for overpaint queries.
    const pdbSerial = parseInt(line.substring(6, 11).trim(), 10);
    const x = parseFloat(line.substring(30, 38));
    const y = parseFloat(line.substring(38, 46));
    const z = parseFloat(line.substring(46, 54));
    // Element symbol: columns 77-78 in standard PDB. For PDBQT, might
    // be in columns 77-78 or can be derived from atom name (cols 13-16).
    // PDBQT atom type normalization map.
    const pdbqtMap: {[k: string]: string} = {
      A: 'C', OA: 'O', NA: 'N', SA: 'S', HD: 'H', HS: 'H',
    };
    let element = line.length >= 78 ? line.substring(76, 78).trim() : '';
    // Normalize PDBQT types that appear in columns 77-78.
    if (element && pdbqtMap[element]) element = pdbqtMap[element];
    if (!element || /\d/.test(element)) {
      // Derive from atom name (cols 13-16): strip digits, take letters.
      const atomName = line.substring(12, 16).trim();
      element = atomName.replace(/[0-9]/g, '').trim();
      if (pdbqtMap[element]) element = pdbqtMap[element];
      if (element.length >= 2)
        element = element[0].toUpperCase() + element[1].toLowerCase();
      else
        element = element.toUpperCase();
    }
    if (isNaN(x) || isNaN(y) || isNaN(z) || !element) continue;
    // Use the actual PDB serial if valid, otherwise fallback to sequential.
    const serial = (!isNaN(pdbSerial) && pdbSerial > 0) ? pdbSerial : atoms.length + 1;
    atoms.push({serial, element, x, y, z});
  }
  return atoms;
}

/** Infers bonds from interatomic distances using covalent radii.
 *  Returns [i, j] pairs (0-based). */
export function inferBonds(atoms: PdbAtom[]): [number, number][] {
  const bonds: [number, number][] = [];
  const defaultR = 0.77; // default covalent radius for unknown elements
  for (let i = 0; i < atoms.length; i++) {
    const ri = COVALENT_RADII[atoms[i].element] ?? defaultR;
    for (let j = i + 1; j < atoms.length; j++) {
      const rj = COVALENT_RADII[atoms[j].element] ?? defaultR;
      const dx = atoms[i].x - atoms[j].x;
      const dy = atoms[i].y - atoms[j].y;
      const dz = atoms[i].z - atoms[j].z;
      const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
      const maxDist = ri + rj + BOND_TOLERANCE;
      if (dist > 0.4 && dist <= maxDist)
        bonds.push([i, j]);
    }
  }
  return bonds;
}

/** Builds a V2000 molblock from parsed PDB atoms.
 *  Atom order is preserved so that molblock index i = PDB atom index i. */
export function pdbAtomsToMolblock(atoms: PdbAtom[]): string | null {
  if (atoms.length === 0) return null;

  const bonds = inferBonds(atoms);

  const lines: string[] = [
    '',
    '     RDKit          3D',
    '',
  ];

  const nAtoms = atoms.length.toString().padStart(3);
  const nBonds = bonds.length.toString().padStart(3);
  lines.push(`${nAtoms}${nBonds}  0  0  0  0  0  0  0  0999 V2000`);

  for (const a of atoms) {
    const xs = a.x.toFixed(4).padStart(10);
    const ys = a.y.toFixed(4).padStart(10);
    const zs = a.z.toFixed(4).padStart(10);
    const el = a.element.padEnd(3);
    lines.push(`${xs}${ys}${zs} ${el} 0  0  0  0  0  0  0  0  0  0  0  0`);
  }

  for (const [i, j] of bonds) {
    const a1 = (i + 1).toString().padStart(3);
    const a2 = (j + 1).toString().padStart(3);
    lines.push(`${a1}${a2}  1  0`);
  }

  lines.push('M  END');
  return lines.join('\n');
}

/**
 * Result of mapping 2D atom indices to 3D atom indices.
 *
 * `mapping[i]` = the 0-based atom index in the 3D molecule that
 * corresponds to 2D atom index `i`. If the mapping couldn't be
 * determined for atom `i`, `mapping[i]` is -1.
 */
export interface AtomIndexMapping {
  /** 2D atom index → 3D atom index (0-based, heavy atoms only).
   *  mapping[i] = 3D index for 2D atom i, or -1 if unmapped. */
  mapping: number[];
  /** Which tier of the fallback chain produced this mapping. */
  method: 'substruct' | 'heavy-atom-order';
  /** Number of successfully mapped atoms. */
  mappedCount: number;
  /** Optional: actual PDB serial numbers from the file (1-based).
   *  If present, use pdbSerials[mapping[i]] instead of mapping[i]+1
   *  to get the Molstar atom 'id'. */
  pdbSerials?: number[];
}

/**
 * Maps 2D heavy-atom indices to 3D heavy-atom indices using RDKit's
 * substructure matching. Falls back through progressively looser
 * matching strategies if strict matching fails.
 *
 * @param rdkit      The RDKit module instance.
 * @param smiles2D   The 2D molecule (SMILES or molblock) — typically
 *                   from the grid cell the user picked atoms in.
 * @param mol3DStr   The 3D molecule (PDB, PDBQT converted to PDB, or
 *                   SDF) — typically from the docked-pose column.
 * @returns          The index mapping, or null if the molecules can't
 *                   be parsed at all.
 */
export function mapAtomIndices2Dto3D(
  rdkit: RDModule,
  smiles2D: string,
  mol3DStr: string,
): AtomIndexMapping | null {
  let mol2D: RDMol | null = null;
  let mol3D: RDMol | null = null;
  // Track whether mol3D was built from a PDB-derived molblock so we know
  // that match indices correspond directly to PDB serial - 1.
  let mol3DFromPdb = false;
  // If we had to strip H from the 3D mol for sanitization, keep a map
  // from the H-stripped index back to the original PDB index.
  let heavyToPdbIdx: number[] | null = null;
  let pdbAtomsParsed: PdbAtom[] | null = null;

  try {
    mol2D = rdkit.get_mol(smiles2D);
    if (!mol2D || !mol2D.is_valid()) return null;

    // Try parsing the 3D string directly first (works for SDF/molblock).
    mol3D = rdkit.get_mol(mol3DStr);

    // If direct parsing fails, try converting PDB/PDBQT to molblock.
    if (!mol3D || !mol3D.is_valid()) {
      mol3D?.delete();
      mol3D = null;
      pdbAtomsParsed = parsePdbAtoms(mol3DStr);
      if (pdbAtomsParsed.length > 0) {
        const molblock = pdbAtomsToMolblock(pdbAtomsParsed);
        if (molblock) {
          // Try parsing with progressively looser options.
          const tryParse = (opts: object): RDMol | null => {
            const m = rdkit.get_mol(molblock, JSON.stringify(opts));
            if (m && !m.is_valid()) {
              m.delete();
              return null;
            }
            return m;
          };
          mol3D = tryParse({sanitize: false, removeHs: false}) ??
            tryParse({sanitize: true, removeHs: false});

          // Last resort: strip H but track the index mapping.
          if (!mol3D) {
            heavyToPdbIdx = [];
            for (let i = 0; i < pdbAtomsParsed.length; i++) {
              if (pdbAtomsParsed[i].element !== 'H' && pdbAtomsParsed[i].element !== 'D')
                heavyToPdbIdx.push(i);
            }
            mol3D = tryParse({sanitize: true, removeHs: true});
          }
          if (mol3D) mol3DFromPdb = true;
        }
      }
    }
    if (!mol3D) return null;

    // Helper: remap match indices back to PDB-serial-compatible indices.
    // If we removed H, map back; otherwise indices are already correct.
    const remapMatch = (match: number[]): number[] => {
      if (!heavyToPdbIdx) return match;
      return match.map((idx) =>
        idx >= 0 && idx < heavyToPdbIdx!.length ? heavyToPdbIdx![idx] : -1);
    };

    // Build PDB serial lookup: pdbSerials[molblockIdx] = actual PDB serial.
    // Reuse pdbAtomsParsed from above (avoid re-parsing).
    const pdbSerials: number[] | undefined = (mol3DFromPdb && pdbAtomsParsed) ?
      pdbAtomsParsed.map((a) => a.serial) : undefined;

    // Tier 1: strict substructure match (works for SDF/molblock 3D files
    // where bond orders are preserved).
    const strictMatch = trySubstructMatch(mol2D, mol3D);
    if (strictMatch) {
      return {
        mapping: remapMatch(strictMatch),
        method: 'substruct',
        mappedCount: strictMatch.filter((x) => x >= 0).length,
        pdbSerials,
      };
    }

    // Tier 2a: SMARTS-based bond-order-agnostic match.
    // Generate SMARTS from 2D mol, replace all bond operators with ~ (any),
    // then match against the 3D mol. This preserves atom types, H counts,
    // and ring info while ignoring bond orders — avoids the false symmetry
    // that flattening all bonds to single creates (e.g., C=O vs C-OH).
    {
      let qmol: RDMol | null = null;
      try {
        const smarts2D = mol2D.get_smarts();
        // Replace explicit bond symbols with ~ (any bond).
        // SMARTS bond chars: - (single), = (double), # (triple), : (aromatic)
        // Also handle @ in ring bonds — only replace standalone bond chars.
        const smartsAnyBond = smarts2D.replace(/((?:\]|[A-Z]|[a-z]|\*))([=\-#:])((?:\[|[A-Z]|[a-z]|\*))/g, '$1~$3');
        qmol = rdkit.get_qmol(smartsAnyBond);
        if (qmol?.is_valid()) {
          const smartsMatch = trySubstructMatch(qmol, mol3D);
          if (smartsMatch) {
            return {
              mapping: remapMatch(smartsMatch),
              method: 'substruct',
              mappedCount: smartsMatch.filter((x) => x >= 0).length,
              pdbSerials,
            };
          }
        }
      } catch {/* SMARTS matching not supported — fall through */} finally {qmol?.delete();}
    }

    // Tier 2b: Flattened molblock fallback.
    // Flatten all bond orders to 1 (single) and match. Less precise than
    // SMARTS but works when SMARTS generation fails.
    {
      let mol2DFlat: RDMol | null = null;
      try {
        const mol2DBlock = mol2D.get_molblock();
        const mol2DFlatBlock = flattenBondOrders(mol2DBlock);
        mol2DFlat = rdkit.get_mol(mol2DFlatBlock, JSON.stringify({sanitize: false, removeHs: false}));
        if (mol2DFlat?.is_valid()) {
          const flatMatch = trySubstructMatch(mol2DFlat, mol3D);
          if (flatMatch) {
            return {
              mapping: remapMatch(flatMatch),
              method: 'substruct',
              mappedCount: flatMatch.filter((x) => x >= 0).length,
              pdbSerials,
            };
          }
        }
      } finally {mol2DFlat?.delete();}
    }

    // Tier 3: heavy-atom serial-order fallback.
    const numHeavy2D = mol2D.get_num_atoms();
    if (mol3DFromPdb && pdbAtomsParsed) {
      const heavyIndices: number[] = [];
      for (let i = 0; i < pdbAtomsParsed.length; i++) {
        if (pdbAtomsParsed[i].element !== 'H' && pdbAtomsParsed[i].element !== 'D')
          heavyIndices.push(i);
      }
      const n = Math.min(numHeavy2D, heavyIndices.length);
      const fallback: number[] = [];
      for (let i = 0; i < numHeavy2D; i++)
        fallback.push(i < n ? heavyIndices[i] : -1);
      return {mapping: fallback, method: 'heavy-atom-order', mappedCount: n, pdbSerials};
    }

    const numHeavy3D = mol3D.get_num_atoms();
    const n = Math.min(numHeavy2D, numHeavy3D);
    const fallback: number[] = [];
    for (let i = 0; i < numHeavy2D; i++)
      fallback.push(i < n ? i : -1);
    return {mapping: fallback, method: 'heavy-atom-order', mappedCount: n};
  } catch (err) {
    return null;
  } finally {
    mol2D?.delete();
    mol3D?.delete();
  }
}

/**
 * Runs a strict substructure match: mol2D as query against mol3D.
 * Returns the match array or null if no match.
 */
function trySubstructMatch(mol2D: RDMol, mol3D: RDMol): number[] | null {
  try {
    const matchJson = mol3D.get_substruct_match(mol2D);
    if (!matchJson) return null;
    const match = JSON.parse(matchJson);

    // RDKit WASM returns different formats depending on version/options:
    // Format 1: plain array [3, 5, 2, ...] — match[i] = 3D idx for query atom i
    if (Array.isArray(match) && match.length > 0)
      return match;

    if (typeof match === 'object' && !Array.isArray(match)) {
      // Format 2: {"atoms": [3, 5, 2, ...], "bonds": [...]}
      // The "atoms" array is match[i] = 3D idx for query atom i.
      if (Array.isArray(match.atoms) && match.atoms.length > 0)
        return match.atoms;

      // Format 3: {3DIdx: 2DIdx, ...} — invert to get 2D→3D mapping
      const numAtoms2D = mol2D.get_num_atoms();
      const result: number[] = new Array(numAtoms2D).fill(-1);
      for (const [key3D, val2D] of Object.entries(match)) {
        const i2D = typeof val2D === 'number' ? val2D : parseInt(val2D as string, 10);
        const i3D = parseInt(key3D, 10);
        if (i2D >= 0 && i2D < numAtoms2D) result[i2D] = i3D;
      }
      return result.some((x) => x >= 0) ? result : null;
    }
    return null;
  } catch {
    return null;
  }
}

/** Replaces all bond types in a V2000 molblock with single bonds (type 1).
 *  This makes the molecule "bond-order-agnostic" for substructure matching
 *  against a 3D structure that has no bond-order information (e.g. PDB). */
export function flattenBondOrders(molblock: string): string {
  const lines = molblock.split('\n');
  // The counts line is at index 3. Parse atom count to know where bonds start.
  if (lines.length < 5) return molblock;
  const countsLine = lines[3];
  const nAtoms = parseInt(countsLine.substring(0, 3).trim(), 10);
  const nBonds = parseInt(countsLine.substring(3, 6).trim(), 10);
  if (isNaN(nAtoms) || isNaN(nBonds)) return molblock;

  // Bond block starts at line 4 + nAtoms.
  const bondStart = 4 + nAtoms;
  for (let i = bondStart; i < bondStart + nBonds && i < lines.length; i++) {
    const line = lines[i];
    if (line.length < 9) continue;
    // Bond type is at positions 6-8 (3 chars). Replace with "  1".
    lines[i] = line.substring(0, 6) + '  1' + line.substring(9);
  }
  return lines.join('\n');
}
