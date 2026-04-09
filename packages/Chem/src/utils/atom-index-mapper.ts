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
 * Three-tier fallback chain:
 *   1. Strict substructure match (handles BRANCH reorder + H addition)
 *   2. Relaxed match (sanitize, ignore chirality — handles protonation)
 *   3. Heavy-atom serial order (H-skip — last resort when topology
 *      matching fails entirely)
 */

import {RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';

// -- PDB/PDBQT → Molblock converter -----------------------------------------

/** Typical single-bond max distance per element pair (Å). */
const BOND_DIST_MAX = 1.9; // covers C-C (1.54), C-N (1.47), C-O (1.43), C-S (1.82)
const BOND_DIST_MIN = 0.5;

interface PdbAtom {
  serial: number;
  element: string;
  x: number; y: number; z: number;
}

/** Parses ATOM/HETATM lines from PDB or PDBQT text. */
function parsePdbAtoms(pdbText: string): PdbAtom[] {
  const atoms: PdbAtom[] = [];
  for (const line of pdbText.split('\n')) {
    const rec = line.substring(0, 6).trim();
    if (rec !== 'ATOM' && rec !== 'HETATM') continue;
    // PDB format: cols 31-38 x, 39-46 y, 47-54 z, 77-78 element
    const x = parseFloat(line.substring(30, 38));
    const y = parseFloat(line.substring(38, 46));
    const z = parseFloat(line.substring(46, 54));
    // Element symbol: columns 77-78 in standard PDB. For PDBQT, might
    // be in columns 77-78 or can be derived from atom name (cols 13-16).
    let element = line.length >= 78 ? line.substring(76, 78).trim() : '';
    if (!element) {
      // Derive from atom name (cols 13-16): strip digits, take first letter(s).
      const atomName = line.substring(12, 16).trim();
      element = atomName.replace(/[0-9]/g, '').substring(0, 2).trim();
      if (element.length === 2)
        element = element[0].toUpperCase() + element[1].toLowerCase();
      else
        element = element.toUpperCase();
    }
    if (isNaN(x) || isNaN(y) || isNaN(z) || !element) continue;
    atoms.push({serial: atoms.length + 1, element, x, y, z});
  }
  return atoms;
}

/** Infers bonds from interatomic distances. Returns [i, j] pairs (0-based). */
function inferBonds(atoms: PdbAtom[]): [number, number][] {
  const bonds: [number, number][] = [];
  for (let i = 0; i < atoms.length; i++) {
    for (let j = i + 1; j < atoms.length; j++) {
      const dx = atoms[i].x - atoms[j].x;
      const dy = atoms[i].y - atoms[j].y;
      const dz = atoms[i].z - atoms[j].z;
      const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
      if (dist >= BOND_DIST_MIN && dist <= BOND_DIST_MAX)
        bonds.push([i, j]);
    }
  }
  return bonds;
}

/** Converts PDB/PDBQT text to a V2000 molblock that RDKit can parse.
 *  Returns null if the text doesn't contain valid ATOM/HETATM records. */
function pdbToMolblock(pdbText: string): string | null {
  const atoms = parsePdbAtoms(pdbText);
  if (atoms.length === 0) return null;

  const bonds = inferBonds(atoms);

  // V2000 header
  const lines: string[] = [
    '', // molecule name
    '     RDKit          3D', // program/timestamp
    '', // comment
  ];

  // Counts line: aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
  const nAtoms = atoms.length.toString().padStart(3);
  const nBonds = bonds.length.toString().padStart(3);
  lines.push(`${nAtoms}${nBonds}  0  0  0  0  0  0  0  0999 V2000`);

  // Atom block
  for (const a of atoms) {
    const xs = a.x.toFixed(4).padStart(10);
    const ys = a.y.toFixed(4).padStart(10);
    const zs = a.z.toFixed(4).padStart(10);
    const el = a.element.padEnd(3);
    lines.push(`${xs}${ys}${zs} ${el} 0  0  0  0  0  0  0  0  0  0  0  0`);
  }

  // Bond block
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
  method: 'substruct' | 'relaxed' | 'heavy-atom-order';
  /** Number of successfully mapped atoms. */
  mappedCount: number;
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

  try {
    mol2D = rdkit.get_mol(smiles2D);
    if (!mol2D || !mol2D.is_valid()) return null;

    // Try parsing the 3D string directly first (works for SDF/molblock).
    mol3D = rdkit.get_mol(mol3DStr);

    // If direct parsing fails, try converting PDB/PDBQT to molblock.
    if (!mol3D || !mol3D.is_valid()) {
      mol3D?.delete();
      mol3D = null;
      const molblock = pdbToMolblock(mol3DStr);
      if (molblock) {
        mol3D = rdkit.get_mol(molblock, JSON.stringify({sanitize: true, removeHs: false}));
        if (mol3D && !mol3D.is_valid()) {mol3D.delete(); mol3D = null;}
        // Retry with removeHs if sanitization with H fails
        if (!mol3D) {
          mol3D = rdkit.get_mol(molblock, JSON.stringify({sanitize: true, removeHs: true}));
          if (mol3D && !mol3D.is_valid()) {mol3D.delete(); mol3D = null;}
        }
      }
    }
    if (!mol3D) return null;

    // Tier 1: strict substructure match.
    // Use 2D as query, 3D as target. The match result is an array where
    // match[i] = the index in mol3D that corresponds to atom i in mol2D.
    const strictMatch = trySubstructMatch(mol2D, mol3D);
    if (strictMatch) {
      return {
        mapping: strictMatch,
        method: 'substruct',
        mappedCount: strictMatch.filter((x) => x >= 0).length,
      };
    }

    // Tier 2: relaxed match — strip stereo, sanitize, retry.
    // Handles protonation-state differences and minor canonicalization
    // discrepancies.
    const relaxedMatch = tryRelaxedMatch(rdkit, smiles2D, mol3DStr);
    if (relaxedMatch) {
      return {
        mapping: relaxedMatch,
        method: 'relaxed',
        mappedCount: relaxedMatch.filter((x) => x >= 0).length,
      };
    }

    // Tier 3: heavy-atom serial-order fallback.
    // Assumes heavy atoms appear in the same order in 2D and 3D.
    // This is wrong for PDBQT (branch-reordered) but better than nothing
    // for simple PDB/SDF where atom order IS preserved.
    const numHeavy2D = mol2D.get_num_atoms(true);
    const numHeavy3D = mol3D.get_num_atoms(true);
    const n = Math.min(numHeavy2D, numHeavy3D);
    const fallback: number[] = [];
    for (let i = 0; i < numHeavy2D; i++)
      fallback.push(i < n ? i : -1);
    return {
      mapping: fallback,
      method: 'heavy-atom-order',
      mappedCount: n,
    };
  } catch {
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
    const match: number[] = JSON.parse(matchJson);
    // get_substruct_match returns {atomIdx3D: atomIdx2D, ...} or an array
    // depending on RDKit version. Handle both formats.
    if (Array.isArray(match) && match.length > 0) {
      // match[i] = index in mol3D for query atom i
      return match;
    }
    // If it's an object, convert: keys = 3D indices, values = 2D indices
    if (typeof match === 'object' && !Array.isArray(match)) {
      const numAtoms2D = mol2D.get_num_atoms(true);
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

/**
 * Relaxed match: re-parse both molecules with sanitization, strip
 * stereochemistry, and retry substructure matching. Handles cases where
 * protonation changes or minor canonical differences block the strict
 * match.
 */
function tryRelaxedMatch(
  rdkit: RDModule,
  smiles2D: string,
  mol3DStr: string,
): number[] | null {
  let m2: RDMol | null = null;
  let m3: RDMol | null = null;
  try {
    // Re-parse with relaxed options
    m2 = rdkit.get_mol(smiles2D, JSON.stringify({sanitize: true, removeHs: true}));
    if (!m2 || !m2.is_valid()) return null;
    m3 = rdkit.get_mol(mol3DStr, JSON.stringify({sanitize: true, removeHs: true}));
    if (!m3 || !m3.is_valid()) {
      m3?.delete();
      m3 = null;
      // Try PDB/PDBQT → molblock conversion
      const molblock = pdbToMolblock(mol3DStr);
      if (molblock) {
        m3 = rdkit.get_mol(molblock, JSON.stringify({sanitize: true, removeHs: true}));
        if (m3 && !m3.is_valid()) {m3.delete(); m3 = null;}
      }
    }
    if (!m3) return null;

    return trySubstructMatch(m2, m3);
  } catch {
    return null;
  } finally {
    m2?.delete();
    m3?.delete();
  }
}
