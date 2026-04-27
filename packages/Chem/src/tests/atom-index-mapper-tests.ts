import {category, test, expect, before} from '@datagrok-libraries/test/src/test';
import * as grok from 'datagrok-api/grok';

import {
  parsePdbAtoms, inferBonds, pdbAtomsToMolblock, flattenBondOrders,
  mapAtomIndices2Dto3D, PdbAtom,
} from '../utils/atom-index-mapper';

// -- Test data ---------------------------------------------------------------

/** Minimal PDB with 3 atoms: N at serial 5, C at serial 10, O at serial 15. */
const PDB_BASIC = [
  'ATOM      5  N   LIG     1       1.000   2.000   3.000  1.00  0.00           N  ',
  'ATOM     10  CA  LIG     1       2.520   2.000   3.000  1.00  0.00           C  ',
  'ATOM     15  O   LIG     1       3.960   2.000   3.000  1.00  0.00           O  ',
].join('\n');

/** PDBQT with non-standard element types (A=aromatic C, OA=O, NA=N, HD=H).
 *  PDBQT format: element type is at columns 77-78 but may contain PDBQT
 *  types like A, OA, NA, HD. The parser normalizes these. */
const PDBQT_ELEMENTS = [
  'ATOM      1  A1  LIG     1       0.000   0.000   0.000  0.00  0.00    -0.069 A ',
  'ATOM      2  OA1 LIG     1       1.400   0.000   0.000  0.00  0.00    -0.200 OA',
  'ATOM      3  NA1 LIG     1       2.800   0.000   0.000  0.00  0.00    -0.150 NA',
  'ATOM      4  HD1 LIG     1       3.500   0.800   0.000  0.00  0.00     0.100 HD',
].join('\n');

/** Ethanol (CCO) as PDB with atoms reordered: O first, then two C's. */
const PDB_ETHANOL_REORDERED = [
  'ATOM      1  O   LIG     1       2.460   0.000   0.000  1.00  0.00           O  ',
  'ATOM      2  C2  LIG     1       1.230   0.000   0.000  1.00  0.00           C  ',
  'ATOM      3  C1  LIG     1       0.000   0.000   0.000  1.00  0.00           C  ',
].join('\n');

const SMILES_ETHANOL = 'CCO';

// -- Tests -------------------------------------------------------------------

category('atom index mapper', () => {
  let rdkitModule: any;

  before(async () => {
    rdkitModule = await grok.functions.call('Chem:getRdKitModule');
  });

  // -- parsePdbAtoms ---------------------------------------------------------

  test('parsePdbAtoms — basic', async () => {
    const atoms = parsePdbAtoms(PDB_BASIC);
    expect(atoms.length, 3);
    expect(atoms[0].element, 'N');
    expect(atoms[0].serial, 5);
    expect(atoms[1].element, 'C');
    expect(atoms[1].serial, 10);
    expect(atoms[2].element, 'O');
    expect(atoms[2].serial, 15);
    // Coordinates
    expect(Math.abs(atoms[0].x - 1.0) < 0.01, true);
    expect(Math.abs(atoms[0].y - 2.0) < 0.01, true);
    expect(Math.abs(atoms[0].z - 3.0) < 0.01, true);
  });

  test('parsePdbAtoms — PDBQT element normalization', async () => {
    const atoms = parsePdbAtoms(PDBQT_ELEMENTS);
    expect(atoms.length, 4);
    expect(atoms[0].element, 'C'); // A → C
    expect(atoms[1].element, 'O'); // OA → O
    expect(atoms[2].element, 'N'); // NA → N
    expect(atoms[3].element, 'H'); // HD → H
  });

  test('parsePdbAtoms — empty input', async () => {
    expect(parsePdbAtoms('').length, 0);
    expect(parsePdbAtoms('REMARK this is not an atom line').length, 0);
  });

  // -- inferBonds ------------------------------------------------------------

  test('inferBonds — covalent radius', async () => {
    // C at origin, C at 1.52Å (typical C-C), O at 2.96Å (C+O = too far from first C)
    const atoms: PdbAtom[] = [
      {serial: 1, element: 'C', x: 0, y: 0, z: 0},
      {serial: 2, element: 'C', x: 1.52, y: 0, z: 0},
      {serial: 3, element: 'O', x: 2.96, y: 0, z: 0},
      {serial: 4, element: 'C', x: 10, y: 10, z: 10}, // distant — no bonds
    ];
    const bonds = inferBonds(atoms);
    // C(0)–C(1) = 1.52Å ≤ 0.76+0.76+0.45 = 1.97Å ✓
    // C(1)–O(2) = 1.44Å ≤ 0.76+0.66+0.45 = 1.87Å ✓
    // C(0)–O(2) = 2.96Å > 1.87Å ✗
    // C(3) at 10,10,10 — too far from all
    expect(bonds.length, 2);
    expect(bonds[0][0], 0);
    expect(bonds[0][1], 1);
    expect(bonds[1][0], 1);
    expect(bonds[1][1], 2);
  });

  // -- pdbAtomsToMolblock ----------------------------------------------------

  test('pdbAtomsToMolblock — valid V2000', async () => {
    const atoms = parsePdbAtoms(PDB_BASIC);
    const molblock = pdbAtomsToMolblock(atoms);
    expect(molblock !== null, true);
    expect(molblock!.includes('V2000'), true);
    expect(molblock!.includes('M  END'), true);
    // 3 atoms, check counts line
    expect(molblock!.includes('  3'), true);
  });

  test('pdbAtomsToMolblock — empty input', async () => {
    expect(pdbAtomsToMolblock([]) === null, true);
  });

  // -- flattenBondOrders -----------------------------------------------------

  test('flattenBondOrders — replaces double bonds', async () => {
    // Minimal molblock with a double bond (type 2)
    const molblock = [
      '',
      '     RDKit          2D',
      '',
      '  2  1  0  0  0  0  0  0  0  0999 V2000',
      '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
      '    1.5000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0',
      '  1  2  2  0',
      'M  END',
    ].join('\n');
    const flat = flattenBondOrders(molblock);
    // Bond line should now be "  1  2  1  0" (type 1)
    expect(flat.includes('  1  2  1'), true);
    // Should NOT contain type 2 in the bond block
    expect(flat.includes('  1  2  2'), false);
  });

  // -- mapAtomIndices2Dto3D --------------------------------------------------

  test('mapAtomIndices2Dto3D — PDB conversion path', async () => {
    const mapping = mapAtomIndices2Dto3D(rdkitModule, SMILES_ETHANOL, PDB_ETHANOL_REORDERED);
    expect(mapping !== null, true);
    // Should succeed via PDB conversion (tier 2), not tier 3 fallback
    expect(mapping!.method, 'substruct');
    expect(mapping!.mappedCount > 0, true);
    // pdbSerials should be populated for PDB input
    expect(Array.isArray(mapping!.pdbSerials), true);
    expect(mapping!.pdbSerials!.length, 3);
  });

  test('mapAtomIndices2Dto3D — mapping correctness', async () => {
    // SMILES CCO: atom 0=C, 1=C, 2=O (RDKit canonical)
    // PDB: serial 1=O, serial 2=C, serial 3=C (reordered)
    const mapping = mapAtomIndices2Dto3D(rdkitModule, SMILES_ETHANOL, PDB_ETHANOL_REORDERED);
    expect(mapping !== null, true);
    // 2D atom 2 (O) should map to PDB atom 0 (the O in the PDB, serial 1)
    // The exact mapping depends on the substruct match result, but:
    // - mapping[2] should point to the O in the 3D mol
    // - pdbSerials[mapping[2]] should be 1 (the O's serial)
    const oIdx2D = 2; // O in SMILES CCO
    const mapped3DIdx = mapping!.mapping[oIdx2D];
    expect(mapped3DIdx >= 0, true);
    if (mapping!.pdbSerials) {
      const serial = mapping!.pdbSerials[mapped3DIdx];
      expect(serial, 1); // O is at PDB serial 1
    }
  });

  test('mapAtomIndices2Dto3D — handles mismatched molecules', async () => {
    // Two completely different molecules — mapping should still return
    // something (fallback tier) without throwing.
    const result = mapAtomIndices2Dto3D(rdkitModule, 'CCO', 'c1ccccc1');
    // Either null or a valid mapping — just verify no crash.
    expect(result === null || Array.isArray(result.mapping), true);
  });
});
