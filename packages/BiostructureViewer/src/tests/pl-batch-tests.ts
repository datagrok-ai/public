import * as grok from 'datagrok-api/grok';

import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {_package} from '../package-test';
import {detectNonWaterHetatmInstances} from '../utils/prolif-panel';


// Tests for the precondition guard used by the per-row PL batch loop.
// The batch handler in package.ts skips rows where hasNonWaterHetatm()
// returns false — these tests lock in that contract so a future change to
// the heuristic doesn't silently change which rows the batch processes.
category('PLBatch', () => {
  test('hasNonWaterHetatm: ligand HETATM passes', async () => {
    // Real PDB needs both ATOM and a non-water HETATM (the heuristic requires
    // both — otherwise the script has nothing meaningful to compute).
    const pdb =
      'ATOM      1  N   ALA A   1      0.000   0.000   0.000  1.00  0.00           N  \n' +
      'HETATM    1  C1  LIG A 100      0.000   0.000   0.000  1.00  0.00           C  ';
    const r = await grok.functions.call(
      `${_package.name}:hasNonWaterHetatm`, {molecule: pdb}) as boolean;
    expect(r, true);
  });

  test('hasNonWaterHetatm: water-only HETATM is rejected', async () => {
    const pdb =
      'ATOM      1  N   ALA A   1      0.000   0.000   0.000  1.00  0.00           N  \n' +
      'HETATM    1  O   HOH A 200      0.000   0.000   0.000  1.00  0.00           O  ';
    const r = await grok.functions.call(
      `${_package.name}:hasNonWaterHetatm`, {molecule: pdb}) as boolean;
    expect(r, false);
  });

  test('hasNonWaterHetatm: empty input is rejected', async () => {
    const r = await grok.functions.call(
      `${_package.name}:hasNonWaterHetatm`, {molecule: ''}) as boolean;
    expect(r, false);
  });

  test('hasNonWaterHetatm: AutoDock-pose marker is rejected (handled by Docking panel)', async () => {
    // The `binding energy` substring check fires before any other guard,
    // so this returns false even if the rest of the input would otherwise
    // qualify. ATOM/HETATM lines are present here only to make the input
    // realistic — they aren't what the assertion exercises.
    const pose =
      'ATOM      1  N   ALA A   1      0.000   0.000   0.000  1.00  0.00           N  \n' +
      'REMARK  Name = test\n' +
      'REMARK  binding energy = -5.0\n' +
      'HETATM    1  C1  LIG A 100      0.000   0.000   0.000  1.00  0.00           C  ';
    const r = await grok.functions.call(
      `${_package.name}:hasNonWaterHetatm`, {molecule: pose}) as boolean;
    expect(r, false);
  });

  // Drives the ligand-picker dropdown in the panel widget. Two contracts
  // matter: (1) common cofactors like HEM are filtered (otherwise they
  // outvote a real inhibitor on the most-common-HETATM heuristic), and
  // (2) the same resname at multiple chain/resid sites yields multiple
  // entries (so the user can pick a specific instance).
  test('detectNonWaterHetatmInstances: filters cofactors, lists ligand instances', async () => {
    // PDB columns: 13–16 atom, 17–20 resname, 21 chain, 22–26 resid (1-indexed).
    // HEM is in PROLIF_SKIP_RESNAMES; STI is a typical small-inhibitor resname.
    const pdb =
      'ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N  \n' +
      'HETATM    1  FE  HEM A 500       0.000   0.000   0.000  1.00  0.00          FE  \n' +
      'HETATM    2  C1  HEM A 500       1.000   0.000   0.000  1.00  0.00           C  \n' +
      'HETATM    3  C1  STI A 401       2.000   0.000   0.000  1.00  0.00           C  \n' +
      'HETATM    4  C1  STI B 402       3.000   0.000   0.000  1.00  0.00           C  ';
    const result = detectNonWaterHetatmInstances(pdb);
    // HEM excluded; both STI instances surface as separate picks.
    expect(result.length, 2);
    expect(result.includes('STI A 401'), true);
    expect(result.includes('STI B 402'), true);
    expect(result.some((s) => s.startsWith('HEM')), false);
  });
});
