/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {before, category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';
import {getRdKitModule} from '@datagrok-libraries/bio/src/chem/rdkit-module';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';

import {
  assembleMolecule,
  ChemEnumModes,
  countForCore,
  enumerate,
  extractRNumbers,
  makeCore,
  makeRGroup,
  moveStartRLabelToBranch,
  normalizeRLabels,
  pickFreeRingDigits,
  remapSingleRLabel,
  substituteRLabelWithRingDigit,
  validateParams,
} from '../polytool/pt-chem-enum';

import {_package} from '../package-test';

/** Parses canonical SMILES via RDKit so two syntactically different SMILES can be compared structurally. */
function canon(smi: string, rdkit: RDModule): string {
  const mol = rdkit.get_mol(smi);
  try {
    expect(mol.is_valid(), true, `RDKit rejected SMILES: ${smi}`);
    return mol.get_smiles();
  } finally { mol.delete(); }
}

// ─── Regex / normalization / remap — no RDKit required ──────────────────────

category('PolyTool: ChemEnum: R-labels', () => {
  test('extract from every supported form', async () => {
    expectArray(extractRNumbers('C[1*]'), [1]);
    expectArray(extractRNumbers('C[*:1]'), [1]);
    expectArray(extractRNumbers('C[*1]'), [1]);
    expectArray(extractRNumbers('C[R1]'), [1]);
    expectArray(extractRNumbers('C[R:1]'), [1]);
  });

  test('extract multi-digit numbers', async () => {
    expectArray(extractRNumbers('C1CC(C[*:100])(N[*:1])C([*:101])N1[*:2]'), [1, 2, 100, 101]);
    expectArray(extractRNumbers('C[R100]N[R:101]'), [100, 101]);
  });

  test('normalize all spellings to [*:N]', async () => {
    expect(normalizeRLabels('C[1*]'), 'C[*:1]');
    expect(normalizeRLabels('C[*:1]'), 'C[*:1]');
    expect(normalizeRLabels('C[*1]'), 'C[*:1]');
    expect(normalizeRLabels('C[R1]'), 'C[*:1]');
    expect(normalizeRLabels('C[R:1]'), 'C[*:1]');
    expect(normalizeRLabels('C[R100]N[R:101]'), 'C[*:100]N[*:101]');
  });

  test('remap the one-and-only R-label to a new slot', async () => {
    expect(remapSingleRLabel('C[1*]', 5), 'C[*:5]');
    expect(remapSingleRLabel('C[R:3]CC', 1), 'C[*:1]CC');
    expect(remapSingleRLabel('[*:7]C', 2), '[*:2]C');
  });

  test('deduplicate repeated R numbers', async () => {
    expectArray(extractRNumbers('C[*:1]C[*:1]'), [1]);
  });

  test('normalize does not touch regular bracket atoms', async () => {
    expect(normalizeRLabels('C[NH3+]C[*:1]'), 'C[NH3+]C[*:1]');
    expect(normalizeRLabels('[13C]C[*:2]'), '[13C]C[*:2]');
  });
});

// ─── SMILES surgery helpers — no RDKit required ─────────────────────────────

category('PolyTool: ChemEnum: SMILES surgery', () => {
  test('moveStartRLabelToBranch: bare atom at start', async () => {
    expect(moveStartRLabelToBranch('[*:1]CC'), 'C([*:1])C');
    expect(moveStartRLabelToBranch('[*:1]Cl'), 'Cl([*:1])');
    expect(moveStartRLabelToBranch('[*:1][NH3+]C'), '[NH3+]([*:1])C');
  });

  test('moveStartRLabelToBranch: with explicit bond', async () => {
    expect(moveStartRLabelToBranch('[*:1]=CC'), 'C(=[*:1])C');
    expect(moveStartRLabelToBranch('[*:1]#C'), 'C(#[*:1])');
  });

  test('moveStartRLabelToBranch: leaves interior unchanged', async () => {
    expect(moveStartRLabelToBranch('CC[*:1]'), 'CC[*:1]');
    expect(moveStartRLabelToBranch('CC([*:1])CC'), 'CC([*:1])CC');
  });

  test('substituteRLabelWithRingDigit collapses branch form', async () => {
    expect(substituteRLabelWithRingDigit('CC([*:1])CC', 1, '9'), 'CC9CC');
    expect(substituteRLabelWithRingDigit('CC([*:2])([*:1])CC', 1, '%50'), 'CC([*:2])%50CC');
  });

  test('substituteRLabelWithRingDigit handles trailing form', async () => {
    expect(substituteRLabelWithRingDigit('CC[*:1]', 1, '7'), 'CC7');
    expect(substituteRLabelWithRingDigit('CC=[*:2]', 2, '%42'), 'CC=%42');
  });

  test('pickFreeRingDigits avoids digits already in use', async () => {
    // 1 is in use in the core → pick 2, 3
    const d = pickFreeRingDigits(['C1CCCCC1[*:1]', '[*:1]C'], 2);
    expect(d.includes(1), false);
    expect(d.length, 2);
  });

  test('pickFreeRingDigits ignores digits inside bracket atoms', async () => {
    // `[*:1]` contains "1" but it's inside brackets — must not count as ring digit 1
    const d = pickFreeRingDigits(['C[*:1]', '[*:1]C'], 1);
    expectArray(d, [1]);
  });
});

// ─── Construction validators — no RDKit required ────────────────────────────

category('PolyTool: ChemEnum: validators', () => {
  test('core without any R-label is rejected', async () => {
    const c = makeCore('CCO', 'row 1');
    expect(c.error != null, true);
  });

  test('r-group with zero R-labels is rejected', async () => {
    const rg = makeRGroup('CC', 1, 'row 1');
    expect(rg.error != null, true);
  });

  test('r-group with multiple R-labels is rejected', async () => {
    const rg = makeRGroup('[*:1]CC[*:2]', 1, 'row 1');
    expect(rg.error != null, true);
  });

  test('r-group with wrong R-number is auto-remapped', async () => {
    const rg = makeRGroup('CC[*:3]', 1, 'row 1');
    expect(rg.error == null, true);
    expect(rg.rNumber, 1);
    expect(rg.sourceRNumber, 3);
    expect(rg.smiles, 'CC[*:1]');
  });
});

// ─── Validation + count prediction — no RDKit required ──────────────────────

category('PolyTool: ChemEnum: count & validate', () => {
  test('zip length mismatch is flagged', async () => {
    const c = makeCore('C[*:1]C[*:2]', 'c1');
    const v = validateParams({
      cores: [c],
      rGroups: new Map([
        [1, [makeRGroup('C[*:1]', 1, 'r1-a'), makeRGroup('N[*:1]', 1, 'r1-b')]],
        [2, [makeRGroup('O[*:2]', 2, 'r2-a')]],
      ]),
      mode: ChemEnumModes.Zip,
    });
    expect(v.ok, false);
    expect(v.errors.some((e) => e.includes('Zip')), true);
  });

  test('cartesian count = product across Rs, summed over cores', async () => {
    const cA = makeCore('C[*:1]', 'a'); // R1 only
    const cB = makeCore('C[*:1]N[*:2]', 'b'); // R1 + R2
    const v = validateParams({
      cores: [cA, cB],
      rGroups: new Map([
        [1, [makeRGroup('C[*:1]', 1, 'a'), makeRGroup('N[*:1]', 1, 'b'), makeRGroup('O[*:1]', 1, 'c')]], // 3
        [2, [makeRGroup('C[*:2]', 2, 'a'), makeRGroup('N[*:2]', 2, 'b')]], // 2
      ]),
      mode: ChemEnumModes.Cartesian,
    });
    // cA contributes 3 (just R1), cB contributes 3*2 = 6 → total 9
    expect(v.predictedCount, 9);
    expect(v.ok, true);
  });

  test('zip count = N per core, summed', async () => {
    const cA = makeCore('C[*:1]', 'a');
    const cB = makeCore('C[*:1]N[*:2]', 'b');
    const v = validateParams({
      cores: [cA, cB],
      rGroups: new Map([
        [1, [makeRGroup('C[*:1]', 1, 'a'), makeRGroup('N[*:1]', 1, 'b')]], // 2
        [2, [makeRGroup('O[*:2]', 2, 'a'), makeRGroup('S[*:2]', 2, 'b')]], // 2
      ]),
      mode: ChemEnumModes.Zip,
    });
    // cA uses only R1 → zip length = 2; cB uses R1 & R2 (both 2) → zip length = 2 → total 4
    expect(v.predictedCount, 4);
  });

  test('uncovered R-number is flagged per core', async () => {
    const c = makeCore('C[*:1]N[*:7]', 'c');
    const v = validateParams({
      cores: [c],
      rGroups: new Map([[1, [makeRGroup('O[*:1]', 1, 'a')]]]),
      mode: ChemEnumModes.Cartesian,
    });
    expect(v.ok, false);
    expect(v.errors.some((e) => e.includes('R7')), true);
  });

  test('hard cap rejects over-limit enumerations', async () => {
    const c = makeCore('C[*:1]N[*:2]', 'c');
    const many1 = Array.from({length: 1000}, (_, i) => makeRGroup('C[*:1]', 1, `r1-${i}`));
    const many2 = Array.from({length: 1001}, (_, i) => makeRGroup('O[*:2]', 2, `r2-${i}`));
    const v = validateParams({
      cores: [c],
      rGroups: new Map([[1, many1], [2, many2]]),
      mode: ChemEnumModes.Cartesian,
      maxResults: 1_000_000,
    });
    expect(v.overCap, true);
    expect(v.ok, false);
  });

  test('countForCore: zip with mismatched per-core lengths → -1', async () => {
    const c = makeCore('C[*:1]N[*:2]', 'c');
    const {count} = countForCore(c, new Map([
      [1, [makeRGroup('C[*:1]', 1, 'a'), makeRGroup('N[*:1]', 1, 'b')]],
      [2, [makeRGroup('O[*:2]', 2, 'a')]],
    ]), ChemEnumModes.Zip);
    expect(count, -1);
  });
});

// ─── Assembly + full enumeration — needs RDKit ──────────────────────────────

category('PolyTool: ChemEnum: assembly', () => {
  let rdkit: RDModule;

  before(async () => { rdkit = await getRdKitModule(); });

  test('single R-group: joins correctly', async () => {
    const res = assembleMolecule('C[*:1]', new Map([[1, 'O[*:1]']]), rdkit);
    expect(res != null, true);
    expect(res, canon('CO', rdkit));
  });

  test('multi R-group: R1 and R2 both joined', async () => {
    const res = assembleMolecule('C[*:1]N[*:2]', new Map([[1, 'O[*:1]'], [2, 'S[*:2]']]), rdkit);
    expect(res != null, true);
    expect(res, canon('OCNS', rdkit));
  });

  test('R-label at SMILES start is handled', async () => {
    const res = assembleMolecule('[*:1]CC', new Map([[1, 'Br[*:1]']]), rdkit);
    expect(res != null, true);
    expect(res, canon('BrCC', rdkit));
  });

  test('core uses existing ring digit 1 — free digit is picked', async () => {
    // Cyclohexyl with an R1 sticking off; joining must not clash with ring-1
    const res = assembleMolecule('C1CCCCC1[*:1]', new Map([[1, 'O[*:1]']]), rdkit);
    expect(res != null, true);
    expect(res, canon('OC1CCCCC1', rdkit));
  });

  test('multi-digit R numbers (R100, R101)', async () => {
    const res = assembleMolecule(
      'C1CC(C[*:100])(N[*:1])C([*:101])N1[*:2]',
      new Map([[1, 'C[*:1]'], [2, 'C[*:2]'], [100, 'C[*:100]'], [101, 'C[*:101]']]),
      rdkit);
    expect(res != null, true);
    // Sanity check: result parses and has no remaining R labels
    const mol = rdkit.get_mol(res!);
    try {
      expect(mol.is_valid(), true);
      expect(res!.includes('[*:'), false);
    } finally { mol.delete(); }
  });
});

// ─── End-to-end enumerate ────────────────────────────────────────────────────

category('PolyTool: ChemEnum: enumerate', () => {
  let rdkit: RDModule;

  before(async () => { rdkit = await getRdKitModule(); });

  test('cartesian across two cores with different R-sets', async () => {
    const cA = makeCore('C[*:1]', 'A');
    const cB = makeCore('C[*:1]N[*:2]', 'B');
    const r1 = [makeRGroup('O[*:1]', 1, 'r1-a'), makeRGroup('S[*:1]', 1, 'r1-b')];
    const r2 = [makeRGroup('C[*:2]', 2, 'r2-a'), makeRGroup('N[*:2]', 2, 'r2-b'), makeRGroup('F[*:2]', 2, 'r2-c')];
    const results = enumerate(
      {cores: [cA, cB], rGroups: new Map([[1, r1], [2, r2]]), mode: ChemEnumModes.Cartesian},
      rdkit)!;
    // cA: |R1|=2 → 2; cB: |R1|*|R2|=6 → total 8
    expect(results.length, 8);
    const first = results.filter((r) => r.coreSmiles === 'C[*:1]');
    const second = results.filter((r) => r.coreSmiles === 'C[*:1]N[*:2]');
    expect(first.length, 2);
    expect(second.length, 6);
  });

  test('zip across two cores with shared R-list length', async () => {
    const cA = makeCore('C[*:1]', 'A');
    const cB = makeCore('C[*:1]N[*:2]', 'B');
    const r1 = [makeRGroup('O[*:1]', 1, 'r1-a'), makeRGroup('S[*:1]', 1, 'r1-b')];
    const r2 = [makeRGroup('C[*:2]', 2, 'r2-a'), makeRGroup('N[*:2]', 2, 'r2-b')];
    const results = enumerate(
      {cores: [cA, cB], rGroups: new Map([[1, r1], [2, r2]]), mode: ChemEnumModes.Zip},
      rdkit)!;
    // cA zip-length = 2, cB zip-length = 2 → total 4
    expect(results.length, 4);
    // each cB result has O+C (i=0) or S+N (i=1)
    const cBResults = results.filter((r) => r.coreSmiles === 'C[*:1]N[*:2]');
    expect(cBResults.length, 2);
    expect(cBResults[0].rGroupSmilesByNum.get(1), 'O[*:1]');
    expect(cBResults[0].rGroupSmilesByNum.get(2), 'C[*:2]');
    expect(cBResults[1].rGroupSmilesByNum.get(1), 'S[*:1]');
    expect(cBResults[1].rGroupSmilesByNum.get(2), 'N[*:2]');
  });

  test('invalid params return null from enumerate()', async () => {
    const c = makeCore('C[*:1]N[*:2]', 'c');
    const res = enumerate({
      cores: [c],
      rGroups: new Map([[1, [makeRGroup('O[*:1]', 1, 'x')]]]), // R2 missing
      mode: ChemEnumModes.Cartesian,
    }, rdkit);
    expect(res, null);
  });
});

// ─── File-based smoke test using the two test CSVs ──────────────────────────

category('PolyTool: ChemEnum: CSV fixtures', () => {
  let rdkit: RDModule;

  before(async () => { rdkit = await getRdKitModule(); });

  async function loadCsv(relPath: string): Promise<DG.DataFrame> {
    const txt = await _package.files.readAsText(relPath);
    return DG.DataFrame.fromCsv(txt);
  }

  test('cores file parses and all rows validate', async () => {
    const df = await loadCsv('tests/chem_enum_cores.csv');
    const col = df.col('Core')!;
    const cores = Array.from({length: col.length}, (_, i) =>
      makeCore(col.get(i), `Cores[${i}]`, rdkit));
    expect(cores.length, 4);
    expect(cores.every((c) => !c.error), true, cores.map((c) => c.error ?? '').filter((e) => e).join('; '));
    // All four reference R1 and R2
    for (const c of cores)
      expect(c.rNumbers.includes(1) && c.rNumbers.includes(2), true, `core ${c.id} R-numbers: ${c.rNumbers}`);
  });

  test('r-groups file: columns R1, R2, R3 populate their slots', async () => {
    const df = await loadCsv('tests/chem_enum_rgroups.csv');
    const buildList = (colName: string, n: number) => {
      const col = df.col(colName)!;
      const out = [];
      for (let i = 0; i < col.length; i++) {
        const v = col.get(i);
        if (v == null || String(v).trim() === '') continue;
        out.push(makeRGroup(String(v), n, `${colName}[${i}]`, rdkit));
      }
      return out;
    };
    const r1 = buildList('R1', 1);
    const r2 = buildList('R2', 2);
    const r3 = buildList('R3', 3);
    expect(r1.length, 4);
    expect(r2.length, 4);
    expect(r3.length, 1);
    for (const list of [r1, r2, r3])
      expect(list.every((rg) => !rg.error), true, list.map((rg) => rg.error ?? '').filter((e) => e).join('; '));
  });

  test('end-to-end cartesian with CSV fixtures produces valid, canonical SMILES', async () => {
    const coresDf = await loadCsv('tests/chem_enum_cores.csv');
    const rgsDf = await loadCsv('tests/chem_enum_rgroups.csv');

    const coreCol = coresDf.col('Core')!;
    const cores = Array.from({length: coreCol.length}, (_, i) =>
      makeCore(coreCol.get(i), `Cores[${i}]`, rdkit));

    const rGroups = new Map<number, any[]>();
    for (const colName of ['R1', 'R2', 'R3']) {
      const n = parseInt(colName.substring(1), 10);
      const col = rgsDf.col(colName)!;
      const list = [];
      for (let i = 0; i < col.length; i++) {
        const v = col.get(i);
        if (v == null || String(v).trim() === '') continue;
        list.push(makeRGroup(String(v), n, `${colName}[${i}]`, rdkit));
      }
      rGroups.set(n, list);
    }

    const results = enumerate({cores, rGroups, mode: ChemEnumModes.Cartesian}, rdkit)!;
    expect(results != null && results.length > 0, true);
    // All outputs must parse and contain no residual R-labels
    for (const r of results) {
      expect(r.smiles.includes('[*:'), false, `residual R-label in result: ${r.smiles}`);
      const m = rdkit.get_mol(r.smiles);
      try {
        expect(m.is_valid(), true, `invalid SMILES: ${r.smiles}`);
      } finally { m.delete(); }
    }
  });
});
