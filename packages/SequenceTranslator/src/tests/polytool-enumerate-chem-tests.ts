/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {before, category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';
import {getRdKitModule} from '@datagrok-libraries/bio/src/chem/rdkit-module';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';

import {
  addRGroupsFromSmiles,
  assembleMolecule,
  BUILTIN_R_GROUP_TEMPLATES,
  buildExportColumns,
  ChemEnumModes,
  copyRGroupList,
  countForCore,
  enumerate,
  extractRNumbers,
  invalidTemplateSmiles,
  makeCore,
  makeRGroup,
  moveStartRLabelToBranch,
  normalizeRLabels,
  parseRGroupTemplates,
  pickDefaultTargetR,
  pickFreeRingDigits,
  remapSingleRLabel,
  rGroupTargetWarnings,
  substituteRLabelWithRingDigit,
  uniqueKeepMask,
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

/** Materializes a BitSet keep-mask into a plain boolean[] for element-wise comparison. */
function maskBools(bs: DG.BitSet): boolean[] {
  return Array.from({length: bs.length}, (_, i) => bs.get(i));
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

// ─── Output dedup keep-mask — no RDKit required ─────────────────────────────

category('PolyTool: ChemEnum: dedup', () => {
  test('keeps first occurrence, drops later duplicates', async () => {
    expectArray(maskBools(uniqueKeepMask(['A', 'B', 'A', 'C', 'B'])), [true, true, false, true, false]);
  });

  test('blank and nullish entries are never collapsed', async () => {
    // Empty/invalid rows (e.g. failed canonicalization) must each be kept, not merged into one.
    expectArray(maskBools(uniqueKeepMask(['', '', 'A', null, 'A', undefined])), [true, true, true, true, false, true]);
  });
});

// ─── Dedup actually shrinks the enumerated output (the "Remove duplicates" option) ──

/** Enumerates, canonicalizes every product, and reports raw vs unique counts — mirrors executeEnumeration. */
function enumeratedVsUnique(
  cores: ReturnType<typeof makeCore>[], rGroups: Map<number, ReturnType<typeof makeRGroup>[]>, rdkit: RDModule,
): {raw: number, unique: number} {
  const results = enumerate({cores, rGroups, mode: ChemEnumModes.Cartesian}, rdkit)!;
  const canonical = results.map((r) => canon(r.smiles, rdkit));
  const unique = uniqueKeepMask(canonical).trueCount;
  return {raw: results.length, unique};
}

category('PolyTool: ChemEnum: dedup reduces count', () => {
  let rdkit: RDModule;

  before(async () => { rdkit = await getRdKitModule(); });

  test('duplicate cores: 6 enumerated collapse to 3 unique', async () => {
    // Two byte-identical cores mean every product is generated twice over.
    const cores = [makeCore('C[*:1]', 'dup-a'), makeCore('C[*:1]', 'dup-b')];
    const r1 = [makeRGroup('O[*:1]', 1, 'a'), makeRGroup('S[*:1]', 1, 'b'), makeRGroup('N[*:1]', 1, 'c')];
    const {raw, unique} = enumeratedVsUnique(cores, new Map([[1, r1]]), rdkit);
    expect(raw, 6); // raw enumerated rows
    expect(unique, 3); // after "Remove duplicates"
    expect(unique < raw, true); // the count actually drops
  });

  test('structural dedup: equivalent R-group SMILES collapse to one molecule', async () => {
    // Phenyl written aromatic vs kekulized is the SAME group — dedup is by canonical SMILES, not string.
    const cores = [makeCore('C[*:1]', 'c')];
    const r1 = [
      makeRGroup('c1ccccc1[*:1]', 1, 'aromatic'),
      makeRGroup('C1=CC=CC=C1[*:1]', 1, 'kekulized'),
      makeRGroup('Cl[*:1]', 1, 'chloro'),
    ];
    const {raw, unique} = enumeratedVsUnique(cores, new Map([[1, r1]]), rdkit);
    expect(raw, 3); // three R-groups → three enumerated rows
    expect(unique, 2); // the two phenyls are one molecule
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

// ─── Copy R-group list to another slot ──────────────────────────────────────

category('PolyTool: ChemEnum: copy R-group list', () => {
  let rdkit: RDModule;
  before(async () => { rdkit = await getRdKitModule(); });

  const rg = (smi: string, n: number) => makeRGroup(smi, n, '', rdkit);

  test('copyRGroupList re-labels to the target R#, appends/replaces, and no-ops on self-copy', async () => {
    const m = new Map([[1, [rg('O[*:1]', 1), rg('N', 1)]]]); // a labeled group + a single atom
    // Append into a new slot R2: the label is remapped [*:1]→[*:2]; the single atom keeps its
    // token but is re-targeted; the source is left untouched; the new slot is created.
    expect(copyRGroupList(m, 1, 2, 'append', rdkit), 2);
    expect(m.get(2)![0].smiles, 'O[*:2]');
    expect(m.get(2)![1].isSingleAtom, true);
    expect(m.get(2)![1].smiles, 'N');
    expect(m.get(2)![1].rNumber, 2);
    expect(m.get(1)!.length, 2, 'source must be unchanged');
    // Appending the same substituents again de-dupes (nothing new to add); replace overwrites.
    expect(copyRGroupList(m, 1, 2, 'append', rdkit), 0);
    expect(m.get(2)!.length, 2);
    expect(copyRGroupList(m, 1, 2, 'replace', rdkit), 2);
    expect(m.get(2)!.length, 2);
    // Copying a slot onto itself does nothing.
    expect(copyRGroupList(m, 1, 1, 'append', rdkit), 0);
    expect(m.get(1)!.length, 2);
  });

  test('shipped r-group-templates.json parses, every SMILES is valid, and inserts re-labeled', async () => {
    // The catalogue ships as a data file; parse it and check every SMILES is a valid R-group.
    const text = await _package.files.readAsText('enumeration/r-group-templates/r-group-templates.json');
    const templates = parseRGroupTemplates(text);
    expect(templates.length > 0, true, 'no templates parsed from r-group-templates.json');
    const bad = invalidTemplateSmiles(templates, rdkit);
    expect(bad.length, 0, `invalid template SMILES: ${bad.join(', ')}`);
    // The in-code fallback must be valid too (it's used when the file is unreadable).
    expect(invalidTemplateSmiles(BUILTIN_R_GROUP_TEMPLATES, rdkit).length, 0);
    // Insert the first template (alkyl C1–C8) into R2 — each lands re-labeled to [*:2].
    const smiles = templates[0].items.map((it) => it.smiles);
    const m = new Map<number, ReturnType<typeof makeRGroup>[]>();
    const n = addRGroupsFromSmiles(m, smiles, 2, 'append', rdkit);
    expect(n, smiles.length);
    expect(m.get(2)!.length, smiles.length);
    expect(m.get(2)!.every((g) => g.rNumber === 2 && g.error == null), true);
    expect(m.get(2)![0].smiles, 'C[*:2]'); // methyl, re-labeled
  });

  test('addRGroupsFromSmiles: replace with an empty list clears the slot; append no-ops', async () => {
    const m = new Map([[2, [rg('O[*:2]', 2)]]]);
    expect(addRGroupsFromSmiles(m, [], 2, 'replace', rdkit), 0);
    expect(m.has(2), false, 'replacing with nothing must delete the slot');
    expect(addRGroupsFromSmiles(m, [], 5, 'append', rdkit), 0);
    expect(m.has(5), false, 'appending nothing must not create an empty slot');
  });

  test('pickDefaultTargetR: free core slot, else lowest other, else exclude+1', async () => {
    expect(pickDefaultTargetR(new Set([1]), new Set([1, 2]), 1), 2); // free core R# wins
    expect(pickDefaultTargetR(new Set([1]), new Set([1, 3]), 1), 3); // non-contiguous core R#s
    expect(pickDefaultTargetR(new Set([1, 2]), new Set([1, 2]), 1), 2); // no free slot → lowest other
    expect(pickDefaultTargetR(new Set([1]), new Set(), 1), 2); // no cores → exclude + 1
    expect(pickDefaultTargetR(new Set([1]), new Set([1, 2])), 2); // template form (no exclude)
    expect(pickDefaultTargetR(new Set(), new Set()), 1); // empty → R1
  });

  test('rGroupTargetWarnings: flags invalid entries and an orphan target R#', async () => {
    expectArray(rGroupTargetWarnings(2, 0, new Set([1, 2])), []); // clean copy into a used slot
    expect(rGroupTargetWarnings(2, 3, new Set([1, 2])).length, 1); // invalid entries only
    expect(rGroupTargetWarnings(5, 0, new Set([1, 2])).length, 1); // orphan target only
    expect(rGroupTargetWarnings(5, 1, new Set([1, 2])).length, 2); // both
    expectArray(rGroupTargetWarnings(5, 0, new Set()), []); // no cores → suppress orphan warning
  });

  test('buildExportColumns: Core + R# columns, errored entries filtered, padded equally', async () => {
    const cores = [makeCore('C[*:1]N[*:2]', 'c', rdkit), makeCore('bad', 'e', rdkit)]; // 2nd errors
    expect(cores[1].error != null, true);
    const m = new Map([
      [1, [rg('O[*:1]', 1), rg('C1CC[*:1]', 1)]], // 2nd is unparseable → filtered out
      [2, [rg('N[*:2]', 2)]],
    ]);
    const cols = buildExportColumns(cores, m);
    expectArray(cols.map((c) => c.name), ['Core', 'R1', 'R2']);
    const valid = (name: string) => cols.find((c) => c.name === name)!.values.filter((v) => v !== '');
    expect(valid('Core').length, 1, 'errored core dropped'); // only C[*:1]N[*:2]
    expect(valid('R1').length, 1, 'invalid R1 entry dropped');
    expect(valid('R1')[0], 'O[*:1]');
    const len = cols[0].values.length;
    expect(cols.every((c) => c.values.length === len), true, 'all columns padded to equal length');
  });

  test('end-to-end: copying R1→R2 matches defining R2 natively', async () => {
    const core = makeCore('C[*:1]N[*:2]', 'c', rdkit);
    const r1 = [rg('O[*:1]', 1), rg('S[*:1]', 1)];
    // Copy path: populate R2 from R1 via the production helper.
    const byNum = new Map([[1, r1]]);
    copyRGroupList(byNum, 1, 2, 'append', rdkit);
    const copied = enumerate({cores: [core], rGroups: byNum, mode: ChemEnumModes.Cartesian}, rdkit)!;
    // Native path: R2 defined directly with [*:2] labels — must produce the same products.
    const r2Native = [rg('O[*:2]', 2), rg('S[*:2]', 2)];
    const native = enumerate(
      {cores: [core], rGroups: new Map([[1, r1], [2, r2Native]]), mode: ChemEnumModes.Cartesian}, rdkit)!;

    expect(copied.length, 4); // 2 × 2
    const canonSet = (rs: typeof copied) => new Set(rs.map((r) => canon(r.smiles, rdkit)));
    const copiedSet = canonSet(copied);
    expect(copiedSet.size, 4, 'copied products should be 4 distinct structures');
    for (const s of canonSet(native))
      expect(copiedSet.has(s), true, `copied products missing native product ${s}`);
    for (const r of copied)
      expect(r.smiles.includes('[*:'), false, `residual R-label in result: ${r.smiles}`);
    // The R2 result column shows each group re-labeled to its slot ([*:2]), not its source [*:1].
    for (const r of copied) {
      const r2 = r.rGroupSmilesByNum.get(2)!;
      expect(r2.includes('[*:2]') && !r2.includes('[*:1]'), true, `R2 column not re-labeled: ${r2}`);
    }
  });
});
