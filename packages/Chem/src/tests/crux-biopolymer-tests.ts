/**
 * Bio's macromolecule assembler on Crux (engine/biopolymer.ts, `Chem:getMacromoleculeAssembler`; Bio-on-Crux PLAN
 * §5.3): the worker jobs (crux/jobs/biopolymer.ts), CruxService.assembleMacromolecules (each worker's library, sent
 * once and again after a restart; the chunks; a long row alone on an idle worker), the main thread's one-row calls.
 * Crux in both engines, as molzip: there is no RDKit path. The batch helpers (crux/biopolymer-batch.ts) need no Crux.
 * The assembly's chemistry is crux-core's (its biopolymer crate's oracles); here a molecule is checked by its canonical
 * SMILES on the main thread's Crux module against a peptide written by hand.
 */
/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {AsmBatch, AsmConnection, AsmLayout, AsmRole, AsmSequence, AssembleResult, packBatch, unpackRow}
  from '@datagrok-libraries/chem-meta/src/macromolecule-assembler';
import type {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {PackageFunctions} from '../package';
import {CruxMacromoleculeAssembler} from '../engine/biopolymer';
import {loadCrux} from '../engine/crux-module';
import {getCruxService} from '../crux/crux-searches';
import {BIOPOLYMER_LONG_ROW} from '../crux/crux-service';
import {pickRows, sliceBatch} from '../crux/biopolymer-batch';
import {inBothEngines} from './engine-utils';

/** HELMCore's A, G and C as Bio's files/monomer-libraries/HELMCoreLibrary.json has them (its other fields left out). */
const TINY_MONOMERS = [
  {symbol: 'A', name: 'Alanine', molfile: '\n     RDKit          2D\n\n  7  6  0  0  0  0  0  0  0  0999 V2000\n    1.6702    1.3929    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.3712    0.6429    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.9279    1.3929    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.2269    0.6429    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n    0.3712   -0.8571    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.9279   -1.6071    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.6702   -1.6071    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n  2  1  1  1\n  2  3  1  0\n  3  4  1  0\n  2  5  1  0\n  5  6  2  0\n  5  7  1  0\nM  RGP  2   4   1   7   2\nM  END\n', smiles: 'C[C@H](N[*:1])C(=O)[*:2]', polymerType: 'PEPTIDE', monomerType: 'Backbone', naturalAnalog: 'A', rgroups: [{alternateId: 'R1-H', capGroupName: 'H', capGroupSMILES: '[*:1][H]', label: 'R1'}, {alternateId: 'R2-OH', capGroupName: 'OH', capGroupSMILES: 'O[*:2]', label: 'R2'}]},
  {symbol: 'G', name: 'Glycine', molfile: '\n     RDKit          2D\n\n  6  5  0  0  0  0  0  0  0  0999 V2000\n    1.0825    1.6250    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0825    0.1250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.2165   -0.6250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.5155    0.1250    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.8146   -0.6250    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n    2.3816   -0.6250    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  2  0\n  2  3  1  0\n  3  4  1  0\n  4  5  1  0\n  2  6  1  0\nM  RGP  2   5   1   6   2\nM  END\n', smiles: 'O=C(CN[*:1])[*:2]', polymerType: 'PEPTIDE', monomerType: 'Backbone', naturalAnalog: 'G', rgroups: [{alternateId: 'R1-H', capGroupName: 'H', capGroupSMILES: '[*:1][H]', label: 'R1'}, {alternateId: 'R2-OH', capGroupName: 'OH', capGroupSMILES: 'O[*:2]', label: 'R2'}]},
  {symbol: 'C', name: 'Cysteine', molfile: '\n     RDKit          2D\n\n  9  8  0  0  0  0  0  0  0  0999 V2000\n    1.4434   -2.1667    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.4434   -0.6667    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.1443    0.0833    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.1547   -0.6667    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.4537    0.0833    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.7528   -0.6667    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n    0.1443    1.5833    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    1.4434    2.3333    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n    2.7424    0.0833    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  2  0\n  2  3  1  0\n  3  4  1  6\n  4  5  1  0\n  5  6  1  0\n  3  7  1  0\n  7  8  1  0\n  2  9  1  0\nM  RGP  3   6   3   8   1   9   2\nM  END\n', smiles: 'O=C([C@H](CS[*:3])N[*:1])[*:2]', polymerType: 'PEPTIDE', monomerType: 'Backbone', naturalAnalog: 'C', rgroups: [{alternateId: 'R1-H', capGroupName: 'H', capGroupSMILES: '[*:1][H]', label: 'R1'}, {alternateId: 'R2-OH', capGroupName: 'OH', capGroupSMILES: 'O[*:2]', label: 'R2'}, {alternateId: 'R3-H', capGroupName: 'H', capGroupSMILES: '[*:3][H]', label: 'R3'}]},
];
const TINY = JSON.stringify(TINY_MONOMERS);
/** The same three with A named Ala: a second library. */
const RENAMED = JSON.stringify(TINY_MONOMERS.map((m) => m.symbol === 'A' ? {...m, symbol: 'Ala'} : m));
const HELM_CORE = 'System:AppData/Bio/monomer-libraries/HELMCoreLibrary.json';

/** The peptides of the tests, by hand (L-amino acids). */
const SMILES = {
  'A-G': 'C[C@H](N)C(=O)NCC(=O)O',
  'G-C': 'NCC(=O)N[C@@H](CS)C(=O)O',
  'A-G-C': 'C[C@H](N)C(=O)NCC(=O)N[C@@H](CS)C(=O)O',
  'G': 'NCC(=O)O',
};

/** One peptide chain of `symbols` ('' a gap), and its `connections`. */
function peptide(symbols: string[], connections: AsmConnection[] = []): AsmSequence {
  return {
    monomers: symbols.map((symbol, position) => symbol === '' ?
      {symbol: '', polymerType: 'PEPTIDE', role: AsmRole.AA, gap: true, position} :
      {symbol, polymerType: 'PEPTIDE', role: AsmRole.AA, position}),
    graphInfo: {connections, disjointSeqStarts: symbols.length > 0 ? [0] : []},
  };
}

let crux: RDModule;

/** The canonical SMILES of `molecule` on the main thread's Crux module ('none' for none). */
function canonical(molecule: string | null | undefined): string {
  if (!molecule)
    return 'none';
  const mol = crux.get_mol(molecule);
  if (!mol)
    return 'unread';
  try {
    return mol.get_smiles();
  } finally {
    mol.delete();
  }
}

/** A MOL block of `smiles` with 2D coordinates (V3000 with `v3`). */
function molblockOf(smiles: string, v3 = false): string {
  const mol = crux.get_mol(smiles)!;
  try {
    mol.set_new_coords();
    return v3 ? mol.get_v3Kmolblock() : mol.get_molblock();
  } finally {
    mol.delete();
  }
}

/** Row `i` of `result` is `smiles` (by canonical SMILES), with no error. */
function expectRow(result: AssembleResult, i: number, smiles: string, what: string): void {
  expect(result.errors[i], null, `${what}: row ${i}'s error`);
  expect(canonical(result.molblocks[i]), canonical(smiles), `${what}: row ${i}`);
}

/** Whether every index of `array` holds a value: none of `toList()`'s holes. */
const dense = (array: unknown[]): boolean => Array.from({length: array.length}, (_, i) => i in array).every((x) => x);

category('crux: biopolymer', () => {
  before(async () => {
    crux = await loadCrux('tests') as unknown as RDModule;
  });

  test('sliceBatch and pickRows: each row reads as in the whole batch', async () => {
    const twoChains: AsmSequence = {
      monomers: [...peptide(['C', 'A']).monomers, ...peptide(['C', 'G']).monomers.map((m) => ({...m, position: m.position + 2}))],
      graphInfo: {disjointSeqStarts: [0, 2],
        connections: [{seqIndex1: 0, monomerIndex1: 0, rGroup1: 3, seqIndex2: 1, monomerIndex2: 0, rGroup2: 3}]},
    };
    const inline: AsmSequence = peptide(['A', 'X', 'G']);
    inline.monomers[1] = {symbol: '', inlineSmiles: '[*:1]NCC(=O)[*:2]', polymerType: 'PEPTIDE', role: AsmRole.AA, position: 1};
    const hbond = peptide(['A', 'G']);
    hbond.graphInfo.connections.push({seqIndex1: 0, monomerIndex1: 0, rGroup1: 3, seqIndex2: 0, monomerIndex2: 1, rGroup2: 3, kind: 'hbond'});
    const seqs = [peptide(['A', 'G', 'C']), null, peptide([]), twoChains, peptide(['G', '', 'C']), inline, hbond,
      peptide(['G'])];
    const batch = packBatch(seqs, seqs.map((_, i) => i === 1 ? 'ambiguity' : null));
    const rowText = (b: AsmBatch, i: number) => JSON.stringify([unpackRow(b, i), b.rowError[i]]);
    const n = seqs.length;
    for (const [start, end] of [[0, n], [0, 1], [1, 4], [3, 3], [5, n], [n - 1, n]]) {
      const part = sliceBatch(batch, start, end);
      expect(part.rowStart.length, end - start + 1, `rows ${start}-${end}`);
      for (let i = 0; i < end - start; i++)
        expect(rowText(part, i), rowText(batch, start + i), `rows ${start}-${end}: row ${i}`);
      const used = new Set(Array.from(part.monSymbol).filter((s) => s >= 0));
      expect(used.size, part.symbols.length, `rows ${start}-${end}: only the symbols they use`);
      expect(part.rowStart.buffer !== batch.rowStart.buffer && part.conn.buffer !== batch.conn.buffer, true,
        'buffers of its own');
    }
    const rows = [0, 2, 3, 5, 7];
    const picked = pickRows(batch, rows);
    rows.forEach((row, k) => expect(rowText(picked, k), rowText(batch, row), `picked row ${row}`));
    const sparse = {...batch, rowError: new Array<string | null>(n)}; // a hole as no error
    expect(sliceBatch(sparse, 0, n).rowError.every((e) => e === null) && dense(sliceBatch(sparse, 0, n).rowError), true,
      'holes in rowError');
    let refused = 0;
    for (const bad of [() => sliceBatch(batch, 2, 1), () => sliceBatch(batch, 0, n + 1), () => pickRows(batch, [1, 1]),
      () => pickRows(batch, [3, 2])]) {
      try {
        bad();
      } catch {
        refused++;
      }
    }
    expect(refused, 4, 'rows out of a batch or out of order');
  });

  test('the library is compiled once per worker', async () => {
    const assembler = new CruxMacromoleculeAssembler();
    await assembler.setLibrary('tiny', TINY);
    const key = assembler.workerKeyForTests('tiny')!;
    const service = await getCruxService();
    const start = await service.biopolymerLibrariesForTests();
    const batch = packBatch([peptide(['A', 'G']), peptide(['G', 'C'])]);
    for (let call = 0; call < 3; call++) {
      const result = await assembler.assemble(batch, 'tiny', {layout: 'linear'});
      expectRow(result, 0, SMILES['A-G'], `call ${call}`);
      expectRow(result, 1, SMILES['G-C'], `call ${call}`);
    }
    const after = await service.biopolymerLibrariesForTests();
    expect(after.length, start.length, 'the pool was restarted');
    after.forEach((worker, w) => {
      expect(worker.compiled, start[w].compiled + 1, `worker ${w}: libraries compiled`);
      expect(worker.keys.includes(key), true, `worker ${w} holds the library`);
    });
    await assembler.setLibrary('tiny', TINY); // the same JSON again: the same key, nothing compiled
    expect(assembler.workerKeyForTests('tiny'), key, 'the key of the same library');
    await assembler.assemble(batch, 'tiny', {layout: 'linear'});
    const again = await service.biopolymerLibrariesForTests();
    expect(again.map((w) => w.compiled).join(), after.map((w) => w.compiled).join(), 'compiled again');
  }, {timeout: 60000});

  test('two libraries, and a key set again with another', async () => {
    const assembler = new CruxMacromoleculeAssembler();
    await assembler.setLibrary('tiny', TINY);
    await assembler.setLibrary('renamed', RENAMED);
    const withA = packBatch([peptide(['A', 'G'])]);
    const withAla = packBatch([peptide(['Ala', 'G'])]);
    expectRow(await assembler.assemble(withA, 'tiny', {layout: 'linear'}), 0, SMILES['A-G'], 'A in tiny');
    expectRow(await assembler.assemble(withAla, 'renamed', {layout: 'linear'}), 0, SMILES['A-G'], 'Ala in renamed');
    const noA = await assembler.assemble(withA, 'renamed', {layout: 'linear'});
    expect(!!noA.errors[0] && noA.molblocks[0] === null, true, 'A is not in renamed');
    const noAla = await assembler.assemble(withAla, 'tiny', {layout: 'linear'});
    expect(!!noAla.errors[0] && noAla.molblocks[0] === null, true, 'Ala is not in tiny');
    const keys = (await (await getCruxService()).biopolymerLibrariesForTests())[0].keys;
    for (const name of ['tiny', 'renamed'])
      expect(keys.includes(assembler.workerKeyForTests(name)!), true, `the first worker holds ${name}`);
    // 'tiny' set to the other library: a new key, the other library's answers
    const oldKey = assembler.workerKeyForTests('tiny');
    await assembler.setLibrary('tiny', RENAMED);
    expect(assembler.workerKeyForTests('tiny') !== oldKey, true, 'a new key for another JSON');
    expectRow(await assembler.assemble(withAla, 'tiny', {layout: 'linear'}), 0, SMILES['A-G'], 'tiny set to renamed');
    await expectRejection(() => assembler.assemble(withA, 'nope', {layout: 'linear'}), 'setLibrary first', 'no library');
  }, {timeout: 60000});

  test('a pool restart, a replaced worker, a worker that lost it: the library is sent again', async () => {
    const assembler = new CruxMacromoleculeAssembler();
    await assembler.setLibrary('tiny', TINY);
    const key = assembler.workerKeyForTests('tiny')!;
    const service = await getCruxService();
    // 4 rows per worker: every worker gets one chunk (a worker that lost the library learns it at its next chunk only)
    const rows = 4 * service.workerCount;
    const batch = packBatch(Array.from({length: rows}, (_, i) => peptide(i % 2 ? ['A', 'G'] : ['G', 'C'])));
    const check = async (what: string) => {
      const result = await assembler.assemble(batch, 'tiny', {layout: 'linear'});
      expectRow(result, 0, SMILES['G-C'], what);
      expectRow(result, rows - 1, SMILES['A-G'], what);
    };
    await check('first');
    const generation = service.generation;
    await service.failJobForTests().catch(() => {}); // a trap: every worker replaced
    expect(service.generation, generation + 1, 'the pool restarted');
    const fresh = await service.biopolymerLibrariesForTests();
    expect(fresh.every((w) => w.compiled === 0 && w.keys.length === 0), true, 'the new workers hold no library');
    await check('after the restart');
    const sent = await service.biopolymerLibrariesForTests();
    expect(sent.every((w) => w.compiled === 1 && w.keys.includes(key)), true, 'each new worker compiled it once');
    service.restartWorkerForTests(0);
    await check('after a worker was replaced');
    const replaced = await service.biopolymerLibrariesForTests();
    expect(replaced[0].compiled, 1, 'the replacement compiled it');
    await service.forgetBiopolymerLibraryForTests(key); // the workers drop it, their clients count it as sent
    await check('after the workers lost it');
    const lost = await service.biopolymerLibrariesForTests();
    lost.forEach((w, i) => expect(w.compiled, replaced[i].compiled + 1, `worker ${i}: compiled again once`));
    expect(service.generation, generation + 1, 'restarted again');
  }, {timeout: 120000});

  test('per-row errors: an unknown monomer, Datagrok\'s own; a library Crux does not read', async () => {
    const assembler = new CruxMacromoleculeAssembler();
    await assembler.setLibrary('tiny', TINY);
    const batch = packBatch([peptide(['A', 'G']), peptide(['A', 'Xyz', 'G']), null, peptide(['G', 'C'])],
      [null, null, 'ambiguity: (A,G)', null]);
    for (const layout of ['linear', 'nonlinear'] as const) {
      const result = await assembler.assemble(batch, 'tiny', {layout});
      expectRow(result, 0, SMILES['A-G'], layout);
      expect(result.errors[1]?.includes('Xyz'), true, `${layout}: the unknown monomer's row: ${result.errors[1]}`);
      expect(result.molblocks[1], null, `${layout}: no MOL block for it`);
      expect(result.layouts[1], AsmLayout.ERROR, `${layout}: its layout`);
      expect(result.errors[2], 'ambiguity: (A,G)', `${layout}: Datagrok's error`);
      expect(result.molblocks[2] === null && result.layouts[2] === AsmLayout.ERROR, true, `${layout}: Datagrok's row`);
      expectRow(result, 3, SMILES['G-C'], layout);
    }
    const service = await getCruxService();
    const generation = service.generation;
    await assembler.setLibrary('bad', '{"not": "a list of monomers"}');
    await expectRejection(() => assembler.assemble(batch, 'bad', {layout: 'linear'}), 'does not read', 'a bad library');
    expect(service.generation, generation, 'a bad library restarted the pool');
  }, {timeout: 60000});

  test('empty cells: empty sequences, holes in the rows\' errors and in a toList() of MOL blocks', async () => {
    const assembler = new CruxMacromoleculeAssembler();
    await assembler.setLibrary('tiny', TINY);
    const batch = packBatch([peptide(['A', 'G']), peptide([]), peptide(['G', 'C']), peptide([])]);
    batch.rowError = new Array(4); // holes: no errors
    const result = await assembler.assemble(batch, 'tiny', {layout: 'linear'});
    expect(dense(result.molblocks) && dense(result.errors) && result.layouts.length === 4, true, 'a value per row');
    expectRow(result, 0, SMILES['A-G'], 'holes');
    expectRow(result, 2, SMILES['G-C'], 'holes');
    for (const i of [1, 3])
      expect(result.errors[i] === null && !result.molblocks[i], true, `the empty sequence ${i}: ${result.errors[i]}`);
    const molblocks = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'molblocks',
      [molblockOf('C[C@H](N)C(=O)O'), null, '', molblockOf('CCO')]).toList() as (string | null)[];
    const marked = await assembler.markStereoAbsolute(molblocks);
    expect(marked.length === 4 && dense(marked), true, 'a value per cell');
    expect(canonical(marked[0]), canonical('C[C@H](N)C(=O)O'), 'the first cell');
    expect(!marked[1] && !marked[2], true, 'the empty cells stay empty');
    expect(canonical(marked[3]), 'CCO', 'the last cell');
  }, {timeout: 60000});

  test('2,000 rows in chunks across the workers, the same in both engines', async () => {
    const assembler = new CruxMacromoleculeAssembler();
    await assembler.setLibrary('tiny', TINY);
    const patterns: (keyof typeof SMILES)[] = ['A-G', 'G-C', 'A-G-C', 'G'];
    const batch = packBatch(Array.from({length: 2000}, (_, i) => peptide(patterns[i % 4].split('-'))));
    const {rdkit, crux: cruxMode} = await inBothEngines(() => assembler.assemble(batch, 'tiny', {layout: 'linear'}));
    for (const result of [rdkit, cruxMode]) {
      expect(result.molblocks.length === 2000 && result.errors.length === 2000 && result.layouts.length === 2000, true,
        'a value per row');
      expect(result.errors.filter((e) => e !== null).length, 0, `errors: ${result.errors.find((e) => e)}`);
      for (const i of [0, 1, 2, 3, 998, 999, 1000, 1001, 1996, 1997, 1998, 1999])
        expectRow(result, i, SMILES[patterns[i % 4]], 'row order');
    }
    expect(JSON.stringify(cruxMode.molblocks) === JSON.stringify(rdkit.molblocks), true, 'the same MOL blocks in both');
    expect(Array.from(cruxMode.layouts).join(), Array.from(rdkit.layouts).join(), 'the same layouts in both');
  }, {timeout: 120000});

  test('a long row goes alone to an idle worker', async () => {
    const assembler = new CruxMacromoleculeAssembler();
    await assembler.setLibrary('tiny', TINY);
    const length = BIOPOLYMER_LONG_ROW + 50;
    const batch = packBatch([peptide(['A', 'G']), peptide(new Array(length).fill('G')), peptide(['G', 'C'])]);
    const service = await getCruxService();
    const generation = service.generation;
    const sent = service.biopolymerLongRowsForTests;
    const result = await assembler.assemble(batch, 'tiny', {layout: 'linear'});
    expect(service.biopolymerLongRowsForTests, sent + 1, 'long rows sent alone');
    expectRow(result, 0, SMILES['A-G'], 'before the long row');
    expectRow(result, 1, 'NCC(=O)' + 'NCC(=O)'.repeat(length - 1) + 'O', 'the long row');
    expectRow(result, 2, SMILES['G-C'], 'after the long row');
    expect(service.generation, generation, 'the pool restarted');
  }, {timeout: 120000});

  test('assembleOne: the molecule and its monomer map, in both engines', async () => {
    const assembler = new CruxMacromoleculeAssembler();
    await assembler.setLibrary('tiny', TINY);
    const seq = peptide(['A', 'G', 'C']);
    const {rdkit, crux: one} = await inBothEngines(() => assembler.assembleOne(seq, 'tiny', {layout: 'linear', map: true}));
    expect(one.error, null, 'error');
    expect(canonical(one.molblock), canonical(SMILES['A-G-C']), 'the molecule');
    expect(one.molblock, rdkit.molblock, 'the same MOL block in both engines');
    expect(one.layout, AsmLayout.LINEAR, 'layout');
    const counts = /M {2}V30 COUNTS (\d+) (\d+)/.exec(one.molblock ?? '');
    expect(!!counts, true, 'a V3000 block');
    const map = one.map!;
    expect(map.atomPosition instanceof Int32Array && map.bondPosition instanceof Int32Array, true, 'Int32Arrays');
    expect(map.atomPosition.length, Number(counts![1]), 'a position per atom');
    expect(map.bondPosition.length, Number(counts![2]), 'a position per bond');
    const positions = new Set([...map.atomPosition, ...map.bondPosition]);
    expect([...positions].every((p) => p >= -1 && p <= 2), true, `positions: ${[...positions]}`);
    expect([0, 1, 2].every((p) => map.atomPosition.includes(p)), true, 'every monomer has atoms');
    expect(JSON.stringify(Array.from(map.atomPosition)), JSON.stringify(Array.from(rdkit.map!.atomPosition)), 'the same map');
    const noMap = await assembler.assembleOne(seq, 'tiny', {layout: 'linear'});
    expect(noMap.map === undefined, true, 'no map unless asked'); // (expect's `expected` defaults to true)
    const inBatch = await assembler.assemble(packBatch([seq]), 'tiny', {layout: 'linear'});
    expect(canonical(inBatch.molblocks[0]), canonical(one.molblock), 'the batch\'s molecule');
    // the workers read the options as the main thread does: 'auto' (CoordGen past a size) is no refusal there, and
    // a row's block is the batch's (the batch once sent `nonlinearEngine` as it came, and Crux's default switch)
    for (const nonlinearEngine of [undefined, 'depictor', 'auto'] as const) {
      const options = {layout: 'nonlinear' as const, nonlinearEngine};
      const batchRow = await assembler.assemble(packBatch([seq]), 'tiny', options);
      const row = await assembler.assembleOne(seq, 'tiny', options);
      expect(batchRow.errors[0], null, `${nonlinearEngine}: the batch's error`);
      expect(batchRow.molblocks[0], row.molblock, `${nonlinearEngine}: the batch's block is the row's`);
      expect(batchRow.layouts[0], row.layout, `${nonlinearEngine}: the same layout`);
    }
    const unknown = await assembler.assembleOne(peptide(['A', 'Xyz']), 'tiny', {layout: 'linear', map: true});
    expect(!!unknown.error && unknown.molblock === null && unknown.map === undefined, true, 'an unknown monomer');
    await expectRejection(() => assembler.assembleOne(seq, 'nope', {layout: 'linear'}), 'setLibrary first', 'no library');
  }, {timeout: 60000});

  test('markStereoAbsolute: one ABS group, V3000; toV2000', async () => {
    const assembler = new CruxMacromoleculeAssembler();
    const inputs = [molblockOf('C[C@H](N)C(=O)O'), molblockOf('N[C@@H](CS)C(=O)O', true), molblockOf('CCO'), null, ''];
    const {rdkit, crux: marked} = await inBothEngines(() => assembler.markStereoAbsolute(inputs));
    expect(JSON.stringify(marked), JSON.stringify(rdkit), 'the same in both engines');
    expect(marked.length, inputs.length, 'a value per input');
    for (const i of [0, 1, 2]) {
      expect(marked[i]!.includes('V3000'), true, `${i}: V3000`);
      expect(canonical(marked[i]), canonical(inputs[i]), `${i}: the molecule`);
      expect(marked[i]!.includes('MASS=1'), false, `${i}: no MASS=1`);
    }
    for (const i of [0, 1])
      expect(marked[i]!.includes('MDLV30/STEABS'), true, `${i}: the ABS group`);
    expect(marked[3], null, 'null');
    expect(marked[4], '', 'an empty MOL block');
    // a column's worth goes to the workers: the same blocks
    const column = await assembler.markStereoAbsolute(new Array(150).fill(inputs[0]));
    expect(column.every((m) => m === marked[0]), true, 'the workers\' blocks');
    const v2000 = await assembler.toV2000(marked[0]!);
    expect(v2000.includes('V2000') && !v2000.includes('V3000'), true, 'V2000');
    expect(canonical(v2000), canonical(inputs[0]), 'the V2000 block\'s molecule');
  }, {timeout: 60000});

  test('HELMCore (Bio\'s library): the same peptides as the inline library', async () => {
    const assembler = new CruxMacromoleculeAssembler();
    await assembler.setLibrary('HELMCore', await grok.dapi.files.readAsText(HELM_CORE));
    await assembler.setLibrary('tiny', TINY);
    const batch = packBatch([peptide(['A', 'G']), peptide(['G', 'C']), peptide(['A', 'G', 'C'])]);
    for (const layout of ['linear', 'nonlinear'] as const) {
      const core = await assembler.assemble(batch, 'HELMCore', {layout});
      const tiny = await assembler.assemble(batch, 'tiny', {layout});
      (['A-G', 'G-C', 'A-G-C'] as const).forEach((name, i) => {
        expectRow(core, i, SMILES[name], `${layout}, HELMCore`);
        expect(core.molblocks[i], tiny.molblocks[i], `${layout}: row ${i} as with the inline library`);
      });
    }
  }, {timeout: 60000});

  test('the synchronous twins: after ready(), as the async calls', async () => {
    const assembler = new CruxMacromoleculeAssembler();
    const seq = peptide(['A', 'G', 'C']);
    let refused = '';
    try {
      assembler.assembleOneSync(seq, 'tiny', {layout: 'linear'});
    } catch (e: any) {
      refused = String(e?.message ?? e);
    }
    expect(refused.includes('setLibrary first') || refused.includes('ready('), true, `before ready: ${refused}`);
    await assembler.ready('tiny', TINY);
    const one = assembler.assembleOneSync(seq, 'tiny', {layout: 'linear', map: true});
    const asyncOne = await assembler.assembleOne(seq, 'tiny', {layout: 'linear', map: true});
    expect(one.molblock, asyncOne.molblock, 'assembleOneSync: the same MOL block');
    expect(JSON.stringify(Array.from(one.map!.atomPosition)), JSON.stringify(Array.from(asyncOne.map!.atomPosition)),
      'assembleOneSync: the same map');
    const batch = packBatch([peptide(['A', 'G']), null, peptide(['A', 'Xyz']), peptide([])],
      [null, 'ambiguity', null, null]);
    const sync = assembler.assembleSync(batch, 'tiny', {layout: 'nonlinear'});
    const fromWorkers = await assembler.assemble(batch, 'tiny', {layout: 'nonlinear'});
    expectRow(sync, 0, SMILES['A-G'], 'assembleSync');
    expect(sync.molblocks[0], fromWorkers.molblocks[0], 'assembleSync: row 0 as the workers');
    expect(!!sync.errors[1] && !!sync.errors[2] && sync.molblocks[1] === null && sync.molblocks[2] === null, true,
      'assembleSync: row errors');
    expect(sync.molblocks[3], '', 'assembleSync: an empty row');
    expect(sync.layouts.length, 4, 'assembleSync: a layout per row');
    const inputs = [molblockOf('C[C@H](N)C(=O)O'), null, ''];
    expect(JSON.stringify(assembler.markStereoAbsoluteSync(inputs)), JSON.stringify(await assembler.markStereoAbsolute(inputs)),
      'markStereoAbsoluteSync');
    expect(assembler.toV2000Sync(inputs[0]!), await assembler.toV2000(inputs[0]!), 'toV2000Sync');
  }, {timeout: 60000});

  test('Chem:getMacromoleculeAssembler: one assembler for the session', async () => {
    const assembler = await PackageFunctions.getMacromoleculeAssembler();
    expect(assembler === await PackageFunctions.getMacromoleculeAssembler(), true, 'the same object');
    for (const method of ['setLibrary', 'assemble', 'assembleOne', 'markStereoAbsolute', 'toV2000', 'ready',
      'assembleOneSync', 'assembleSync', 'markStereoAbsoluteSync', 'toV2000Sync'] as const)
      expect(typeof assembler[method], 'function', method);
  });
});

/** `call()` rejects with a message holding `text`. */
async function expectRejection(call: () => Promise<unknown>, text: string, what: string): Promise<void> {
  const failure = await call().then(() => 'no failure', (e: any) => String(e?.message ?? e));
  expect(failure.includes(text), true, `${what}: ${failure}`);
}
