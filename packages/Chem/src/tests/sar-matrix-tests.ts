import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, test, expect, expectFloat, before, delay} from '@datagrok-libraries/test/src/test';

import {_package} from '../package-test';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {MmpFragments} from '../analysis/molecular-matched-pairs/mmp-analysis/mmpa-misc';
import {buildMatchedSeries, clusterRelatedCores} from '../analysis/sar-matrix/build/sar-matrix-clustering';
import {assembleSinglePositionMatrix, fitAdditiveModel,
  linkStaged} from '../analysis/sar-matrix/build/sar-matrix-assemble';
import {checkAgainstMolecules, cutWarning, decomposeByColumns, defaultAxis, holdsFragments, SarFragmentColumns,
  standardizeFragment}
  from '../analysis/sar-matrix/build/sar-matrix-columns';
import {cellPossible, LinkStages, planLink} from '../analysis/sar-matrix/build/sar-matrix-link';
import {computeMatrixConfidence} from '../analysis/sar-matrix/build/sar-matrix-confidence';
import {matrixCore} from '../analysis/sar-matrix/render/sar-matrix-depict';
import {SarRankScheme} from '../analysis/sar-matrix/build/sar-matrix-ranking';
import {fitRoleEffects, RoleDesign, RoleFit, RoleSummary, ROLE_FIT_MAX_SWEEPS}
  from '../analysis/sar-matrix/summary/sar-matrix-role-fit';
import {SummaryHost} from '../analysis/sar-matrix/summary/sar-matrix-summary-data';
import {SummaryPanel} from '../analysis/sar-matrix/summary/sar-matrix-summary-panel';
import {runSarMatrix, SarGrouping, SarMatrixParams} from '../analysis/sar-matrix/build/sar-matrix-run';
import {SCALING_METHODS} from '../analysis/molecular-matched-pairs/mmp-viewer/mmp-constants';
import {computeAllTransfers, spearman, transferStats} from '../analysis/sar-matrix/transfer/sar-matrix-transfer';
import {CoreCluster, MatchedSeries, SarMatrix, SarMatrixCell, SarMatrixColumn, SarMatrixRow}
  from '../analysis/sar-matrix/sar-matrix-types';

/** Minimal fake fragmentation: ids 1=CoreA, 2=CoreB, 3=Me, 4=Et (0 is the empty fragment). */
function fakeFrags(): MmpFragments {
  const idToName = ['', 'CoreA', 'CoreB', 'Me', 'Et'];
  const sizes = Uint32Array.from(idToName.map((n) => n.length));
  const fragCodes: [number, number][][] = [
    [[1, 3]], // mol 0: CoreA + Me
    [[1, 4]], // mol 1: CoreA + Et
    [[2, 3]], // mol 2: CoreB + Me
    [[2, 4]], // mol 3: CoreB + Et
  ];
  return {fragCodes, idToName, sizes};
}

function realCell(value: number, molIdx = 0, smiles: string | null = null): SarMatrixCell {
  return {kind: 'real', value, molIdx, smiles};
}

/**
 * Distinct probe compounds for the two sides of a transfer. These are fingerprinted, so they have to
 * parse; the side-'b' entries are close analogs of their side-'a' counterparts but never equal to them,
 * because one compound on both sides is skipped rather than matched.
 *
 * The R-groups come from `makeMatrix`, which labels column i 'Si' on every matrix — identical across
 * the two sides, which is what a pairing requires. Tests that are not about the compound-similarity
 * floor pass a threshold of 0, so their outcome cannot turn on a fingerprint value.
 */
const XFER_MOLS_A = ['Cc1ccccc1', 'CCCCCCCC', 'c1ccncc1', 'OC(=O)CCl', 'C1CCOC1', 'c1ccc2ccccc2c1',
  'c1ccc(F)cc1', 'C1CCCCC1'];
const XFER_MOLS_B = ['CCc1ccccc1', 'CCCCCCCCC', 'Cc1ccncc1', 'OC(=O)CCCl', 'C1CCOCC1', 'Cc1ccc2ccccc2c1',
  'c1ccc(Cl)cc1', 'C1CCCCCC1'];

/** A row of observed cells for one side of a transfer, column i carrying that side's probe i. */
function xferRow(values: number[], side: 'a' | 'b' = 'a'): SarMatrixCell[] {
  const mols = side === 'a' ? XFER_MOLS_A : XFER_MOLS_B;
  if (values.length > mols.length)
    throw new Error(`xferRow: ${values.length} columns asked for, only ${mols.length} probes defined`);
  return values.map((v, i) => realCell(v, (side === 'a' ? 0 : 100) + i, mols[i]));
}

function virtualCell(value: number, support = 1): SarMatrixCell {
  return {kind: 'virtual', value, molIdx: null, smiles: null, support};
}

function emptyCell(): SarMatrixCell {
  return {kind: 'empty', value: null, molIdx: null, smiles: null};
}

/** A bare-bones matrix wrapper around a hand-built cells grid, for testing the pure
 *  confidence/ranking functions without going through assembly. */
function makeMatrix(cells: SarMatrixCell[][], positions: string[] = ['R1']): SarMatrix {
  const rows: SarMatrixRow[] = cells.map((_, i) =>
    ({coreSmiles: `Core${i}`, keySmiles: `Core${i}`, label: `Core ${i}`, foldedValues: {}}));
  const columns: SarMatrixColumn[] = cells[0].map((_, i) => ({position: positions[0], substSmiles: `S${i}`}));
  let realCount = 0;
  for (const row of cells) {
    for (const cell of row) {
      if (cell.kind === 'real')
        realCount++;
    }
  }
  return {
    id: 'm', label: '', rows, columns, cells, minActivity: 0, maxActivity: 0, realCount, virtualCount: 0,
    scores: {}, positions, refValues: {},
    siteKey: '',
    level: 2,
  };
}

/** A fully-populated additive matrix: value[r][c] = rowEffect[r] + colEffect[c]. */
function additiveMatrix(rowEffect: number[], colEffect: number[]): SarMatrix {
  const cells = rowEffect.map((re) => colEffect.map((ce) => realCell(re + ce)));
  return makeMatrix(cells);
}

/** A 3x2 grid: three substituents at one ring position crossed with two at another, so the set has a
 *  shared scaffold and genuinely varies at two sites. */
function twoSiteGrid(): {molecules: DG.Column, activity: DG.Column<number>} {
  const smiles = ['Cc1ccc(F)cc1', 'CCc1ccc(F)cc1', 'Clc1ccc(F)cc1',
    'Cc1ccc(OC)cc1', 'CCc1ccc(OC)cc1', 'Clc1ccc(OC)cc1'];
  const molecules = DG.Column.fromStrings('smiles', smiles);
  molecules.semType = DG.SEMTYPE.MOLECULE;
  return {molecules, activity: DG.Column.fromList('double', 'activity', [5.1, 5.8, 6.4, 6.0, 6.7, 7.2])};
}

/** Four analogs of one para-disubstituted core, split the way an R-group decomposition leaves them:
 *  a core carrying every attachment point, and one column per R position. */
function fragmentTable(): {molecules: DG.Column, activity: DG.Column<number>, core: DG.Column,
  r1: DG.Column, r2: DG.Column} {
  const molecules = DG.Column.fromStrings('smiles',
    ['Cc1ccc(F)cc1', 'CCc1ccc(F)cc1', 'Cc1ccc(Cl)cc1', 'CCc1ccc(Cl)cc1']);
  molecules.semType = DG.SEMTYPE.MOLECULE;
  return {
    molecules,
    activity: DG.Column.fromList('double', 'activity', [5.1, 5.8, 6.4, 7.0]) as DG.Column<number>,
    core: DG.Column.fromStrings('Core', new Array(4).fill('[*:1]c1ccc([*:2])cc1')),
    r1: DG.Column.fromStrings('R1', ['C[*:1]', 'CC[*:1]', 'C[*:1]', 'CC[*:1]']),
    r2: DG.Column.fromStrings('R2', ['F[*:2]', 'F[*:2]', 'Cl[*:2]', 'Cl[*:2]']),
  };
}

function e2eParams(useMcsAnchors: boolean): SarMatrixParams {
  return {
    scaling: SCALING_METHODS.NONE, fragmentCutoff: 1, predictVirtual: true, grouping: SarGrouping.Site,
    fragmentationLevels: 2, higherIsBetter: true, threshold: 0.4, useMcsAnchors,
    // Below the viewer's default, so the fixture's exact matrix set survives the floor.
    rankScheme: SarRankScheme.Potency, minCompounds: 1, predictUnmeasured: true,
  };
}

/** Everything a rerun must reproduce: which matrices exist, and the cores/columns each is built from. */
function matrixShape(matrices: SarMatrix[]): string {
  return matrices.map((m) => [m.label, m.level, m.siteKey, m.positions.join('+'),
    m.rows.map((r) => r.coreSmiles).join(','), m.columns.map((c) => c.substSmiles).join(',')].join('|')).join('\n');
}

/** Two grid chemotypes the core grouping handles, plus isolated ring series it cannot: only the
 *  second kind is left for the MCS to pool. */
function mixedCoverage(): {molecules: DG.Column, activity: DG.Column<number>} {
  const smiles = ['Cc1ccc(F)cc1', 'CCc1ccc(F)cc1', 'Clc1ccc(F)cc1',
    'Cc1ccc(OC)cc1', 'CCc1ccc(OC)cc1', 'Clc1ccc(OC)cc1',
    'Cc1ccc2ccccc2c1', 'CCc1ccc2ccccc2c1', 'Clc1ccc2ccccc2c1', 'Brc1ccc2ccccc2c1',
    'Cc1cccc2ccccc12', 'CCc1cccc2ccccc12', 'Clc1cccc2ccccc12',
    'FC1CCCCC1', 'ClC1CCCCC1', 'BrC1CCCCC1'];
  const molecules = DG.Column.fromStrings('smiles', smiles);
  molecules.semType = DG.SEMTYPE.MOLECULE;
  return {molecules, activity: DG.Column.fromList('double', 'activity',
    smiles.map((_s, i) => 5 + (i % 5) * 0.4))};
}

/** Three cores x two substituents with every slot occupied but CoreC/Et, so the additive model can
 *  predict both blanks. A-Et is in the set and never assayed; only CoreC/Et is genuinely unmade. */
function unmeasuredMatrix(predictUnmeasured = false):
  {matrix: SarMatrix, etIdx: number, row: (core: string) => number} {
  const cluster: CoreCluster = {
    id: 'c0', siteKey: '', level: 2,
    series: [
      {coreSmiles: 'CoreA', members: [{molIdx: 0, substSmiles: 'Me'}, {molIdx: 1, substSmiles: 'Et'}]},
      {coreSmiles: 'CoreB', members: [{molIdx: 2, substSmiles: 'Me'}, {molIdx: 3, substSmiles: 'Et'}]},
      {coreSmiles: 'CoreC', members: [{molIdx: 4, substSmiles: 'Me'}]},
    ],
  };
  const matrix = assembleSinglePositionMatrix(cluster, ['A-Me', 'A-Et', 'B-Me', 'B-Et', 'C-Me'],
    Float32Array.from([10, NaN, 8, 6, 4]), true, predictUnmeasured);
  return {matrix, etIdx: matrix.columns.findIndex((c) => c.substSmiles === 'Et'),
    row: (core) => matrix.rows.findIndex((r) => r.coreSmiles === core)};
}

category('SAR Matrix', () => {
  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  test('series by shared core', async () => {
    const series = buildMatchedSeries(fakeFrags(), 10);
    expect(series.length, 2);
    for (const s of series) {
      expect(s.members.length, 2);
      const subs = s.members.map((m) => m.substSmiles).sort();
      expect(subs[0], 'Et');
      expect(subs[1], 'Me');
    }
  });

  test('columns from substituents', async () => {
    const cluster: CoreCluster = {
      id: 'c0',
      siteKey: '',
      level: 2,
      series: [
        {coreSmiles: 'CoreA', members: [
          {molIdx: 0, substSmiles: 'Me'},
          {molIdx: 1, substSmiles: 'Et'},
          {molIdx: 2, substSmiles: 'Pr'},
        ]},
        {coreSmiles: 'CoreB', members: [
          {molIdx: 3, substSmiles: 'Me'},
        ]},
      ],
    };
    const molecules = ['A-Me', 'A-Et', 'A-Pr', 'B-Me'];
    const activities = Float32Array.from([1, 2, 3, 4]);

    const noPredict = assembleSinglePositionMatrix(cluster, molecules, activities, false);
    expect(noPredict.columns.map((c) => c.substSmiles).sort().join(','), 'Et,Me,Pr');
    expect(noPredict.realCount, 4);
    expect(noPredict.virtualCount, 0);
    const hasVirtual = noPredict.cells.some((row) => row.some((cell) => cell.kind === 'virtual'));
    expect(hasVirtual, false, 'predict:false must leave missing combinations empty, not virtual');

    const predicted = assembleSinglePositionMatrix(cluster, molecules, activities, true);
    const etIdx = predicted.columns.findIndex((c) => c.substSmiles === 'Et');
    const prIdx = predicted.columns.findIndex((c) => c.substSmiles === 'Pr');
    const filled = [predicted.cells[1][etIdx], predicted.cells[1][prIdx]];
    expect(filled.some((c) => c.kind === 'virtual' && c.value !== null), true,
      'predict:true fills at least one missing combination');
  });

  test('NaN activity is unmeasured', async () => {
    const cluster: CoreCluster = {
      id: 'c0',
      siteKey: '',
      level: 2,
      series: [
        {coreSmiles: 'CoreA', members: [
          {molIdx: 0, substSmiles: 'Me'},
          {molIdx: 1, substSmiles: 'Et'}, // activity missing (NaN)
        ]},
        // B-Et is what keeps the Et column observed, and so unpruned: a column whose only compound is
        // the unmeasured one is an axis the data never touched, and is dropped before this assertion.
        {coreSmiles: 'CoreB', members: [
          {molIdx: 2, substSmiles: 'Me'},
          {molIdx: 3, substSmiles: 'Et'},
        ]},
      ],
    };
    const molecules = ['A-Me', 'A-Et', 'B-Me', 'B-Et'];
    const activities = Float32Array.from([10, NaN, 5, 7]);
    const matrix = assembleSinglePositionMatrix(cluster, molecules, activities, false);

    expect(matrix.realCount, 3, 'the NaN-activity member must not count as a real observation');
    const etIdx = matrix.columns.findIndex((c) => c.substSmiles === 'Et');
    const aEt = matrix.cells[0][etIdx];
    expect(aEt.kind, 'unmeasured', 'the compound is in the set, only its activity is missing');
    expect(aEt.value, null);
    expect(aEt.molIdx, 1, 'an unmeasured cell keeps the compound it stands for');
    // A poisoned min (e.g. 0 from an unguarded NaN comparison) would break the potency color scale.
    expect(matrix.minActivity, 5, 'minActivity must be the min of the finite activities only');
    expect(matrix.maxActivity, 10, 'maxActivity must be the max of the finite activities only');
  });

  test('untested is not an analog', async () => {
    const {matrix, etIdx, row} = unmeasuredMatrix();
    expect(matrix.cells[row('CoreA')][etIdx].kind, 'unmeasured', 'a compound the set holds is not one to make');
    expect(matrix.cells[row('CoreA')][etIdx].value, null, 'with prediction off, no number may be reported');
    expect(matrix.cells[row('CoreC')][etIdx].kind, 'virtual', 'an unmade slot must still be predicted');
    expect(matrix.cells.flat().filter((c) => c.kind === 'virtual').length, 1, 'only the unmade slot is proposed');
  });

  test('prune unobserved core', async () => {
    // CoreC bears only a compound with no activity, so its row can never be predicted from anything.
    const cluster: CoreCluster = {
      id: 'c0', siteKey: '', level: 2,
      series: [
        {coreSmiles: 'CoreA', members: [{molIdx: 0, substSmiles: 'Me'}, {molIdx: 1, substSmiles: 'Et'}]},
        {coreSmiles: 'CoreB', members: [{molIdx: 2, substSmiles: 'Me'}, {molIdx: 3, substSmiles: 'Et'}]},
        {coreSmiles: 'CoreC', members: [{molIdx: 4, substSmiles: 'Me'}]},
      ],
    };
    const matrix = assembleSinglePositionMatrix(cluster, ['A-Me', 'A-Et', 'B-Me', 'B-Et', 'C-Me'],
      Float32Array.from([10, 9, 8, 6, NaN]), true);
    expect(matrix.rows.length, 2, 'the core with nothing measured on it is not an axis');
    expect(matrix.rows.some((r) => r.coreSmiles === 'CoreC'), false);
    expect(matrix.rows.map((r) => r.label).join(','), 'Core 1,Core 2', 'labels renumber over the gap');
    expect(matrix.cells.flat().filter((c) => c.kind === 'empty').length, 0,
      'every surviving cell is measured or predictable');
    expect(matrix.realCount, 4, 'dropping an unobserved row loses no measured compound');
  });

  test('predict untested compound', async () => {
    const {matrix, etIdx, row} = unmeasuredMatrix(true);
    const aEt = matrix.cells[row('CoreA')][etIdx];
    expect(aEt.kind, 'unmeasured', 'a predicted number must not promote the compound to an unmade analog');
    expect(aEt.molIdx, 1, 'the cell still stands for the compound the set holds');
    // Four measurements against four identifiable effects, so the fit reproduces every one of them
    // exactly. Et reads -2 against Me where both are measured (CoreB: 6 vs 8), and CoreA measures 10 at
    // Me, so A-Et is 8 — the only value with no residual anywhere.
    expect(aEt.value !== null && Math.abs(aEt.value - 8) < 1e-6, true,
      'CoreA at Me (10) shifted by the Et-vs-Me step (-2)');
    expect(matrix.cells.flat().filter((c) => c.kind === 'virtual').length, 1,
      'the count of analogs to synthesize is unchanged');
  });

  test('additive model fit', async () => {
    // 3x3 grid: row2 and col2 have no observation at all.
    const cells: SarMatrixCell[][] = [
      [realCell(1), realCell(3), emptyCell()],
      [realCell(5), emptyCell(), emptyCell()],
      [emptyCell(), emptyCell(), emptyCell()],
    ];
    const predict = fitAdditiveModel(cells, 3, 3);

    // Column 1 reads +2 against column 0 where both are measured (row 0: 3 vs 1), and row 1 measures 5
    // at column 0, so the additive answer at (1,1) is 7 — the one value that reproduces all three
    // measurements exactly. Averaging the margins in a single pass answers 5 instead, because row 0's
    // mean folds column 1's +2 into row 0's own effect.
    const p11 = predict(1, 1);
    expect(p11 !== null, true);
    expectFloat(p11!.value, 7);
    expect(p11!.support, 1, 'support = min(rowN=1, colN=1)');

    expect(predict(0, 2), null, 'column 2 has no observation — unpredictable');
    expect(predict(2, 0), null, 'row 2 has no observation — unpredictable');
  });

  test('confidence needs 4 cells', async () => {
    const cells: SarMatrixCell[][] = [
      [realCell(1), realCell(2), emptyCell()],
      [realCell(3), emptyCell(), emptyCell()],
    ];
    expect(computeMatrixConfidence(makeMatrix(cells)), null);
  });

  test('confidence on additive', async () => {
    const matrix = additiveMatrix([0, 1, 2, 3], [0, 10, 20, 30]);
    const conf = computeMatrixConfidence(matrix);
    expect(conf !== null, true);
    expect(conf!.n, 16);
    expect(conf!.r2 > 0.85, true, `expected r2 close to 1, got ${conf!.r2}`);
    expect(conf!.rmse < 5, true, `expected a small rmse, got ${conf!.rmse}`);
  });

  test('confidence on non-additive', async () => {
    // Diagonal spikes: no row/column effect explains this pattern.
    const cells: SarMatrixCell[][] = Array.from({length: 4}, (_, ri) =>
      Array.from({length: 4}, (_, ci) => realCell(ri === ci ? 100 : 0)));
    const conf = computeMatrixConfidence(makeMatrix(cells));
    expect(conf !== null, true);
    expect(conf!.r2 < 0, true, `expected a negative r2 for a non-additive matrix, got ${conf!.r2}`);
  });

  test('confidence per R-position', async () => {
    // Two position groups (R1 = cols 0-2, R2 = cols 3-5), each additive on its own.
    const cells: SarMatrixCell[][] = [
      [realCell(1), realCell(2), realCell(3), realCell(10), realCell(20), realCell(30)],
      [realCell(2), realCell(3), realCell(4), realCell(20), realCell(30), realCell(40)],
    ];
    const mat = makeMatrix(cells, ['R1', 'R2']);
    mat.columns.forEach((c, i) => c.position = i < 3 ? 'R1' : 'R2');
    const conf = computeMatrixConfidence(mat);
    expect(conf !== null, true);
    expect(conf!.total, 12, 'all 12 observed cells across both position slices are counted');
    expect(conf!.n, 12, 'every cell is cross-validatable within its own slice');
    expect(conf!.r2 > 0.85, true, `both slices are additive, so r2 should be high, got ${conf!.r2}`);
  });

  test('spearman correlation', async () => {
    // A monotone but non-linear relation is a perfect rank correlation even though Pearson is < 1.
    expectFloat(spearman([1, 2, 3, 4], [1, 4, 9, 16])!, 1, 0.001);
    expectFloat(spearman([1, 2, 3, 4], [4, 3, 2, 1])!, -1, 0.001);
    // One extreme value can dominate Pearson but not the ranks — the ordering is still 1:1.
    expectFloat(spearman([1, 2, 3, 100], [1, 2, 3, 4])!, 1, 0.001);
    // Guards: too few shared points and a constant side both give null.
    expect(spearman([1, 2], [1, 2]), null, 'fewer than the minimum shared points is null');
    expect(spearman([1, 1, 1, 1], [1, 2, 3, 4]), null, 'a constant side has no ordering to correlate');
  });

  test('transfer across matrices', async () => {
    // Two matrices whose columns carry the same R-groups (S0..S3); a row in each tracks the other.
    const matA = makeMatrix([
      xferRow([1, 2, 3, 4]),
      xferRow([5, 4, 3, 2]), // anti-correlated with matB's row, so it isn't the pick
    ]);
    const matB = makeMatrix([xferRow([2, 3, 4, 5], 'b')]);
    // Threshold 0 opens the compound-similarity gate, so the outcome rests on the R-groups alone.
    const transfers = await computeAllTransfers([matA, matB], 0);
    expect(transfers.length, 1, 'a transfer between the two matrices must be found');
    expect(transfers[0].a.matrixIndex !== transfers[0].b.matrixIndex, true,
      'the two cores come from different matrices');
    expect(transfers[0].substituents.length, 4, 'all four shared R-groups are used');
    expectFloat(transfers[0].correlation, 1, 0.01);
  });

  test('transfer correlation floor', async () => {
    const matA = makeMatrix([xferRow([1, 2, 3, 4])]);
    // Ranks 2,1,4,3 against 1,2,3,4 is a spearman of 0.6 — a real trend, under the floor for four points.
    const below = makeMatrix([xferRow([2, 1, 4, 3], 'b')]);
    expect((await computeAllTransfers([matA, below], 0)).length, 0, 'a correlation under the floor is no transfer');
    // Same shapes and the same open similarity gate, so a passing correlation must still come through:
    // without this the assertion above would hold even if the pairing had failed for some other reason.
    const above = makeMatrix([xferRow([2, 3, 4, 5], 'b')]);
    expect((await computeAllTransfers([matA, above], 0)).length, 1, 'the pair is otherwise eligible');
  });

  test('transfer dedupe by position', async () => {
    // The same two cores pair at both R1 (ρ=1) and R2 (ρ=0.8) — only R1 survives.
    const positions = ['R1', 'R2'];
    const twoPosition = (values: number[], side: 'a' | 'b'): SarMatrix => {
      const m = makeMatrix([xferRow(values, side)], positions);
      m.columns.forEach((c, i) => c.position = i < 4 ? 'R1' : 'R2');
      return m;
    };
    const matA = twoPosition([1, 2, 3, 4, 1, 2, 4, 3], 'a');
    const matB = twoPosition([2, 3, 4, 5, 1, 2, 3, 4], 'b');
    const transfers = await computeAllTransfers([matA, matB], 0);
    expect(transfers.length, 1, 'the two R-positions collapse to one entry for this core pair');
    expect(transfers[0].a.position, 'R1', 'R1 is the more-correlated position, so it is kept');
  });

  test('foldMatch on equal steps',
    async () => {
      // matB = matA − 9, so the trends are identical: correlation 1 and every step delta matches.
      const matA = makeMatrix([xferRow([10, 12, 14, 20])]); // deltas +2, +2, +6
      const matB = makeMatrix([xferRow([1, 3, 5, 11], 'b')]); // same deltas
      const transfers = await computeAllTransfers([matA, matB], 0);
      expect(transfers.length, 1, 'one transfer between the two single-row matrices');
      const stats = transferStats(transfers[0], true);
      expectFloat(stats.correlation, 1, 0.01);
      expect(stats.foldMatch !== null, true);
      expectFloat(stats.foldMatch!, 1, 0.01);
      expect(stats.benefiting, null, 'no virtual cells in either core — nothing to benefit from the transfer');
    });

  test('foldMatch on unequal steps', async () => {
    const matA = makeMatrix([xferRow([10, 12, 14, 20])]); // deltas +2, +2, +6
    const matB = makeMatrix([xferRow([1, 2, 3, 4], 'b')]); // same direction, +1 each
    const transfers = await computeAllTransfers([matA, matB], 0);
    expect(transfers.length, 1);
    const stats = transferStats(transfers[0], true);
    expect(stats.foldMatch !== null, true);
    // steps: min(2,1)/max=0.5, min(2,1)/max=0.5, min(6,1)/max=1/6  ->  mean ~0.389
    expectFloat(stats.foldMatch!, 0.389, 0.01);
    expect(stats.foldMatch! < 1, true);
  });

  test('transferStats both sides', async () => {
    // Each core is untested at the R-group the other one measured, so the pairing argues for an analog
    // on both sides. The stronger is the leader's, and the leader is scanned second, so stopping at the
    // first side with a candidate would report the weaker analog as the one to make.
    const matA = makeMatrix([[...xferRow([1, 2, 3, 4]), virtualCell(20),
      realCell(6, 5, XFER_MOLS_A[5])]]);
    const matB = makeMatrix([[...xferRow([2, 3, 4, 5], 'b'), realCell(7, 104, XFER_MOLS_B[4]),
      virtualCell(10)]]);
    const transfers = await computeAllTransfers([matA, matB], 0);
    expect(transfers.length, 1, 'the four measured pairs correlate, so the pair is a transfer');
    const stats = transferStats(transfers[0], true);
    expect(stats.benefiting !== null, true, 'the pairing argues for an analog on each side');
    expect(stats.benefiting!.side, 'a', 'the strongest analog wins wherever it sits');
    expectFloat(stats.benefiting!.value, 20, 0.01);
  });

  test('transfer carries analogs', async () => {
    // The leader measured all five R-groups; the follower has the fifth only as a prediction. That
    // column is what the transfer argues for, so it comes back alongside the four that are matched —
    // and separately from them, since it has no second observation to correlate.
    const matA = makeMatrix([xferRow([1, 2, 3, 4, 5])]);
    const matB = makeMatrix([[...xferRow([2, 3, 4, 5], 'b'), virtualCell(10)]]);
    const transfers = await computeAllTransfers([matA, matB], 0);
    expect(transfers.length, 1);
    expect(transfers[0].substituents.length, 4, 'four measured pairs carry the correlation');
    expect(transfers[0].predictedSubstituents.length, 1, 'the predicted analog is carried separately');
    expect(transfers[0].predictedSubstituents[0], 'S4', 'the R-group the follower has not made');
    expect(transfers[0].predictedACols[0], 4, 'the leader measured it at its fifth column');
    expect(transfers[0].predictedBCols[0], 4, 'the follower has it predicted at its fifth column');
  });

  test('transfer ignores shared cpd', async () => {
    // The overlapping cover puts one compound in several matrices at once. Matched against itself it
    // would put its own potency on both axes and report a perfect correlation that means nothing, so
    // these two series — same compounds, same potencies — must yield no transfer at all.
    const matA = makeMatrix([xferRow([1, 2, 3, 4])]);
    const matB = makeMatrix([xferRow([1, 2, 3, 4])]);
    expect((await computeAllTransfers([matA, matB], 0)).length, 0,
      'a series compared against its own compounds is not a transfer');
  });

  test('transfer similarity floor', async () => {
    // Identical R-groups and a perfect trend either way, so only the similarity gate decides. At 0
    // every pair passes; at 1 none can, since the two sides are different compounds by construction.
    const matA = makeMatrix([xferRow([1, 2, 3, 4])]);
    const matB = makeMatrix([xferRow([2, 3, 4, 5], 'b')]);
    expect((await computeAllTransfers([matA, matB], 0)).length, 1, 'an open gate admits the transfer');
    expect((await computeAllTransfers([matA, matB], 1)).length, 0,
      'demanding all but identical compounds admits none');
  });

  test('link core and fragment', async () => {
    const svc = await chemCommonRdKit.getRdKitService();
    const linked = await svc.linkRGroupFragments(['c1ccc([*:1])cc1'], [['C[*:1]']], [1]);
    expect(linked.length, 1);
    expect(linked[0].includes('*'), false, 'no stray attachment point should remain');
    const mol = chemCommonRdKit.checkMoleculeValid(linked[0]);
    expect(mol !== null, true, 'the assembled SMILES must parse as a valid molecule');
    expect(mol.get_num_atoms(), 7, 'toluene: 6 ring carbons + 1 methyl carbon');
    mol.delete();
  });

  test('link skips missing position', async () => {
    const svc = await chemCommonRdKit.getRdKitService();
    // Rows of one matrix can have different cores, so a core need not carry every decomposed position.
    const linked = await svc.linkRGroupFragments(
      ['c1ccc([*:1])cc1'], [['C[*:1]'], ['Cl[*:2]']], [1, 2]);
    expect(linked.length, 1);
    expect(linked[0] !== '', true, 'an absent attachment point must not discard the whole structure');
    const mol = chemCommonRdKit.checkMoleculeValid(linked[0]);
    expect(mol !== null, true, 'the assembled SMILES must parse as a valid molecule');
    expect(mol.get_num_atoms(), 7, 'toluene: R1 attached, R2 skipped because the core has no [*:2]');
    mol.delete();
  });

  // An R-group written with an isotope dummy comes back from RDKit with the label FIRST, so the
  // linker has to move it into a branch — and that move rewrites the neighbour order the chiral tag
  // was measured against. Getting it wrong draws the mirror image of the compound that was assayed.
  test('link preserves stereochemistry', async () => {
    const svc = await chemCommonRdKit.getRdKitService();
    const rdkit = chemCommonRdKit.getRdKitModule();
    // Substituting the label in place keeps the neighbour order the fragment wrote, so this is the
    // molecule the join must reproduce, whatever route it takes to build it.
    const truth = (fragment: string): string => {
      const mol = rdkit.get_mol(fragment.split('[*:1]').join('c%11ccccc%11'));
      const smiles = mol.get_smiles();
      mol.delete();
      return smiles;
    };
    const fragments = ['[*:1][C@@H](C)O', '[*:1][C@H](C)O', '[*:1]/C=C/C', '[*:1][C@](C)(F)Cl'];
    const linked = await svc.linkRGroupFragments(
      new Array(fragments.length).fill('c%11ccccc%11[*:1]'), [fragments], [1]);
    for (let i = 0; i < fragments.length; i++) {
      const mol = chemCommonRdKit.checkMoleculeValid(linked[i]);
      expect(mol !== null, true, `${fragments[i]} must assemble into a valid molecule`);
      expect(mol.get_smiles(), truth(fragments[i]), `${fragments[i]} must keep its configuration`);
      mol.delete();
    }
  });

  test('run is deterministic', async () => {
    const {molecules, activity} = twoSiteGrid();
    const first = await runSarMatrix(molecules, activity, e2eParams(false));
    const second = await runSarMatrix(molecules, activity, e2eParams(false));
    expect(first.length > 0, true, 'the fragmentation/clustering/assembly pipeline must yield matrices');
    expect(matrixShape(first), matrixShape(second), 'two runs over one input must give identical matrices');
  });

  // The option only ever adds: it pools the series no shared core could group and asks an MCS for
  // their common core. Whatever the core grouping found must survive it untouched.
  // Compared by cluster id, not by content: an anchor the MCS finds can reshape a matrix's rows, and
  // that is the option working. What must never happen is a cluster losing its matrix, or a compound
  // the core grouping placed ending up in none.
  test('MCS costs no matrix', async () => {
    const {molecules, activity} = mixedCoverage();
    const params = {...e2eParams(false), grouping: SarGrouping.Similarity};
    const off = await runSarMatrix(molecules, activity, params);
    const on = await runSarMatrix(molecules, activity, {...params, useMcsAnchors: true});
    const placed = (matrices: SarMatrix[]): Set<number> => {
      const mols = new Set<number>();
      for (const matrix of matrices) {
        for (const row of matrix.cells) {
          for (const cell of row) {
            if (cell.kind === 'real' && cell.molIdx !== null)
              mols.add(cell.molIdx);
          }
        }
      }
      return mols;
    };
    const onIds = new Set(on.map((matrix) => matrix.id));
    const vanished = off.filter((matrix) => !onIds.has(matrix.id));
    const onMols = placed(on);
    const droppedMols = [...placed(off)].filter((molIdx) => !onMols.has(molIdx));
    expect(off.length > 0, true, 'the core grouping must produce matrices on its own');
    expect(vanished.length, 0, 'every cluster with a matrix must still have one once the MCS is on');
    expect(droppedMols.length, 0, 'no compound the core grouping placed may end up in no matrix');
    expect(on.length >= off.length, true, 'pooling the leftovers cannot yield fewer matrices');
  });

  // A decomposed matrix still shows one column axis, so `positions` is length 1 either way; what
  // separates it from the single-position fallback is that the decomposition named every position it
  // found, which is what `refValues` carries.
  test('site key anchors matrix', async () => {
    const {molecules, activity} = twoSiteGrid();
    const matrices = await runSarMatrix(molecules, activity, e2eParams(false));
    const decomposed = matrices.filter((m) => Object.keys(m.refValues).length > 1);
    expect(decomposed.length > 0, true, 'a cluster grouped by site carries its shared scaffold, so it ' +
      'must decompose against it with no MCS rather than falling back to a single position');
  });

  // Fragment ids are minted as the workers discover fragments, so the same data reaches this stage
  // under different ids from one run to the next. Anything downstream resolving a tie by "whichever
  // came first" then yields a different set of matrices for identical input.
  test('series order by chemistry', async () => {
    // Same four molecules and the same two cores, with the core ids swapped and the rows reversed.
    const asDiscovered: MmpFragments = {
      idToName: ['', 'CoreA', 'CoreB', 'Me', 'Et'],
      sizes: Uint32Array.from(['', 'CoreA', 'CoreB', 'Me', 'Et'].map((n) => n.length)),
      fragCodes: [[[1, 3]], [[1, 4]], [[2, 3]], [[2, 4]]],
    };
    const asRediscovered: MmpFragments = {
      idToName: ['', 'CoreB', 'CoreA', 'Et', 'Me'],
      sizes: Uint32Array.from(['', 'CoreB', 'CoreA', 'Et', 'Me'].map((n) => n.length)),
      fragCodes: [[[2, 4]], [[2, 3]], [[1, 4]], [[1, 3]]],
    };
    const shape = (frags: MmpFragments): string => buildMatchedSeries(frags, 1)
      .map((s) => `${s.coreSmiles}:${s.members.map((m) => m.molIdx).join(',')}`).join(' | ');
    expect(shape(asDiscovered), shape(asRediscovered),
      'the same molecules must give the same series whichever ids the workers assigned');
  });

  // Clustering by core similarity alone will put cores whose substituents hang off different places
  // into one matrix. Their substituent vocabularies are disjoint, so the column axis pools both and a
  // column means one position on some rows and another on the rest.
  test('clusters share a site', async () => {
    // One triazole, attachments at three different ring positions — near-identical by fingerprint.
    const series: MatchedSeries[] = [
      {coreSmiles: 'CCc1c([*:1])nnn1-c1ccc(F)cc1', members: [{molIdx: 0, substSmiles: 'C[*:1]'}]},
      {coreSmiles: 'COc1c([*:1])nnn1-c1ccc(F)cc1', members: [{molIdx: 1, substSmiles: 'C[*:1]'}]},
      {coreSmiles: 'COC(=O)c1nnn(-c2ccc(F)cc2)c1[*:1]', members: [{molIdx: 2, substSmiles: 'C[*:1]'}]},
      {coreSmiles: 'N#Cc1nnn(-c2ccc(F)cc2)c1[*:1]', members: [{molIdx: 3, substSmiles: 'CC[*:1]'}]},
    ];
    const clusters = await clusterRelatedCores(series, 0.3);
    const grouped = clusters.filter((c) => c.series.length > 1);
    // Without this the assertion below is vacuous: all-singleton output would pass it having tested
    // nothing, and singletons are exactly what a too-tight threshold produces.
    expect(grouped.length > 0, true, 'these cores must cluster at all for the check to mean anything');
    for (const cluster of grouped) {
      expect(cluster.siteKey !== '', true,
        'a cluster holding several cores must name the site they share, or its columns mix positions');
      for (const matched of cluster.series) {
        expect(matched.coreSmiles.includes('[*:'), true,
          'every core in a clustered series must carry the attachment its substituents hang off');
      }
    }
  });
});

/** Building the matrices from columns that already hold a decomposition, instead of fragmenting the
 *  structures. Kept apart from the pipeline category: nothing here fragments anything. */
category('SAR Matrix: R-group columns', () => {
  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  const col = (name: string, values: string[]): DG.Column => DG.Column.fromStrings(name, values);
  const measured = (n: number): Float32Array => new Float32Array(n).fill(1);

  /** Decomposition columns for a table. A plain string is the same fragment in every row; the last
   *  group runs across the top unless `axis` names another. */
  const spec = (core: string, groups: {[name: string]: string[] | string},
    axis?: string): SarFragmentColumns => {
    const names = Object.keys(groups);
    const rows = Math.max(...names.map((k) => Array.isArray(groups[k]) ? groups[k].length : 0));
    const pick = (k: string): DG.Column =>
      col(k, Array.isArray(groups[k]) ? groups[k] : new Array(rows).fill(groups[k]));
    const on = axis ?? names[names.length - 1];
    return {core: col('Core', new Array(rows).fill(core)),
      rows: names.filter((k) => k !== on).map(pick), column: pick(on)};
  };

  /** The decomposition run end to end, so an assertion lands on the assembled matrix. */
  const run = (smiles: string[], activity: (number | null)[], columns: SarFragmentColumns,
    extra: Partial<SarMatrixParams> = {}): Promise<SarMatrix[]> => {
    const molecules = col('smiles', smiles);
    molecules.semType = DG.SEMTYPE.MOLECULE;
    return runSarMatrix(molecules,
      DG.Column.fromList('double', 'activity', activity) as DG.Column<number>,
      {...e2eParams(false), fragmentColumns: columns, ...extra});
  };

  /** The plan one record gets, the way assembly builds it: the whole cell, or the row key with the
   *  axis point left open. */
  const recordPlan = (columns: SarFragmentColumns, rowCount: number, row = false,
    idx = 0): LinkStages | null => {
    const decomp = decomposeByColumns(columns, null, measured(rowCount)).decomps[0];
    const record = decomp.records[idx];
    const links = decomp.links!;
    const positions = row ? decomp.positions.filter((p) => p !== columns.column.name) : decomp.positions;
    return planLink(record.coreSmiles, record.values, positions, links,
      row ? links.sites[columns.column.name] : [], !row);
  };

  test('fragment columns build the matrix', async () => {
    const {molecules, activity, core, r1, r2} = fragmentTable();
    const matrices = await runSarMatrix(molecules, activity,
      {...e2eParams(false), fragmentColumns: {core, rows: [r1], column: r2}});
    expect(matrices.length, 1, 'no series column, so one matrix');
    const matrix = matrices[0];
    expect(matrix.positions.join(','), 'R2', 'the axis is the column named as the axis');
    expect(matrix.rows.length, 2, 'rows are the distinct core + R1 pairs');
    expect(matrix.columns.length, 2);
    expect(matrix.realCount, 4, 'every compound reaches a cell');
    // Every row draws the core carrying its own substituents, with the axis point left open.
    const keys = matrix.rows.map((row) => row.keySmiles);
    expect(keys.every((k) => k.includes('c1ccc') && k.includes('[*:2]')), true,
      'each row is the core with its R1 attached and the axis still open');
    expect(new Set(keys).size, 2, 'and the rows differ by the R1 that distinguishes them');
  });

  // A connector drawn alone is a bare chain with nothing to recognise it by, so it stays on its core.
  test('a bridging row fragment stays on its core', async () => {
    const matrices = await run(['CC(=O)NCCOCC(=O)Nc1ccccc1', 'CC(=O)NCCCCC(=O)Nc1ccccc1',
      'CC(=O)NCCOCC(=O)Nc1ccncc1'], [5.1, 5.8, 6.4],
    spec('CC(=O)N[*:1]', {
      'Linker': ['[*:1]CCOCC(=O)[*:2]', '[*:1]CCCCC(=O)[*:2]', '[*:1]CCOCC(=O)[*:2]'],
      'E3 ligand': ['[*:2]NC1=CC=CC=C1', '[*:2]NC1=CC=CC=C1', '[*:2]NC1=CC=NC=C1'],
    }), {minCompounds: 1});
    const keys = matrices[0].rows.map((row) => row.keySmiles);
    expect(keys.every((k) => k.includes('C(=O)N') && k.includes('[*:2]')), true,
      'the connector is drawn on its core, with the column attachment left open');
    const virtual = matrices[0].cells.flat().filter((c) => c.kind === 'virtual');
    expect(virtual.length > 0, true, 'the fixture must leave something to propose');
    expect(virtual.every((c) => c.smiles !== null && !c.smiles.includes('[*:')), true,
      'a chained decomposition still assembles its proposals whole');
  });

  // R-Group Analysis writes its core as a molblock: empty title line, coordinates differing per compound.
  test('a molblock core keys as one scaffold', async () => {
    const {r1, r2} = fragmentTable();
    const mol = chemCommonRdKit.getRdKitModule().get_mol('[*:1]c1ccc([*:2])cc1');
    mol.set_new_coords();
    const block = mol.get_molblock();
    mol.delete();
    expect(block.startsWith('\n'), true, 'the fixture must carry the empty title line');
    const core = col('Core', new Array(4).fill(block));
    const {decomps} = decomposeByColumns({core, rows: [r1], column: r2}, null, measured(4));
    expect(new Set(decomps[0].records.map((r) => r.coreSmiles)).size, 1, 'one scaffold keys as one');
    expect(decomps[0].links!.fills['R2'].join(','), '2', 'and still says what each fragment fills');
    expect(recordPlan({core, rows: [r1], column: r2}, 4)!.length, 1, 'so its fragments recombine');
  });

  // The linker skips an attachment the core does not carry and returns what is left, so a core that
  // does not fit its fragments yields a real, parseable, WRONG molecule rather than nothing.
  test('attachments come from the fragments', async () => {
    const {core, r1, r2} = fragmentTable();
    const renamed = col('E3 ligand', r2.toList());
    const named = col('E3 ligand', ['VHL', 'CRBN', 'VHL', 'CRBN']);
    const fills = decomposeByColumns({core, rows: [r1], column: r2}, null, measured(4)).decomps[0].links!.fills;
    expect(fills['R1'].join(',') + '/' + fills['R2'].join(','), '1/2', 'read off the fragments');
    expect(decomposeByColumns({core, rows: [r1], column: renamed}, null, measured(4)).decomps[0].links!
      .fills['E3 ligand'].join(','), '2', 'the column name is not what says it');
    expect(recordPlan({core, rows: [r1], column: named}, 4), null,
      'and a fragment carrying no attachment is never recombined');
  });

  // Two positions fold into the row identity and the third runs across the top.
  test('a three-point core folds two and enumerates one', async () => {
    const matrices = await run(['CNc1cc(OC)nc(-c2ccccc2)n1', 'CCNc1cc(OC)nc(-c2ccccc2)n1',
      'CNc1cc(Cl)nc(-c2ccccc2)n1'], [5.1, 5.8, 6.4],
    spec('[*:1]c1nc([*:2])nc([*:3])c1', {
      R1: ['CN[*:1]', 'CCN[*:1]', 'CN[*:1]'], R2: '[*:2]c1ccccc1',
      R3: ['[*:3]OC', '[*:3]OC', '[*:3]Cl'],
    }), {minCompounds: 1});
    const matrix = matrices[0];
    expect(matrix.positions.join(','), 'R3', 'the last attachment runs across the top');
    expect(`${matrix.rows.length}x${matrix.columns.length}`, '2x2', 'R1+R2 down the side, R3 across');
    expect(matrix.rows.every((row) => row.keySmiles.includes('[*:3]') &&
      !row.keySmiles.includes('[*:1]') && !row.keySmiles.includes('[*:2]')), true,
    'a row is the folded fragments joined onto the core, with only the axis left open');
    const virtual = matrix.cells.flat().filter((c) => c.kind === 'virtual');
    expect(virtual.length, 1, 'the one unmade R1 x R3 combination');
    expect(virtual[0].smiles !== null && !virtual[0].smiles.includes('[*:'), true,
      'and it assembles whole, with every attachment filled');
  });

  // Two R-group runs on one table give `R2` and `R2_1`, both carrying `[*:2]`, and a macrocycle closes
  // between two fragments rather than onto the core. Neither is a walk the linker can make.
  test('a plan the linker cannot walk is refused', async () => {
    const claimed = spec('[*:1]c1ccc([*:2])cc1', {
      R2_1: ['[*:2]CC[*:7]', '[*:2]CCC[*:7]', '[*:2]CC[*:7]', '[*:2]CCC[*:7]'],
      R2: ['F[*:2]', 'F[*:2]', 'Cl[*:2]', 'Cl[*:2]'],
    });
    expect(recordPlan(claimed, 4, true), null, 'the axis point is not the row plan to spend');
    expect(recordPlan(claimed, 4), null, 'and two fragments cannot both fill it');
    expect(recordPlan(spec('[*:1]c1ccc([*:2])cc1', {
      R1: ['[*:1]C(=O)N[*:3]', '[*:1]CCN[*:3]', '[*:1]C(=O)N[*:3]', '[*:1]CCN[*:3]'],
      R2: ['[*:2]CCO[*:3]', '[*:2]CCO[*:3]', '[*:2]CCCO[*:3]', '[*:2]CCCO[*:3]'],
    }), 4), null, 'the bond a macrocycle closes is one the linker never makes');
  });

  // Text RDKit cannot sanitize re-parses as SMARTS on the way back, and a piece carrying one point
  // twice closes on itself — both join into something that reads as a real proposed compound.
  test('unreadable input is never recombined', async () => {
    const smarts = spec('[*:1]c1ccc([*:2])cc1', {
      R1: ['[*:1][C,N]CC', '[*:1][C,N]CC', 'CC[*:1]', 'CC[*:1]'],
      R2: ['F[*:2]', 'Cl[*:2]', 'F[*:2]', 'Cl[*:2]'],
    });
    expect(recordPlan(smarts, 4), null, 'a SMARTS atom list carries `[*:1]` but is not a structure');
    expect(recordPlan(smarts, 4, false, 2) !== null, true, 'only that record is refused');
    expect(recordPlan(spec('O=C([*:1])N([*:2])C(=O)[*:1]',
      {R1: '[*:1]CCCCCC[*:1]', R2: ['[*:2]C', '[*:2]CC', '[*:2]C', '[*:2]CC']}), 4, true), null,
    'nor is a core carrying [*:1] twice');
    // Component names on the axis are a supported layout, and must not cost the rows their plan.
    const labels = spec('[*:1]c1ccc([*:2])cc1', {'R1': ['C[*:1]', 'CC[*:1]', 'C[*:1]', 'CC[*:1]'],
      'E3 ligand': ['VHL', 'VHL', 'CRBN', 'CRBN']});
    expect(recordPlan(labels, 4, true)!.length, 1, 'R1 still meets the core whatever the axis holds');
    expect(recordPlan(labels, 4), null, 'but a name is not a fragment, so no cell can be built');
  });

  // Half an R-group table is blank and the join already erases the point a blank fills, so a blank
  // must not cost the series its plan. A blank CONNECTOR is different: erasing its point strands
  // whatever the next pass would have joined there.
  test('a blank substituent is hydrogen, a blank connector is not', async () => {
    const stages = recordPlan(spec('[*:1]c1cc([*:2])cc([*:3])c1', {
      R1: ['C[*:1]', 'C[*:1]', 'CC[*:1]', 'CC[*:1]'], R2: ['CO[*:2]', '', 'CO[*:2]', ''],
      R3: ['Cl[*:3]', 'F[*:3]', 'Cl[*:3]', 'F[*:3]'],
    }), 4, false, 1)!;
    expect(stages.length, 1, 'the substituted rows still say where every position attaches');
    expect(Object.keys(stages[0]).sort().join(','), 'R1,R2,R3', 'including the one that is blank');
    const connector = spec('[*:1]NC(=O)c1ccccc1', {
      R1: ['[*:1]CC([*:2])C', '', '[*:1]CC([*:2])C', ''],
      R2: ['[*:2]c1ccccc1', '[*:2]C1CC1', '[*:2]C1CC1', '[*:2]c1ccccc1'],
    });
    expect(recordPlan(connector, 4, false, 1), null, 'R2 has nowhere to go once R1 is gone');
    expect(recordPlan(connector, 4, false, 0) !== null, true, 'rows that carry it still assemble');
  });

  // The unsubstituted parent is the reference the substituted columns are read against. Asserted
  // through the matrix, not the records: it can survive decomposition and still never reach a cell.
  test('a blank on the column axis is the parent, not a discard', async () => {
    const columns = spec('[*:1]c1ccc([*:2])cc1', {R1: ['C[*:1]', 'CC[*:1]', 'C[*:1]', 'CC[*:1]'],
      R2: ['', '', 'Cl[*:2]', 'Cl[*:2]']});
    expect(decomposeByColumns(columns, null, measured(4)).decomps[0].records.length, 4, 'none dropped');
    const matrix = (await run(['Cc1ccccc1', 'CCc1ccccc1', 'Cc1ccc(Cl)cc1', 'CCc1ccc(Cl)cc1'],
      [5.2, 5.6, 6.1, 6.5], columns))[0];
    expect(matrix.columns.length, 2, 'the parent takes a column beside the chloro one');
    expect(matrix.columns.some((c) => c.substSmiles === ''), true, 'and it is the blank fragment');
    expect(matrix.realCount, 4, 'every measured compound reaches a cell');
  });

  // Warhead - linker - spacer - cap: the spacer never touches the core, so the row key folds one
  // link at a time. With the core at one end instead, the COLUMN fragment goes on first.
  test('a chain folds transitively, from whichever end the core sits', async () => {
    const chain = spec('CC(=O)N[*:1]', {Linker: '[*:1]CCOCC[*:2]',
      Spacer: ['[*:2]NCC[*:3]', '[*:2]NCCC[*:3]', '[*:2]NCC[*:3]'],
      Cap: ['[*:3]C', '[*:3]C', '[*:3]CC']});
    const rows = recordPlan(chain, 3, true)!;
    expect(rows.length, 2, 'two links means two passes');
    expect(`${rows[0]['Linker']}/${rows[1]['Spacer']}`, '1/2', 'each meets what the last exposed');
    expect(recordPlan(chain, 3)!.length, 3, 'a whole cell takes one more pass for the cap');
    const fromLigand = recordPlan(spec('[*:2]NC1=CC=CC=C1', {
      Warhead: ['CC(=O)N[*:1]', 'CCC(=O)N[*:1]', 'CC(=O)N[*:1]', 'CCC(=O)N[*:1]'],
      Linker: ['[*:1]CCOCC(=O)[*:2]', '[*:1]CCOCC(=O)[*:2]', '[*:1]CCCCC(=O)[*:2]',
        '[*:1]CCCCC(=O)[*:2]'],
    }), 4)!;
    expect(`${fromLigand[0]['Linker']}/${fromLigand[1]['Warhead']}`, '2/1',
      'the column fragment meets the core first, the row fragment what it exposed');
  });

  test('each core is its own series', async () => {
    const {r1, r2} = fragmentTable();
    const twoCores = col('Core', ['[*:1]c1ccc([*:2])cc1', '[*:1]c1ccc([*:2])cc1',
      '[*:1]C1CCC([*:2])CC1', '[*:1]C1CCC([*:2])CC1']);
    expect(decomposeByColumns({core: twoCores, rows: [r1], column: r2}, null, measured(4)).clusters.length, 2,
      'two scaffolds make two matrices, not one of twice the height');
    const pooled = decomposeByColumns({core: twoCores, rows: [], column: r2}, null, measured(4));
    expect(pooled.clusters.length, 1, 'with nothing on the row axis the cores become the rows');
    expect(new Set(pooled.decomps[0].records.map((r) => r.coreSmiles)).size, 2, 'and both reach it');
  });

  /** A furan (terminal) and a pyrimidine (connector) in one R1 column, plus the des-substituted
   *  parent, which is what joins the two topologies into one comparable block. */
  const mixed = (): {smiles: string[], activity: number[], columns: SarFragmentColumns} => {
    const pyr = '[*:1]c1ncc([*:2])cn1';
    return {
      smiles: ['c1ccc(-c2ccoc2)cc1', 'c1ccc(-c2ncccn2)cc1', 'Cc1cnc(-c2ccccc2)nc1',
        'Clc1cnc(-c2ccccc2)nc1', 'COc1cnc(-c2ccccc2)nc1'],
      activity: [5.1, 5.4, 6.0, 6.3, 6.6],
      columns: spec('[*:1]c1ccccc1', {R1: ['[*:1]c1ccoc1', pyr, pyr, pyr, pyr],
        R2: ['', '', 'C[*:2]', 'Cl[*:2]', 'CO[*:2]']}),
    };
  };

  // Cores of different shapes in one table: three carry every attachment themselves, the fourth
  // carries two and lets R3 hang off R2 and R4 off R3. Each shape is read per compound, so the
  // chained one folds transitively and still leaves the axis point open.
  test('cores of different shapes each fold their own way', async () => {
    const direct = ['[*:1]C1CCC([*:2])C([*:3])C1[*:4]', '[*:1]C1CC([*:2])C([*:3])C1[*:4]',
      '[*:1]C1CCCC([*:2])C([*:3])C1[*:4]'];
    const chained = '[*:1]c1ccc([*:2])cc1';
    const cores: string[] = [];
    const r1: string[] = [];
    const r2: string[] = [];
    const r3: string[] = [];
    const r4: string[] = [];
    const mols: string[] = [];
    const block = (core: string, second: string, third: string): void => {
      for (const first of ['C[*:1]', 'CC[*:1]']) {
        for (const last of ['[*:4]Cl', '[*:4]F']) {
          cores.push(core); r1.push(first); r2.push(second); r3.push(third); r4.push(last);
          mols.push('C'.repeat(1 + mols.length % 6) + 'O'.repeat(1 + Math.floor(mols.length / 6)));
        }
      }
    };
    for (const core of direct)
      block(core, '[*:2]C', '[*:3]OC');
    // The same R3 column holds a terminal fragment on the cores above and a bridge on this one.
    block(chained, '[*:2]c1ccc([*:3])cc1', '[*:3]CC[*:4]');
    // One combination left unmade, so its prediction has to be built through the whole chain.
    for (const arr of [cores, r1, r2, r3, r4, mols])
      arr.pop();

    const matrices = await run(mols, mols.map((_v, i) => 5 + i * 0.1), {
      core: col('Core', cores), rows: [col('R1', r1), col('R2', r2), col('R3', r3)],
      column: col('R4', r4)}, {minCompounds: 1});

    expect(matrices.length, 4, 'each core is its own series, whatever shape it is');
    const chain = matrices.find((m) => m.rows[0].coreSmiles.includes('c1'))!;
    const flat = matrices.filter((m) => m !== chain);
    expect(chain.rows.every((row) => row.keySmiles.includes('[*:4]') &&
      !row.keySmiles.includes('[*:3]')), true,
    'the chained row folds R2 and R3 onto the core and leaves only the axis open');
    expect(flat.every((m) => m.rows.every((row) => row.keySmiles.includes('C1'))), true,
      'a row of terminal fragments still carries its own core');
    const predicted = chain.cells.flat().filter((c) => c.kind === 'virtual');
    expect(predicted.length > 0, true, 'the unmade combination is offered');
    expect(predicted.every((c) => c.smiles !== null && !c.smiles.includes('[*:')), true,
      'and assembles whole, through both connectors');
    // No point on the chained core is where the columns hang — the axis is two links out — so
    // nothing may be marked as one.
    expect(matrixCore(chain).split('\n').some((l) => l.slice(31, 34).trim() === 'R'), false,
      'an axis the core does not carry is not drawn on it');
  });

  test('a column mixing shapes keeps both kinds of row', async () => {
    const {smiles, activity, columns} = mixed();
    const matrix = (await run(smiles, activity, columns))[0];
    expect(matrix.rows.length, 2, 'the furan and the pyrimidine are both rows');
    expect(matrix.realCount, 5, 'no compound is dropped for disagreeing with the other topology');
    expect(matrix.columns.some((c) => c.substSmiles === ''), true, 'the parent keeps its own column');
    expect(new Set(matrix.rows.map((row) => row.keySmiles)).size, 2,
      'each row is drawn from its own fragments, not the core they share');
  });

  // The PROTAC shape: the core is the bridge, a terminal arm hangs off each end, and the axis is one
  // of the arms. The row is the bridge carrying the other arm, with the axis end still open.
  test('a connector core carries an arm on each side', async () => {
    const linkers = ['[*:1]CCOCC[*:2]', '[*:1]CCCC[*:2]'];
    const warheads = ['[*:1]c1ccc(Cl)cc1', '[*:1]c1ccc(F)cc1'];
    const ligands = ['[*:2]N1CCCCC1', '[*:2]C1CCC(=O)NC1=O'];
    const cores: string[] = [];
    const wh: string[] = [];
    const e3: string[] = [];
    const mols: string[] = [];
    for (const linker of linkers) {
      for (const warhead of warheads) {
        for (const ligand of ligands) {
          cores.push(linker); wh.push(warhead); e3.push(ligand);
          mols.push('C'.repeat(1 + mols.length) + 'O');
        }
      }
    }
    // One pairing left unmade, so it has to be predicted across the whole bridge.
    for (const arr of [cores, wh, e3, mols])
      arr.pop();

    const matrices = await run(mols, mols.map((_v, i) => 7 + i * 0.1), {
      core: col('Linker', cores), rows: [col('Warhead', wh)],
      column: col('E3 ligand', e3)}, {minCompounds: 1});

    expect(matrices.length, 2, 'each linker is its own series');
    expect(matrices.every((m) => m.rows.every((row) => row.keySmiles.includes('[*:2]') &&
      !row.keySmiles.includes('[*:1]'))), true,
    'a row is the linker carrying its warhead, with the ligand end still open');
    const predicted = matrices.flatMap((m) => m.cells.flat()).filter((c) => c.kind === 'virtual');
    expect(predicted.length > 0, true, 'the unmade pairing is offered');
    expect(predicted.every((c) => c.smiles !== null && !c.smiles.includes('[*:')), true,
      'and assembles across the bridge with both ends filled');
  });

  // The furan carries no second attachment, so predicting there asks for a compound the linker would
  // build as something else — the substituent is dropped and an already-measured molecule returns.
  test('a cell with no attachment point is greyed, not predicted', async () => {
    const {smiles, activity, columns} = mixed();
    const matrix = (await run(smiles, activity, columns))[0];
    const furan = matrix.rows.findIndex((row) => row.keySmiles.includes('o1') ||
      row.keySmiles.includes('co'));
    expect(furan >= 0, true, 'the furan row must survive');
    matrix.columns.forEach((column, ci) => expect(matrix.cells[furan][ci].kind,
      column.substSmiles === '' ? 'real' : 'impossible', `furan x "${column.substSmiles}"`));
    const impossible = matrix.cells.flat().filter((cell) => cell.kind === 'impossible');
    expect(impossible.length, 3, 'exactly the three substituted furan combinations');
    expect(impossible.every((c) => c.value === null && c.smiles === null), true,
      'carrying neither a potency nor a structure');
    expect(matrix.cells.flat().every((c) => c.kind !== 'virtual' ||
      (c.smiles !== null && !c.smiles.includes('[*:'))), true, 'and every proposal assembles whole');
  });

  // Where the rule contradicts a compound that was measured, it is the rule misreading that row.
  test('a measured compound is never called impossible', async () => {
    const matrices = await run(['CC(CNC(=O)c1ccccc1)c1ccccc1', 'C1CC1CC(=O)Nc1ccccc1',
      'CC(CNC(=O)c1ccccc1)C1CC1', 'CC(=O)Nc1ccccc1'], [5.5, 5.9, 6.2, 6.4],
    spec('[*:1]NC(=O)c1ccccc1', {R1: ['[*:1]CC([*:2])C', '', '[*:1]CC([*:2])C', ''],
      R2: ['[*:2]c1ccccc1', '[*:2]C1CC1', '[*:2]C1CC1', '[*:2]c1ccccc1']}), {minCompounds: 1});
    expect(matrices.every((m) => m.cells.flat().every((c) => c.kind !== 'impossible')), true,
      'every cell of that row holds a compound, so the row is read wrongly rather than impossible');
  });

  // A core point outside the picked columns is the series' problem, already reported; letting it
  // condemn cells would grey every hole in a table where one R column was simply left out.
  test('a point no column fills is not an impossible cell', async () => {
    const core = '[*:1]c1cc([*:2])cc([*:3])c1';
    const links = decomposeByColumns(spec(core, {R1: ['C[*:1]', 'CC[*:1]', 'C[*:1]', 'CC[*:1]'],
      R2: ['Cl[*:2]', 'Cl[*:2]', 'F[*:2]', 'F[*:2]']}), null, measured(4)).decomps[0].links!;
    expect(cellPossible(core, {R1: 'C[*:1]', R2: 'Cl[*:2]'}, ['R2', 'R1'], links), true,
      '[*:3] is outside the decomposition, so it proves nothing');
  });

  test('every R-group notation reads as [*:n]', async () => {
    for (const written of ['C[*:1]', '[1*]C', '[*1]C', '[R1]C', '[R:1]C', '[2*:1]C', '*C |$_R1;$|'])
      expect(standardizeFragment(written).value, 'C[*:1]', written);
    const rdkit = chemCommonRdKit.getRdKitModule();
    const mol = rdkit.get_mol('[1*]c1ccc([2*])cc1');
    mol.set_new_coords();
    // R-Groups Analysis writes its core this way.
    const block = mol.get_molblock().replace('ISO', 'RGP');
    mol.delete();
    expect(standardizeFragment(block).value, standardizeFragment('[*:1]c1ccc([*:2])cc1').value, 'an RGP molblock');
    expect(standardizeFragment('VHL').structure, false, 'a name stays a name');
  });

  test('an unnumbered attachment point takes its column\'s number', async () => {
    const core = '[*:1]c1ccc([*:2])cc1';
    const {decomps} = decomposeByColumns(spec(core, {R1: ['C[*:1]', 'CC[*:1]', 'C[*:1]', 'CC[*:1]'],
      R2: ['*F', '*F', '*Cl', '*Cl']}), null, measured(4));
    expect(decomps[0].records.map((r) => r.values['R2']).join(','), 'F[*:2],F[*:2],Cl[*:2],Cl[*:2]');

    // Nothing numbers "Aryl", so nothing about it can be claimed impossible.
    const matrix = (await run(['Cc1ccccc1', 'CCc1ccccc1', 'CCc1ccc(F)cc1'], [5.1, 5.6, 6.2],
      spec(core, {R1: ['C[*:1]', 'CC[*:1]', 'CC[*:1]'], Aryl: ['', '', '*F']})))[0];
    expect(matrix.cells.flat().some((c) => c.kind === 'impossible'), false, 'no cell is struck through');
    expect(matrix.cells.flat().some((c) => c.kind === 'virtual'), true, 'the hole is predicted instead');
  });

  test('a repeated compound keeps its measured copy', async () => {
    const matrix = (await run(['Cc1ccc(F)cc1', 'CCc1ccc(F)cc1', 'CCc1ccc(F)cc1', 'Cc1ccc(Cl)cc1', 'CCc1ccc(Cl)cc1'],
      [5.1, null, 5.8, 6.4, 7.0], spec('[*:1]c1ccc([*:2])cc1', {
        R1: ['C[*:1]', 'CC[*:1]', 'CC[*:1]', 'C[*:1]', 'CC[*:1]'],
        R2: ['F[*:2]', 'F[*:2]', 'F[*:2]', 'Cl[*:2]', 'Cl[*:2]']})))[0];
    const cells = matrix.cells.flat();
    expect(cells.some((c) => c.kind === 'real' && c.molIdx === 2), true, 'the measured copy takes the cell');
    expect(cells.some((c) => c.molIdx === 1), false, 'the unmeasured copy does not');
    expect(matrix.realCount, 4);
  });

  test('R-groups that are not the table\'s molecules are caught', async () => {
    const {molecules, core, r1, r2} = fragmentTable();
    const {decomps} = decomposeByColumns({core, rows: [r1], column: r2}, null, measured(4));
    const matching = await checkAgainstMolecules(decomps, molecules.toList());
    expect(`${matching.mismatched}/${matching.checked}`, '0/4', 'each compound rebuilds from its own R-groups');
    const shifted = await checkAgainstMolecules(decomps, [...molecules.toList()].reverse());
    expect(`${shifted.mismatched}/${shifted.checked}`, '4/4', 'rows out of step with the molecules are all reported');
  });

  test('a bridge written in both of its columns closes one ring', async () => {
    const bridge = 'C(C[*:101])[*:100]';
    const columns = spec('C1CC(C[*:100])(N[*:1])C([*:101])N1[*:2]',
      {R1: ['CC(=O)[*:1]', 'CC(=O)[*:1]'], R100: bridge, R101: bridge, R2: ['C[*:2]', 'CC[*:2]']});
    const {decomps} = decomposeByColumns(columns, null, measured(2));
    const record = decomps[0].records[0];
    const plan = planLink(record.coreSmiles, record.values, decomps[0].positions, decomps[0].links!, [], true);
    expect(plan !== null, true, 'the bridge is one piece, not two claims on one point');
    const [built] = await linkStaged([record.coreSmiles], [plan], (_i, p) => record.values[p] ?? '');
    const canonical = (smiles: string): string => {
      const mol = chemCommonRdKit.getRdKitModule().get_mol(smiles);
      const out = mol.get_smiles();
      mol.delete();
      return out;
    };
    expect(canonical(built!), canonical('CC(=O)NC12CCN(C)C1CCC2'), 'the bicyclic compound it describes');
  });

  /** SPGI's compounds that carry a decomposition, with Average Mass as the value. */
  let spgi: DG.DataFrame | null = null;
  const spgiWithCore = async (): Promise<DG.DataFrame> => {
    if (spgi === null) {
      const df = await grok.dapi.files.readCsv('System:DemoFiles/chem/SPGI.csv');
      const core = df.col('Core')!;
      spgi = df.clone(DG.BitSet.create(df.rowCount, (i) => !core.isNone(i) && core.get(i) !== ''));
    }
    return spgi;
  };
  const spgiColumns = (df: DG.DataFrame): SarFragmentColumns => ({core: df.col('Core')!,
    rows: ['R1', 'R3', 'R100', 'R101'].map((name) => df.col(name)!), column: df.col('R2')!});

  test('SPGI: the dialog offers the R-group columns and picks R2', async () => {
    const df = await spgiWithCore();
    const offered = df.columns.toList().filter(holdsFragments).map((c) => c.name);
    expect(offered.join(','), 'Core,R1,R2,R3,R100,R101', 'only columns holding attachment points');
    expect(defaultAxis(['R1', 'R2', 'R3', 'R100', 'R101'].map((name) => df.col(name)!))?.name, 'R2',
      'the R-group every compound has, not the bridge');
  });

  test('SPGI: a series too large for one matrix is named, with what it shows and leaves out', async () => {
    const df = await spgiWithCore();
    const activities = new Float32Array(df.rowCount).map((_v, i) => df.col('Average Mass')!.get(i));
    expect(decomposeByColumns(spgiColumns(df), null, activities).cuts.length, 0, 'R2 across fits every series');
    const byR1: SarFragmentColumns = {core: df.col('Core')!,
      rows: ['R2', 'R3', 'R100', 'R101'].map((name) => df.col(name)!), column: df.col('R1')!};
    const {cuts} = decomposeByColumns(byR1, null, activities);
    expect(cuts.length, 1, 'with R1 across, only the triazoles exceed one matrix');
    const cut = cuts[0];
    expect(cut.keptRows === cut.rows && cut.keptColumns < cut.columns, true, 'the columns are what is cut');
    expect(cut.keptRows * cut.keptColumns <= 250000 && cut.leftOut > 0, true, JSON.stringify(cut));
    const message = cutWarning(cuts, () => 'Series 4')!;
    const n = (value: number): string => value.toLocaleString('en-US');
    for (const part of [`Series 4 is too large for one matrix (${n(cut.rows)} × ${n(cut.columns)} R-group`,
      `shows only the ${n(cut.keptColumns)} columns with the most compounds`,
      `leaves out ${n(cut.leftOut)} of its ${n(cut.compounds)} compounds`, 'fewer values for Matrix columns'])
      expect(message.includes(part), true, `"${part}" in: ${message}`);
  });

  test('SPGI: every compound rebuilds from its core and R-groups', async () => {
    const df = await spgiWithCore();
    const activities = new Float32Array(df.rowCount).map((_v, i) => df.col('Average Mass')!.get(i));
    const {decomps} = decomposeByColumns(spgiColumns(df), null, activities);
    expect(decomps.length, 6, 'one series per core');
    const {checked, mismatched} = await checkAgainstMolecules(decomps, df.col('Structure')!.toList());
    expect(checked > 250, true, `${checked} compounds rebuilt`);
    expect(mismatched, 0, 'and each is the compound on its row');
  }, {timeout: 120000});

  test('SPGI: rows and proposals come out whole', async () => {
    const df = await spgiWithCore();
    const matrices = await runSarMatrix(df.col('Structure')!, df.col('Average Mass') as DG.Column<number>,
      {...e2eParams(false), minCompounds: 3, fragmentColumns: spgiColumns(df)});
    expect(matrices.length >= 5, true, `${matrices.length} series`);
    // The two compounds whose R1 and R2 are one ring cannot vary R2, so only their rows stay open.
    const ringRow = (row: SarMatrixRow): boolean => (row.foldedValues['R1'] ?? '').includes('[*:2]');
    for (const matrix of matrices) {
      for (const row of matrix.rows.filter((r) => !ringRow(r))) {
        expect(row.keySmiles.replace('[*:2]', '').includes('[*:'), false,
          `${matrix.label}: row ${row.keySmiles} leaves only R2 open`);
      }
      matrix.rows.forEach((row, ri) => {
        if (!ringRow(row))
          expect(matrix.cells[ri].some((c) => c.kind === 'impossible'), false, `${matrix.label} ${row.label}`);
      });
    }
    const predicted = matrices.flatMap((m) => m.cells.flat()).filter((c) => c.kind === 'virtual');
    const drawn = predicted.filter((c) => c.smiles !== null);
    expect(drawn.every((c) => !c.smiles!.includes('[*:')), true, 'no proposal keeps an attachment point');
    expect(drawn.length >= predicted.length * 0.95, true, `${drawn.length} of ${predicted.length} proposals drawn`);
    const bridged = matrices.find((m) => m.rows[0].coreSmiles.includes('[*:100]'));
    expect(bridged !== undefined && bridged.rows.length >= 2, true, 'the bridged pyrrolidines make a series');
    expect(bridged!.rows.every((row) => !row.keySmiles.includes('[*:100]') && !row.keySmiles.includes('[*:101]')),
      true, 'their bridge is closed in every row');
  }, {timeout: 300000});
  // ---- One additive fit over every role column ---------------------------------------------------

  const CORE_EFFECT: {[value: string]: number} = {C1: 0.3, C2: 0.0, C3: -0.3};
  const WARHEAD_EFFECT: {[value: string]: number} = {W1: 0.2, W2: 0.0, W3: -0.2};
  const LIGAND_EFFECT: {[value: string]: number} = {L1: 0.1, L2: -0.1};

  /** One compound of a role design: its value in each role, its activity and its source row. */
  type DesignRow = [string[], number, number];
  const design = (names: string[], rows: DesignRow[]): RoleDesign => ({names,
    values: names.map((_name, r) => rows.map(([values]) => values[r])),
    activity: rows.map(([, activity]) => activity), molIdx: rows.map(([, , mol]) => mol),
    minSupport: 3, higherIsBetter: true});

  /**
   * 3 cores x 3 warheads x 2 ligands, `y = 6.0 + core + warhead + ligand`, no noise, molIdx 0..17.
   * Balanced and count-centred, so the least-squares solution is exactly those effects, and the
   * count-prior of 1 then shrinks each by `n/(n+1)`: 6/7 for a core or warhead, 9/10 for a ligand.
   */
  const roleFixture = (drop: number[] = [], ligandNames = ['L1', 'L2']): RoleDesign => {
    const rows: DesignRow[] = [];
    const ligands = Object.keys(LIGAND_EFFECT);
    let k = 0;
    for (const core of Object.keys(CORE_EFFECT)) {
      for (const warhead of Object.keys(WARHEAD_EFFECT)) {
        for (let li = 0; li < ligands.length; li++, k++) {
          if (!drop.includes(k)) {
            rows.push([[core, warhead, ligandNames[li]],
              6.0 + CORE_EFFECT[core] + WARHEAD_EFFECT[warhead] + LIGAND_EFFECT[ligands[li]], k]);
          }
        }
      }
    }
    return design(['Core', 'Warhead', 'Ligand'], rows);
  };

  const roleOf = (fit: RoleFit, name: string): RoleSummary => fit.roles.find((r) => r.name === name)!;
  const coefOf = (fit: RoleFit, name: string, value: string): number =>
    roleOf(fit, name).levels.find((l) => l.value === value)!.coef;
  const meanOf = (design: RoleDesign): number =>
    design.activity.reduce((sum, v) => sum + v, 0) / design.activity.length;

  test('fitRoleEffects recovers a balanced design', async () => {
    const raw = fitRoleEffects(roleFixture());
    expect(raw !== null, true, 'three roles with three, three and two levels are rankable');
    const fit = raw!;
    expect(fit.converged, true);
    expect(fit.compounds, 18);
    expect(fit.dropped, 0);
    // Raw group means answer 0.3 exactly and centring inside a core block answers 0.15; neither is the
    // penalised least-squares solution that reproduces all 18 values.
    for (const [value, effect] of Object.entries(CORE_EFFECT))
      expectFloat(coefOf(fit, 'Core', value), effect * 6 / 7, 1e-5, `Core ${value}`);
    for (const [value, effect] of Object.entries(WARHEAD_EFFECT))
      expectFloat(coefOf(fit, 'Warhead', value), effect * 6 / 7, 1e-5, `Warhead ${value}`);
    for (const [value, effect] of Object.entries(LIGAND_EFFECT))
      expectFloat(coefOf(fit, 'Ligand', value), effect * 9 / 10, 1e-5, `Ligand ${value}`);
    expectFloat(fit.mean, 6.0, 1e-6, 'the library average every offset is read against');
    expect(fit.roles.map((r) => r.name).join(','), 'Core,Warhead,Ligand',
      'roles come back ordered by noise-corrected spread, widest first');
    expect(roleOf(fit, 'Core').repeat! > 0.99, true,
      `a noiseless design has to refit to the same offsets on each half, got ${roleOf(fit, 'Core').repeat}`);
    // Two points correlate at exactly ±1 whatever they are, so a two-value role reporting r 1.00 would
    // wear the strongest repeatability chip on the screen on no evidence at all.
    expect(roleOf(fit, 'Ligand').repeat, null, 'two comparable levels are not a repeatability check');
  });

  /** The same penalised backfit written the slow way and run far past the module's own tolerance: the
   *  fixed point the shipped fit has to land on, however many sweeps that takes. */
  const longRunFit = (design: RoleDesign, passes: number): {[value: string]: number} => {
    const levels = design.values.map((column) => Array.from(new Set(column)));
    const codes = design.values.map((column, r) => column.map((v) => levels[r].indexOf(v)));
    const counts = levels.map((values, r) => {
      const n = values.map(() => 0);
      for (const code of codes[r])
        n[code]++;
      return n;
    });
    const theta = levels.map((values) => values.map(() => 0));
    let mu = meanOf(design);
    for (let pass = 0; pass < passes; pass++) {
      for (let r = 0; r < levels.length; r++) {
        const acc = levels[r].map(() => 0);
        for (let k = 0; k < design.activity.length; k++) {
          let others = mu;
          for (let q = 0; q < levels.length; q++) {
            if (q !== r)
              others += theta[q][codes[q][k]];
          }
          acc[codes[r][k]] += design.activity[k] - others;
        }
        const next = acc.map((a, v) => a / (counts[r][v] + 1));
        const total = counts[r].reduce((sum, c) => sum + c, 0);
        const shift = next.reduce((sum, t, v) => sum + counts[r][v] * t, 0) / total;
        mu += shift;
        for (let v = 0; v < next.length; v++)
          theta[r][v] = next[v] - shift;
      }
    }
    const out: {[value: string]: number} = {};
    levels.forEach((values, r) => values.forEach((v, i) => out[v] = theta[r][i]));
    return out;
  };

  test('fitRoleEffects recovers an unbalanced design', async () => {
    const design = roleFixture([4, 9, 15]);
    const fit = fitRoleEffects(design)!;
    expect(fit.converged, true, 'a fit that stops short suppresses the whole breakdown silently');
    // The prior leaves each role off centre by a constant on every pass. Measured before that shift is
    // taken back out, the movement plateaus at the constant and no sweep cap ever clears it.
    expect(fit.sweeps < ROLE_FIT_MAX_SWEEPS / 4, true, `${fit.sweeps} sweeps`);
    const reference = longRunFit(design, 4000);
    const mean = meanOf(design);
    let offFixedPoint = 0;
    let offOnePass = 0;
    for (let r = 0; r < design.names.length; r++) {
      for (const level of roleOf(fit, design.names[r]).levels) {
        let sum = 0;
        let n = 0;
        for (let k = 0; k < design.activity.length; k++) {
          if (design.values[r][k] === level.value) {
            sum += design.activity[k];
            n++;
          }
        }
        offFixedPoint = Math.max(offFixedPoint, Math.abs(level.coef - reference[level.value]));
        offOnePass = Math.max(offOnePass, Math.abs(level.coef - (sum / n - mean)));
      }
    }
    expect(offFixedPoint < 1e-5, true, `off the long-run fixed point by ${offFixedPoint}`);
    // With holes a level's own group mean carries whatever it happened to be paired with: W2 is
    // measured mostly on the strongest core and reads +0.105 that way against +0.005 adjusted.
    // Averaging the margins in one pass is a different and worse estimator, not a shortcut to this one.
    expect(offOnePass > 0.05, true, `one-pass margins differ by only ${offOnePass}`);
    // Σ n_v·θ_v = 0 per role is what makes an offset a difference from the library average, and what
    // makes a level no training fold saw predict at that average instead of at nothing.
    for (const role of fit.roles) {
      let weighted = 0;
      for (const level of role.levels)
        weighted += level.n * level.coef;
      expectFloat(weighted, 0, 1e-6, `${role.name} is not count-centred`);
    }
  });

  /**
   * Sixteen values of A crossed with three of P, and a B that is A with one pair of compounds
   * exchanged — so A and B are near-aliased and the credit between them only separates slowly.
   */
  const slowFixture = (): RoleDesign => {
    const rows: DesignRow[] = [];
    for (let i = 0; i < 16; i++) {
      for (let p = 0; p < 3; p++) {
        for (let rep = 0; rep < 5; rep++) {
          const k = rows.length;
          rows.push([[`A${i}`, `B${i === 0 && rep === 0 ? 1 : i === 1 && rep === 0 ? 0 : i}`, `P${p}`],
            6.0 + (i - 8) * 0.05 + (p - 1) * 0.05 + ((k % 5) - 2) * 0.01, k]);
        }
      }
    }
    return design(['A', 'B', 'P'], rows);
  };

  test('a near-aliased design is carried to its fixed point', async () => {
    const fit = fitRoleEffects(slowFixture())!;
    // A cap sized for one small two-factor matrix stops here, and a fit that stops short is refused
    // outright — the breakdown disappears rather than reporting anything.
    expect(fit.converged, true);
    expect(fit.sweeps > 50, true, `two near-aliased roles need ${fit.sweeps} passes`);
    // Testing the sweep's movement on the last role alone stops after two passes with A0 at −0.33: the
    // last role can be still while an earlier one is a third of its final value away.
    expectFloat(coefOf(fit, 'A', 'A0'), -0.1921, 1e-3);
  });

  test('levels below the support floor are fitted but not ranked', async () => {
    const design = roleFixture();
    const plain = fitRoleEffects(design)!;
    // Two compounds carrying a fourth warhead a full log unit above everything else, on opposite cores
    // and opposite ligands so the pair adds no core or ligand imbalance of its own.
    design.values[0].push('C1', 'C3');
    design.values[1].push('W4', 'W4');
    design.values[2].push('L1', 'L2');
    design.activity.push(6.0 + 0.3 + 0.1 + 1.0, 6.0 - 0.3 - 0.1 + 1.0);
    design.molIdx.push(18, 19);
    const fit = fitRoleEffects(design)!;
    const warhead = roleOf(fit, 'Warhead');
    expect(warhead.fitted, 4, 'the thin level stays in the fit — its compounds inform the other roles');
    expect(warhead.thin, 1);
    expect(warhead.levels.some((l) => l.value === 'W4'), false, 'but carries no readable offset');
    // Ungated it absorbs its own two rows whole and tops the leaderboard at three times the widest
    // offset any well-sampled warhead has.
    expect(roleOf(fitRoleEffects({...design, minSupport: 2})!, 'Warhead').levels[0].value, 'W4');
    // The other roles only feel the warhead role being re-centred over one more level; the pair's own
    // residual does not reach them.
    let moved = 0;
    for (const name of ['Core', 'Ligand']) {
      for (const level of roleOf(fit, name).levels)
        moved = Math.max(moved, Math.abs(level.coef - coefOf(plain, name, level.value)));
    }
    expect(moved < 0.05, true, `the other roles moved by ${moved}`);
  });

  test('levels are ordered by the activity direction', async () => {
    const design = roleFixture();
    const up = roleOf(fitRoleEffects(design)!, 'Warhead').levels;
    const down = roleOf(fitRoleEffects({...design, higherIsBetter: false})!, 'Warhead').levels;
    expect(up.length > 1, true, 'the role has levels to order');
    expect(up[0].coef > up[up.length - 1].coef, true, 'higher is better puts the largest offset first');
    expect(down[0].value, up[up.length - 1].value,
      'lower is better puts the same role\'s worst-on-the-other-reading level first');
    // The fit itself is in the column's own units either way: only the order changes.
    expect(down.map((l) => l.value).sort().join(), up.map((l) => l.value).sort().join());
    for (const level of down)
      expectFloat(level.coef, up.find((l) => l.value === level.value)!.coef, 1e-9);
  });

  test('a blank role value is a level', async () => {
    const fit = fitRoleEffects(roleFixture([], ['L1', '']))!;
    const blank = roleOf(fit, 'Ligand').levels.find((l) => l.value === '');
    expect(blank !== undefined, true, 'the unsubstituted parent is a value, not a missing one');
    expect(blank!.n, 9, 'dropping its rows would exclude compounds the tab counts as placed');
    expectFloat(blank!.coef, -0.1 * 9 / 10, 1e-5, 'and it is fitted like every other level');
  });

  const PARTNER_EFFECT = [-0.15, -0.05, 0.05, 0.15];

  /** One role of `effects.length` values fully crossed with a four-value partner, plus a three-cycle
   *  residual belonging to neither role so the fit has an error the split rule has to clear. */
  const splitFixture = (effects: number[]): RoleDesign => {
    const rows: DesignRow[] = [];
    for (let a = 0; a < effects.length; a++) {
      for (let b = 0; b < PARTNER_EFFECT.length; b++) {
        const k = rows.length;
        rows.push([[`V${a}`, `P${b}`], 6.0 + effects[a] + PARTNER_EFFECT[b] + ((k % 3) - 1) * 0.06, k]);
      }
    }
    return design(['Component', 'Partner'], rows);
  };

  test('the two-group split fires only when it should', async () => {
    const twoGroups = fitRoleEffects(splitFixture([0.2, 0.2, 0.2, 0.2, -0.3, -0.3, -0.3]))!;
    const split = roleOf(twoGroups, 'Component').split;
    expect(split !== null, true, 'a 0.5 step against a ±0.10 cross-validated error is a partition');
    expect(`${split!.hiCount}|${split!.loCount}`, '4|3');
    expect(split!.gap > twoGroups.cvRmse!, true, `gap ${split!.gap} under the error`);
    expect(split!.hiMean > 0 && split!.loMean < 0, true, 'each side reports its own mean offset');
    expect(`${split!.hiN}|${split!.loN}`, '16|12', 'and the compounds behind it, not its level count');

    const even = fitRoleEffects(splitFixture([-0.06, -0.04, -0.02, 0, 0.02, 0.04, 0.06]))!;
    expect(roleOf(even, 'Component').split, null,
      'evenly spaced values hold no gap the fit can resolve');

    // The widest gap here is 0.37, far over the error, and the split is still refused: one value is not
    // a group, and reporting it as one would call a single outlier half the chemistry.
    const outlier = fitRoleEffects(splitFixture([0.5, 0.02, 0.01, 0, -0.01, -0.02, -0.03]))!;
    expect(roleOf(outlier, 'Component').split, null);
  });

  test('disconnected blocks are not ranked together', async () => {
    const design = roleFixture();
    let k = 18;
    // Two cores x three warheads on one ligand, sharing no value in any role with the block above.
    for (const core of ['X1', 'X2']) {
      for (const warhead of ['Y1', 'Y2', 'Y3']) {
        design.values[0].push(core);
        design.values[1].push(warhead);
        design.values[2].push('Z1');
        design.activity.push(9.0 + (k % 2) * 0.2);
        design.molIdx.push(k++);
      }
    }
    const fit = fitRoleEffects(design)!;
    expect(fit.compounds, 18, 'only the largest connected block is fitted');
    expect(fit.dropped, 6);
    const listed = fit.roles.flatMap((role) => role.levels.map((level) => level.value));
    expect(listed.some((v) => v.startsWith('X') || v.startsWith('Y') || v === 'Z1'), false,
      'an offset from a separate block rests on an unrelated baseline, so it is never listed beside these');
    for (const role of fit.roles)
      expect(role.fitted, role.name === 'Ligand' ? 2 : 3, `${role.name} counts only the kept set`);
  });

  const SOLID_EFFECT = [0.035, 0.0105, -0.0105, -0.035];

  /** Four well-sampled values carrying a real effect, crossed with sixty values seen three times each
   *  carrying none — the shape where a raw count-weighted sd ranks the roles the wrong way round. */
  const spreadFixture = (): RoleDesign => {
    const rows: DesignRow[] = [];
    let seed = 4242;
    const rnd = (): number => {
      seed = (seed * 1103515245 + 12345) & 0x7fffffff;
      return seed / 0x7fffffff - 0.5;
    };
    for (let t = 0; t < 60; t++) {
      for (let c = 0; c < 3; c++)
        rows.push([[`S${(t * 3 + c) % 4}`, `T${t}`], 6.0 + SOLID_EFFECT[(t * 3 + c) % 4] + 0.5 * rnd(), rows.length]);
    }
    return design(['Solid', 'Thin'], rows);
  };

  test('the spread subtracts the estimation noise it is comparing across', async () => {
    const fit = fitRoleEffects(spreadFixture())!;
    const rawSd = (role: RoleSummary): number => {
      let weighted = 0;
      let total = 0;
      for (const level of role.levels) {
        weighted += level.n * level.coef * level.coef;
        total += level.n;
      }
      return Math.sqrt(weighted / total);
    };
    const solid = roleOf(fit, 'Solid');
    const thin = roleOf(fit, 'Thin');
    // Sixty values seen three times each carry no effect here, and still scatter: fitting each from
    // three compounds leaves a σ²·levels/observations term in the raw count-weighted sd that a
    // four-value role does not pay. Raw, that term alone outranks the role with the real effect.
    expect(rawSd(thin) > rawSd(solid), true,
      `raw sds are ${rawSd(thin).toFixed(4)} thin against ${rawSd(solid).toFixed(4)} solid`);
    expect(fit.roles[0].name, 'Solid',
      `noise put ${fit.roles[0].name} first: ${fit.roles.map((r) => `${r.name} ${r.spread.toFixed(4)}`)}`);
    expect(thin.spread < 0.01, true, `a role with no effect reports spread ${thin.spread}`);
  });

  /** Everything a `SummaryPanel` reads, over matrices built by fragmentation: no fragment column ran
   *  across the top, so `axisRole` and `coreRole` are null and navigation is inert. */
  const summaryHost = (matrices: SarMatrix[], rowCount: number): SummaryHost => ({
    matrices,
    matrixRoots: matrices.map((m) => m.parentId ?? m.id),
    matrixTiers: matrices.map((m) => m.level),
    transferSummary: {scanned: false, count: 0},
    assayedCount: matrices.length ? 1 : 0,
    selectRoleValue: () => {},
    higherIsBetter: true,
    scalingLabel: 'none',
    activityIsLog: true,
    activityColumnName: 'activity',
    hostRowCount: rowCount,
    unscalableCount: 0,
    axisRole: null,
    coreRole: null,
    // Read the same way the viewer reads it, so the stub cannot disagree with the host it stands in for.
    coresAreSeries: matrices.length > 0 &&
      matrices.every((m) => new Set(m.rows.map((r) => r.coreSmiles)).size === 1),
    roleColumns: [],
    setColumnAxis: () => {},
    predictVirtual: true,
    predictUnmeasured: true,
    computing: false,
    noMatricesMessage: () => 'no matrices',
    formatActivity: (value: number) => value.toFixed(2),
    cellIdText: () => null,
    observedNeighbours: () => 0,
    revealCell: () => {},
    revealMatrix: () => {},
    showTab: () => {},
    addCellsToMakeList: () => {},
  });

  test('the role fit is absent in fragmentation mode', async () => {
    const {molecules, activity} = mixedCoverage();
    const matrices = await runSarMatrix(molecules, activity, e2eParams(false));
    expect(matrices.length > 0, true, 'the fragmentation pipeline must produce matrices to summarise');
    const panel = new SummaryPanel(summaryHost(matrices, molecules.length));
    try {
      panel.activateSummaryTab();
      // Collection runs off the activation stack so the loader can paint.
      await delay(50);
      const segment = Array.from(panel.root.querySelectorAll<HTMLElement>('.chem-sar-sum-seg'))
        .find((el) => el.textContent === 'Effects');
      expect(segment !== undefined, true, 'the segment bar must offer the Effects segment');
      segment!.click();
      const pane = panel.root.querySelector('.chem-sar-sum-pane')!;
      const titles = Array.from(pane.querySelectorAll('.chem-sar-sum-card-title'))
        .map((el) => el.textContent ?? '');
      // A substituent label discovered by fragmentation is local to its own series, so there is no
      // scale on which one global offset could be read and the cards are absent rather than empty.
      expect(titles.some((t) => t.includes('offsets from the additive fit')), false,
        `a role card rendered without a fragment column: ${titles.join(' | ')}`);
      expect(titles.length, 3, `Effects holds its three original cards, got ${titles.join(' | ')}`);
    } finally {
      panel.release();
    }
  });

  test('cross-validation reacts to non-additivity', async () => {
    const clean = fitRoleEffects(roleFixture())!;
    expect(clean.cvR2! > 0.9, true, `an additive table cross-validates, got ${clean.cvR2}`);
    const spiked = (size: number): RoleFit => {
      const design = roleFixture();
      design.activity[0] += size;
      return fitRoleEffects(design)!;
    };
    const mild = spiked(0.5);
    // One cell no combination of role offsets can reach: the fit absorbs part of it in sample and pays
    // for the rest on the held-out compound, so the in-sample error stays under the cross-validated one
    // and the cross-validated R² is the number that moves.
    expect(mild.residualSd < mild.cvRmse!, true, 'an in-sample error cannot be the gate');
    expect(mild.cvR2! < clean.cvR2! - 0.1, true, `cvR2 barely moved: ${mild.cvR2}`);
    expect(mild.roles.map((r) => r.name).join(','), 'Core,Warhead,Ligand',
      'the leaderboard still looks orderly, which is why the gate cannot be the ranking itself');
    // The Summary tab withholds the ranking below 0.5 rather than ranking a table the additive reading
    // does not describe.
    expect(spiked(1.5).cvR2! < 0.5, true);
  });
});
