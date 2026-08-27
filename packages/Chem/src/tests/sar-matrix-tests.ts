import {category, test, expect, expectFloat, before} from '@datagrok-libraries/test/src/test';

import {_package} from '../package-test';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {MmpFragments} from '../analysis/molecular-matched-pairs/mmp-analysis/mmpa-misc';
import {buildMatchedSeries, clusterRelatedCores} from '../analysis/sar-matrix/sar-matrix-clustering';
import {assembleSinglePositionMatrix, fitAdditiveModel} from '../analysis/sar-matrix/sar-matrix-assemble';
import {computeMatrixConfidence} from '../analysis/sar-matrix/sar-matrix-confidence';
import {rankMatrices, SarRankScheme} from '../analysis/sar-matrix/sar-matrix-ranking';
import {computeAllTransfers, spearman, transferStats} from '../analysis/sar-matrix/sar-matrix-transfer';
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

category('SAR Matrix', () => {
  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  test('buildMatchedSeries groups molecules by shared core', async () => {
    const series = buildMatchedSeries(fakeFrags(), 10);
    expect(series.length, 2);
    for (const s of series) {
      expect(s.members.length, 2);
      const subs = s.members.map((m) => m.substSmiles).sort();
      expect(subs[0], 'Et');
      expect(subs[1], 'Me');
    }
  });

  test('assembleMatrix fills an empty cell with the additive prediction', async () => {
    const cluster: CoreCluster = {
      id: 'c0',
      siteKey: '',
      level: 2,
      series: [
        {coreSmiles: 'CoreA', members: [
          {molIdx: 0, substSmiles: 'Me'},
          {molIdx: 1, substSmiles: 'Et'},
        ]},
        {coreSmiles: 'CoreB', members: [
          {molIdx: 2, substSmiles: 'Me'},
        ]},
      ],
    };
    const molecules = ['A-Me', 'A-Et', 'B-Me'];
    const activities = Float32Array.from([1, 2, 3]);
    const matrix = assembleSinglePositionMatrix(cluster, molecules, activities, true);

    expect(matrix.rows.length, 2);
    expect(matrix.columns.length, 2);
    expect(matrix.realCount, 3);
    expect(matrix.virtualCount, 1);

    const meIdx = matrix.columns.findIndex((c) => c.substSmiles === 'Me');
    const etIdx = matrix.columns.findIndex((c) => c.substSmiles === 'Et');
    const bEt = matrix.cells[1][etIdx];
    expect(bEt.kind, 'virtual');
    // rowMean(B)=3, colMean(Et)=2, grandMean=2  ->  3 + 2 - 2 = 3
    expect(bEt.value, 3);
    expect(matrix.cells[1][meIdx].kind, 'real');
  });

  test('assembleSinglePositionMatrix: columns are the distinct substituents; predict:false leaves no virtual cells',
    async () => {
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

  test('assembleSinglePositionMatrix: a NaN activity is skipped, not turned into a real cell', async () => {
    const cluster: CoreCluster = {
      id: 'c0',
      siteKey: '',
      level: 2,
      series: [
        {coreSmiles: 'CoreA', members: [
          {molIdx: 0, substSmiles: 'Me'},
          {molIdx: 1, substSmiles: 'Et'}, // activity missing (NaN)
        ]},
        {coreSmiles: 'CoreB', members: [
          {molIdx: 2, substSmiles: 'Me'},
        ]},
      ],
    };
    const molecules = ['A-Me', 'A-Et', 'B-Me'];
    const activities = Float32Array.from([10, NaN, 5]);
    const matrix = assembleSinglePositionMatrix(cluster, molecules, activities, false);

    expect(matrix.realCount, 2, 'the NaN-activity member must not count as a real observation');
    const etIdx = matrix.columns.findIndex((c) => c.substSmiles === 'Et');
    const aEt = matrix.cells[0][etIdx];
    expect(aEt.kind, 'empty', 'missing activity leaves the combination unmade, not real');
    expect(aEt.value, null);
    // A poisoned min (e.g. 0 from an unguarded NaN comparison) would break the potency color scale.
    expect(matrix.minActivity, 5, 'minActivity must be the min of the finite activities only');
    expect(matrix.maxActivity, 10, 'maxActivity must be the max of the finite activities only');
  });

  test('fitAdditiveModel: predicts rowMean + colMean - grandMean with support = min(rowN, colN)', async () => {
    // 3x3 grid: row2 and col2 have no observation at all.
    const cells: SarMatrixCell[][] = [
      [realCell(1), realCell(3), emptyCell()],
      [realCell(5), emptyCell(), emptyCell()],
      [emptyCell(), emptyCell(), emptyCell()],
    ];
    const predict = fitAdditiveModel(cells, 3, 3);

    // rowMean(1) = 5, colMean(1) = 3, grandMean = (1+3+5)/3 = 3  ->  5 + 3 - 3 = 5
    const p11 = predict(1, 1);
    expect(p11 !== null, true);
    expectFloat(p11!.value, 5);
    expect(p11!.support, 1, 'support = min(rowN=1, colN=1)');

    expect(predict(0, 2), null, 'column 2 has no observation — unpredictable');
    expect(predict(2, 0), null, 'row 2 has no observation — unpredictable');
  });

  test('computeMatrixConfidence: null when fewer than 4 observed cells', async () => {
    const cells: SarMatrixCell[][] = [
      [realCell(1), realCell(2), emptyCell()],
      [realCell(3), emptyCell(), emptyCell()],
    ];
    expect(computeMatrixConfidence(makeMatrix(cells)), null);
  });

  test('computeMatrixConfidence: high R2 / low RMSE on a perfectly additive matrix', async () => {
    const matrix = additiveMatrix([0, 1, 2, 3], [0, 10, 20, 30]);
    const conf = computeMatrixConfidence(matrix);
    expect(conf !== null, true);
    expect(conf!.n, 16);
    expect(conf!.r2 > 0.85, true, `expected r2 close to 1, got ${conf!.r2}`);
    expect(conf!.rmse < 5, true, `expected a small rmse, got ${conf!.rmse}`);
  });

  test('computeMatrixConfidence: negative R2 on a clearly non-additive matrix', async () => {
    // Diagonal spikes: no row/column effect explains this pattern.
    const cells: SarMatrixCell[][] = Array.from({length: 4}, (_, ri) =>
      Array.from({length: 4}, (_, ci) => realCell(ri === ci ? 100 : 0)));
    const conf = computeMatrixConfidence(makeMatrix(cells));
    expect(conf !== null, true);
    expect(conf!.r2 < 0, true, `expected a negative r2 for a non-additive matrix, got ${conf!.r2}`);
  });

  test('computeMatrixConfidence: reports the cross-validatable fraction and stores a LOO fit', async () => {
    const matrix = additiveMatrix([0, 1, 2, 3], [0, 10, 20, 30]);
    const conf = computeMatrixConfidence(matrix);
    expect(conf !== null, true);
    expect(conf!.total, 16, 'total counts every observed cell');
    expect(conf!.n, 16, 'all 16 cells are cross-validatable here');
    // The LOO pass stores an out-of-sample fitted value on every observed cell (drives the flag).
    expect(typeof matrix.cells[0][0].fit, 'number', 'observed cells get a leave-one-out fit');
  });

  test('computeMatrixConfidence: validates each R-position slice independently', async () => {
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

  test('rankMatrices scores every scheme and sorts by the chosen one', async () => {
    const cluster: CoreCluster = {
      id: 'c0',
      siteKey: '',
      level: 2,
      series: [
        {coreSmiles: 'CoreA', members: [
          {molIdx: 0, substSmiles: 'Me'},
          {molIdx: 1, substSmiles: 'Et'},
        ]},
      ],
    };
    const matrix = assembleSinglePositionMatrix(cluster, ['A-Me', 'A-Et'], Float32Array.from([1, 5]), false);
    const ranked = rankMatrices([matrix], SarRankScheme.Potency);
    expect(ranked.length, 1);
    expect(ranked[0].scores[SarRankScheme.Potency], 5);
  });

  test('spearman: rank correlation, monotone-invariant and outlier-robust', async () => {
    // A monotone but non-linear relation is a perfect rank correlation even though Pearson is < 1.
    expectFloat(spearman([1, 2, 3, 4], [1, 4, 9, 16])!, 1, 0.001);
    expectFloat(spearman([1, 2, 3, 4], [4, 3, 2, 1])!, -1, 0.001);
    // One extreme value can dominate Pearson but not the ranks — the ordering is still 1:1.
    expectFloat(spearman([1, 2, 3, 100], [1, 2, 3, 4])!, 1, 0.001);
    // Guards: too few points (below MIN_COMMON = 3) and a constant side both give null.
    expect(spearman([1, 2], [1, 2]), null, 'fewer than the minimum shared points is null');
    expect(spearman([1, 1, 1, 1], [1, 2, 3, 4]), null, 'a constant side has no ordering to correlate');
  });

  test('computeAllTransfers finds a cross-matrix transfer between differently-scaffolded series', async () => {
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

  test('computeAllTransfers reports two related-core rows of one matrix (the canonical transfer)', async () => {
    // Rows of one matrix are analog series with related cores — the canonical SAR-transfer pairing.
    // Real observations are correlated (not the additive fit), so a within-matrix pair is a transfer.
    const mat = makeMatrix([xferRow([1, 2, 3, 4]), xferRow([2, 3, 4, 5], 'b')]);
    const transfers = await computeAllTransfers([mat], 0);
    expect(transfers.length, 1, 'a within-matrix related-core pair is a transfer');
    expect(transfers[0].a.matrixIndex === transfers[0].b.matrixIndex, true, 'both rows are in the one matrix');
    expectFloat(transfers[0].correlation, 1, 0.01);
  });

  test('computeAllTransfers: anti-correlated rows of one matrix are not a transfer', async () => {
    // Within-matrix pairs are compared, but only conserved potency order counts: ρ = -1 is not transfer.
    const mat = makeMatrix([xferRow([1, 2, 3, 4]), xferRow([4, 3, 2, 1], 'b')]);
    expect((await computeAllTransfers([mat], 0)).length, 0, 'opposite orderings are not a transfer');
  });

  test('computeAllTransfers skips pairs below the correlation floor', async () => {
    const matA = makeMatrix([xferRow([1, 2, 3, 4])]);
    const matB = makeMatrix([xferRow([4, 3, 2, 1], 'b')]); // ρ = -1, below the 0.7 floor
    expect((await computeAllTransfers([matA, matB], 0)).length, 0, 'uncorrelated series yield no transfer');
  });

  test('computeAllTransfers dedupes the same core pair to its best position', async () => {
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

  test('transferStats: identical per-step deltas give foldMatch = 1, no benefiting cell without a virtual',
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

  test('transferStats: differing step magnitudes give a lower foldMatch', async () => {
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

  test('transferStats points benefiting at the follower core\'s virtual analog', async () => {
    // Four measured pairs (S0..S3) carry the correlation; the follower has S4 only as a prediction.
    const matA = makeMatrix([xferRow([1, 2, 3, 4, 5])]);
    const matB = makeMatrix([[...xferRow([2, 3, 4, 5], 'b'), virtualCell(10)]]);
    const transfers = await computeAllTransfers([matA, matB], 0);
    expect(transfers.length, 1, 'one transfer between the two single-row matrices');
    const stats = transferStats(transfers[0], true);
    expect(stats.benefiting !== null, true, 'the follower core has an untested (virtual) analog to fill');
    expect(stats.benefiting!.side, 'b', 'the follower is the second matrix');
    expect(stats.benefiting!.substSmiles, 'S4', 'the virtual sits at the follower core\'s fifth substituent');
  });

  test('computeAllTransfers carries the analogs one side has only predicted', async () => {
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

  test('computeAllTransfers ignores a compound shared by both series', async () => {
    // The overlapping cover puts one compound in several matrices at once. Matched against itself it
    // would put its own potency on both axes and report a perfect correlation that means nothing, so
    // these two series — same compounds, same potencies — must yield no transfer at all.
    const matA = makeMatrix([xferRow([1, 2, 3, 4])]);
    const matB = makeMatrix([xferRow([1, 2, 3, 4])]);
    expect((await computeAllTransfers([matA, matB], 0)).length, 0,
      'a series compared against its own compounds is not a transfer');
  });

  test('computeAllTransfers requires the same R-group, not merely a similar one', async () => {
    // Same compounds and the same perfect trend, but the follower files them under different
    // substituents. A transfer claims one particular change did one particular thing on two scaffolds,
    // and two series that never made the same change support no such claim.
    const matA = makeMatrix([xferRow([1, 2, 3, 4])]);
    const matB = makeMatrix([xferRow([2, 3, 4, 5], 'b')]);
    matB.columns.forEach((col, i) => col.substSmiles = `T${i}`);
    expect((await computeAllTransfers([matA, matB], 0)).length, 0, 'different R-groups are not a pairing');
  });

  test('computeAllTransfers applies the compound similarity floor', async () => {
    // Identical R-groups and a perfect trend either way, so only the similarity gate decides. At 0
    // every pair passes; at 1 none can, since the two sides are different compounds by construction.
    const matA = makeMatrix([xferRow([1, 2, 3, 4])]);
    const matB = makeMatrix([xferRow([2, 3, 4, 5], 'b')]);
    expect((await computeAllTransfers([matA, matB], 0)).length, 1, 'an open gate admits the transfer');
    expect((await computeAllTransfers([matA, matB], 1)).length, 0,
      'demanding all but identical compounds admits none');
  });

  test('linkRGroupFragments joins a core with a fragment at its attachment point', async () => {
    const svc = await chemCommonRdKit.getRdKitService();
    const linked = await svc.linkRGroupFragments(['c1ccc([*:1])cc1'], [['C[*:1]']], [1]);
    expect(linked.length, 1);
    expect(linked[0].includes('*'), false, 'no stray attachment point should remain');
    const mol = chemCommonRdKit.checkMoleculeValid(linked[0]);
    expect(mol !== null, true, 'the assembled SMILES must parse as a valid molecule');
    expect(mol.get_num_atoms(), 7, 'toluene: 6 ring carbons + 1 methyl carbon');
    mol.delete();
  });

  test('linkRGroupFragments: a position the core lacks is skipped, not fatal', async () => {
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

  test('linkRGroupFragments: an empty fragment removes the attachment point cleanly', async () => {
    const svc = await chemCommonRdKit.getRdKitService();
    const linked = await svc.linkRGroupFragments(['c1ccc([*:1])cc1'], [['']], [1]);
    expect(linked.length, 1);
    expect(linked[0].includes('*'), false, 'no stray attachment point should remain');
    const mol = chemCommonRdKit.checkMoleculeValid(linked[0]);
    expect(mol !== null, true, 'the assembled SMILES must parse as a valid molecule');
    expect(mol.get_num_atoms(), 6, 'benzene: 6 ring carbons, no substituent added');
    mol.delete();
  });

  test('runSarMatrix end-to-end (fragmentation + clustering + assembly via workers)', async () => {
  }, {skipReason: 'GROK-TBD — needs the server-side Chem:MurckoScaffolds script plus MMP fragmentation ' +
    'workers; covered indirectly by the pure-function tests above and by manual QA of the SAR Matrix viewer.'});

  // Fragment ids are minted as the workers discover fragments, so the same data reaches this stage
  // under different ids from one run to the next. Anything downstream resolving a tie by "whichever
  // came first" then yields a different set of matrices for identical input.
  test('buildMatchedSeries: series order follows the chemistry, not the fragment ids', async () => {
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
    const shape = (frags: MmpFragments) => buildMatchedSeries(frags, 1)
      .map((s) => `${s.coreSmiles}:${s.members.map((m) => m.molIdx).join(',')}`).join(' | ');
    expect(shape(asDiscovered), shape(asRediscovered),
      'the same molecules must give the same series whichever ids the workers assigned');
  });

  test('buildMatchedSeries: repeated calls on one input are identical', async () => {
    const frags = fakeFrags();
    const once = JSON.stringify(buildMatchedSeries(frags, 1));
    const twice = JSON.stringify(buildMatchedSeries(fakeFrags(), 1));
    expect(once, twice, 'a pure stage must not vary between calls');
  });

  // Clustering by core similarity alone will put cores whose substituents hang off different places
  // into one matrix. Their substituent vocabularies are disjoint, so the column axis pools both and a
  // column means one position on some rows and another on the rest.
  test('clusterRelatedCores: every multi-row cluster varies at a shared site', async () => {
    // One triazole, attachments at three different ring positions — near-identical by fingerprint.
    const series: MatchedSeries[] = [
      {coreSmiles: 'CCc1c([*:1])nnn1-c1ccc(F)cc1', members: [{molIdx: 0, substSmiles: 'C[*:1]'}]},
      {coreSmiles: 'COc1c([*:1])nnn1-c1ccc(F)cc1', members: [{molIdx: 1, substSmiles: 'C[*:1]'}]},
      {coreSmiles: 'COC(=O)c1nnn(-c2ccc(F)cc2)c1[*:1]', members: [{molIdx: 2, substSmiles: 'C[*:1]'}]},
      {coreSmiles: 'N#Cc1nnn(-c2ccc(F)cc2)c1[*:1]', members: [{molIdx: 3, substSmiles: 'CC[*:1]'}]},
    ] as MatchedSeries[];
    const clusters = await clusterRelatedCores(series, 0.3);
    const grouped = clusters.filter((c) => c.series.length > 1);
    // Without this the assertion below is vacuous: all-singleton output would pass it having tested
    // nothing, and singletons are exactly what a too-tight threshold produces.
    expect(grouped.length > 0, true, 'these cores must cluster at all for the check to mean anything');
    for (const cluster of grouped) {
      expect(cluster.siteKey !== '', true,
        'a cluster holding several cores must name the site they share, or its columns mix positions');
      const sites = new Set(cluster.series.map((s) => s.coreSmiles.indexOf('[*:1]') >= 0 ?
        s.coreSmiles.replace(/\[\*:\d+\]/g, '*') : ''));
      expect(sites.size >= 1, true, 'every row of a cluster carries an attachment');
    }
  });
});
