import {category, test, expect, expectFloat, before} from '@datagrok-libraries/test/src/test';

import {_package} from '../package-test';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {MmpFragments} from '../analysis/molecular-matched-pairs/mmp-analysis/mmpa-misc';
import {buildMatchedSeries} from '../analysis/sar-matrix/sar-matrix-clustering';
import {assembleSinglePositionMatrix, fitAdditiveModel} from '../analysis/sar-matrix/sar-matrix-assemble';
import {computeMatrixConfidence} from '../analysis/sar-matrix/sar-matrix-confidence';
import {rankMatrices, SarRankScheme} from '../analysis/sar-matrix/sar-matrix-ranking';
import {computeAllTransfers, transferStats} from '../analysis/sar-matrix/sar-matrix-transfer';
import {CoreCluster, SarMatrix, SarMatrixCell, SarMatrixColumn, SarMatrixRow}
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

function realCell(value: number, molIdx = 0): SarMatrixCell {
  return {kind: 'real', value, molIdx, smiles: null};
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
  const rows: SarMatrixRow[] = cells.map((_, i) => ({coreFragId: i, coreSmiles: `Core${i}`, label: `Core ${i}`}));
  const columns: SarMatrixColumn[] = cells[0].map((_, i) => ({position: positions[0], substSmiles: `S${i}`, count: 0}));
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
      series: [
        {coreFragId: 1, coreSmiles: 'CoreA', members: [
          {molIdx: 0, substFragId: 3, substSmiles: 'Me'},
          {molIdx: 1, substFragId: 4, substSmiles: 'Et'},
        ]},
        {coreFragId: 2, coreSmiles: 'CoreB', members: [
          {molIdx: 2, substFragId: 3, substSmiles: 'Me'},
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
        series: [
          {coreFragId: 1, coreSmiles: 'CoreA', members: [
            {molIdx: 0, substFragId: 3, substSmiles: 'Me'},
            {molIdx: 1, substFragId: 4, substSmiles: 'Et'},
            {molIdx: 2, substFragId: 5, substSmiles: 'Pr'},
          ]},
          {coreFragId: 2, coreSmiles: 'CoreB', members: [
            {molIdx: 3, substFragId: 3, substSmiles: 'Me'},
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
      series: [
        {coreFragId: 1, coreSmiles: 'CoreA', members: [
          {molIdx: 0, substFragId: 3, substSmiles: 'Me'},
          {molIdx: 1, substFragId: 4, substSmiles: 'Et'}, // activity missing (NaN)
        ]},
        {coreFragId: 2, coreSmiles: 'CoreB', members: [
          {molIdx: 2, substFragId: 3, substSmiles: 'Me'},
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

  test('rankMatrices scores every scheme and sorts by the chosen one', async () => {
    const cluster: CoreCluster = {
      id: 'c0',
      series: [
        {coreFragId: 1, coreSmiles: 'CoreA', members: [
          {molIdx: 0, substFragId: 3, substSmiles: 'Me'},
          {molIdx: 1, substFragId: 4, substSmiles: 'Et'},
        ]},
      ],
    };
    const matrix = assembleSinglePositionMatrix(cluster, ['A-Me', 'A-Et'], Float32Array.from([1, 5]), false);
    const ranked = rankMatrices([matrix], SarRankScheme.Potency);
    expect(ranked.length, 1);
    expect(ranked[0].scores[SarRankScheme.Potency], 5);
  });

  test('computeAllTransfers finds a cross-matrix transfer between differently-scaffolded series', async () => {
    // Two matrices whose columns share substituent SMILES (S0..S2); a row in each tracks the other.
    const matA = makeMatrix([
      [realCell(1), realCell(2), realCell(3)],
      [realCell(5), realCell(4), realCell(3)], // uncorrelated with matB's row, so it isn't the pick
    ]);
    const matB = makeMatrix([[realCell(2), realCell(3), realCell(4)]]);
    const transfers = computeAllTransfers([matA, matB]);
    const cross = transfers.find((t) => t.crossMatrix);
    expect(cross !== undefined, true, 'a cross-series transfer between the two matrices must be found');
    expect(cross!.a.matrixIndex !== cross!.b.matrixIndex, true, 'the two cores come from different matrices');
    expect(cross!.substituents.length, 3, 'all three shared substituents are used');
    expectFloat(cross!.correlation, 1, 0.01);
  });

  test('computeAllTransfers marks a within-matrix pair as not cross-matrix', async () => {
    const mat = makeMatrix([
      [realCell(1), realCell(2), realCell(3)],
      [realCell(2), realCell(3), realCell(4)],
    ]);
    const transfers = computeAllTransfers([mat]);
    expect(transfers.length, 1, 'one within-matrix transfer');
    expect(transfers[0].crossMatrix, false, 'same matrix — not a cross-series transfer');
    expectFloat(transfers[0].correlation, 1, 0.01);
  });

  test('computeAllTransfers skips pairs below the correlation floor', async () => {
    const mat = makeMatrix([
      [realCell(1), realCell(2), realCell(3)],
      [realCell(3), realCell(1), realCell(2)], // ρ ≈ -0.5, below the 0.7 floor
    ]);
    expect(computeAllTransfers([mat]).length, 0, 'uncorrelated rows yield no transfer');
  });

  test('computeAllTransfers dedupes the same core pair to its best position', async () => {
    // Same two cores share substituents at both R1 (ρ=1) and R2 (ρ≈0.98) — only R1 survives.
    const cells: SarMatrixCell[][] = [
      [realCell(1), realCell(2), realCell(3), realCell(1), realCell(2), realCell(4)],
      [realCell(2), realCell(3), realCell(4), realCell(1), realCell(2), realCell(3)],
    ];
    const positions = ['R1', 'R2'];
    const mat = makeMatrix(cells, positions);
    mat.columns.forEach((c, i) => c.position = i < 3 ? 'R1' : 'R2');
    const transfers = computeAllTransfers([mat]);
    expect(transfers.length, 1, 'the two R-positions collapse to one entry for this core pair');
    expect(transfers[0].a.position, 'R1', 'R1 is the more-correlated position, so it is kept');
  });

  test('transferStats: identical per-step deltas give foldMatch = 1, no benefiting cell without a virtual',
    async () => {
      // matB = matA − 9, so the trends are identical: correlation 1 and every step delta matches.
      const matA = makeMatrix([[realCell(10), realCell(12), realCell(14), realCell(20)]]); // deltas +2, +2, +6
      const matB = makeMatrix([[realCell(1), realCell(3), realCell(5), realCell(11)]]); // same deltas
      const transfers = computeAllTransfers([matA, matB]);
      expect(transfers.length, 1, 'one cross-series transfer between the two single-row matrices');
      const stats = transferStats(transfers[0], true);
      expectFloat(stats.correlation, 1, 0.01);
      expect(stats.foldMatch !== null, true);
      expectFloat(stats.foldMatch!, 1, 0.01);
      expect(stats.benefiting, null, 'no virtual cells in either core — nothing to benefit from the transfer');
    });

  test('transferStats: differing step magnitudes give a lower foldMatch', async () => {
    const matA = makeMatrix([[realCell(10), realCell(12), realCell(14), realCell(20)]]); // deltas +2, +2, +6
    const matB = makeMatrix([[realCell(1), realCell(2), realCell(3), realCell(4)]]); // same direction, +1 each
    const transfers = computeAllTransfers([matA, matB]);
    expect(transfers.length, 1);
    const stats = transferStats(transfers[0], true);
    expect(stats.foldMatch !== null, true);
    // steps: min(2,1)/max=0.5, min(2,1)/max=0.5, min(6,1)/max=1/6  ->  mean ~0.389
    expectFloat(stats.foldMatch!, 0.389, 0.01);
    expect(stats.foldMatch! < 1, true);
  });

  test('transferStats points benefiting at the follower core\'s virtual analog', async () => {
    const matA = makeMatrix([[realCell(1), realCell(2), realCell(3), realCell(4)]]);
    const matB = makeMatrix([[realCell(2), realCell(3), realCell(4), virtualCell(10)]]);
    const transfers = computeAllTransfers([matA, matB]);
    expect(transfers.length, 1, 'one cross-series transfer between the two single-row matrices');
    const stats = transferStats(transfers[0], true);
    expect(stats.benefiting !== null, true, 'the follower core has an untested (virtual) analog to fill');
    expect(stats.benefiting!.side, 'b', 'the follower is the second matrix');
    expect(stats.benefiting!.substIndex, 3, 'the virtual sits at the fourth shared substituent');
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
});
