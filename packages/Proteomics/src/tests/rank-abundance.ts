import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {findAbundanceByCondition} from '../utils/abundance-detection';
import {dockRankAbundanceCharts} from '../viewers/rank-abundance';
import {setGroups} from '../analysis/experiment-setup';
import {SEMTYPE} from '../utils/proteomics-types';

/** Candidates-shaped mock: protein id + the two AVG Group Quantity columns +
 * named groups (numerator / denominator). */
function makeCandidatesDf(num: number[], den: number[], g1 = 'DMD', g2 = 'WT'): DG.DataFrame {
  const prot = num.map((_, i) => `P${i + 1}`);
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('ProteinGroups', prot),
    DG.Column.fromFloat32Array('AVG Group Quantity Numerator', new Float32Array(num)),
    DG.Column.fromFloat32Array('AVG Group Quantity Denominator', new Float32Array(den)),
  ]);
  df.col('ProteinGroups')!.semType = SEMTYPE.PROTEIN_ID;
  setGroups(df, {group1: {name: g1, columns: []}, group2: {name: g2, columns: []}});
  return df;
}

category('Rank-Abundance', () => {
  test('resolver reads per-condition abundance from Candidates group quantities', async () => {
    const df = makeCandidatesDf([100, 1000, 10], [50, 500, 5]);
    const ab = findAbundanceByCondition(df);
    expect(ab !== null, true);
    expect(ab!.group1.name, 'DMD');
    expect(ab!.group2.name, 'WT');
    // log10(100) = 2, log10(1000) = 3
    expect(Math.abs((ab!.group1.log10Intensity[0] as number) - 2) < 1e-4, true);
    expect(Math.abs((ab!.group1.log10Intensity[1] as number) - 3) < 1e-4, true);
    expect(Math.abs((ab!.group2.log10Intensity[0] as number) - Math.log10(50)) < 1e-4, true);
  });

  test('resolver returns null with no abundance columns and no group columns', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('ProteinGroups', ['P1', 'P2']),
      DG.Column.fromFloat32Array('log2FC', new Float32Array([1.0, -1.0])),
    ]);
    expect(findAbundanceByCondition(df), null);
  });

  test('resolver maps non-positive abundance to null', async () => {
    const df = makeCandidatesDf([100, 0, -5], [10, 10, 10]);
    const ab = findAbundanceByCondition(df)!;
    expect(ab.group1.log10Intensity[1], null); // 0 → null
    expect(ab.group1.log10Intensity[2], null); // negative → null
  });

  test('dockRankAbundanceCharts adds rank column (1 = highest abundance) and docks', async () => {
    const df = makeCandidatesDf([100, 1000, 10], [50, 500, 5]);
    df.name = 'Rank-abundance test';
    const tv = grok.shell.addTableView(df);
    try {
      const docked = dockRankAbundanceCharts(tv, df);
      expect(docked, true);
      const rankCol = df.col('rank: DMD');
      expect(rankCol !== null, true);
      // abundances [100,1000,10] → highest (1000) is rank 1, 100 → rank 2, 10 → rank 3
      expect(rankCol!.get(1), 1);
      expect(rankCol!.get(0), 2);
      expect(rankCol!.get(2), 3);
      // two rank-abundance scatter plots docked
      let scatters = 0;
      for (const v of tv.viewers) if (v.type === DG.VIEWER.SCATTER_PLOT) scatters++;
      expect(scatters, 2);
    } finally {
      tv.close();
    }
  });

  test('dockRankAbundanceCharts returns false when no abundance', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('ProteinGroups', ['P1', 'P2']),
    ]);
    df.name = 'No-abundance test';
    const tv = grok.shell.addTableView(df);
    try {
      expect(dockRankAbundanceCharts(tv, df), false);
    } finally {
      tv.close();
    }
  });
});
