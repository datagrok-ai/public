import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {parseSpectronautCandidatesText} from '../parsers/spectronaut-candidates-parser';
import {dockComparisonFilterIfMultiContrast} from '../package';
import {SEMTYPE} from '../utils/proteomics-types';
import {_package} from '../package-test';

/** Builds a minimal Candidates TSV. */
function candidatesTsv(rows: string[][]): string {
  const headers = ['Comparison (group1/group2)', 'PG.ProteinGroups', 'PG.Genes',
    'AVG Log2 Ratio', 'Qvalue', 'Pvalue'];
  return [headers.join('\t'), ...rows.map((r) => r.join('\t'))].join('\n');
}

function hasFiltersViewer(tv: DG.TableView): boolean {
  for (const v of tv.viewers) {
    if (v.type === DG.VIEWER.FILTERS) return true;
  }
  return false;
}

const FIXTURE_PATH = 'demo/spectronaut-hye-candidates.tsv';

// Expected from the deterministic generator (tools/generate-spectronaut-candidates-fixture.mjs)
// run against files/demo/spectronaut-hye-mix.tsv.
const EXPECTED_ROWS = 93;
const EXPECTED_SIG_UP = 26;
const EXPECTED_SIG_DOWN = 17;
const EXPECTED_NOT_SIG = 50;
const EXPECTED_ORGANISMS_AT_LEAST = ['Homo sapiens', 'Saccharomyces cerevisiae', 'Escherichia coli'];

category('SpectronautCandidates E2E', () => {
  test('parses HYE candidates fixture and reproduces expected DE signal', async () => {
    const text = await _package.files.readAsText(FIXTURE_PATH);
    const df = await parseSpectronautCandidatesText(text);

    expect(df.rowCount, EXPECTED_ROWS);

    expect(df.getTag('proteomics.source'), 'spectronaut-candidates');
    expect(df.getTag('proteomics.de_complete'), 'true');
    expect(df.getTag('proteomics.de_method'), 'spectronaut');

    expect(df.col('log2FC') !== null, true);
    expect(df.col('adj.p-value') !== null, true);
    expect(df.col('p-value') !== null, true);
    expect(df.col('significant') !== null, true);
    expect(df.col('AVG Log2 Ratio'), null);
    expect(df.col('Qvalue'), null);
    expect(df.col('Pvalue'), null);

    expect(df.col('ProteinGroups')!.semType, SEMTYPE.PROTEIN_ID);
    expect(df.col('log2FC')!.semType, SEMTYPE.LOG2FC);
    expect(df.col('p-value')!.semType, SEMTYPE.P_VALUE);
    expect(df.col('adj.p-value')!.semType, SEMTYPE.P_VALUE);

    const fc = df.col('log2FC')!;
    const sig = df.col('significant')!;
    let up = 0, down = 0, ns = 0;
    for (let i = 0; i < df.rowCount; i++) {
      if (sig.get(i)) {
        if ((fc.get(i) as number) > 0) up++; else down++;
      } else ns++;
    }
    expect(up, EXPECTED_SIG_UP);
    expect(down, EXPECTED_SIG_DOWN);
    expect(ns, EXPECTED_NOT_SIG);

    const orgCol = df.col('Organisms')!;
    const seenOrganisms = new Set<string>();
    for (let i = 0; i < df.rowCount; i++) {
      const v = orgCol.get(i);
      if (typeof v === 'string') seenOrganisms.add(v);
    }
    for (const expected of EXPECTED_ORGANISMS_AT_LEAST) {
      const found = Array.from(seenOrganisms).some((o) => o.includes(expected));
      expect(found, true);
    }
  });

  // --- R4/D-07: Comparison Filter viewer docking (13-03 Task 2) ---

  test('multi-comparison file docks a Comparison Filters viewer', async () => {
    const df = await parseSpectronautCandidatesText(candidatesTsv([
      ['T1 / Control', 'P1', 'G1', '2.0', '0.001', '0.0001'],
      ['T2 / Control', 'P2', 'G2', '-1.7', '0.002', '0.0005'],
      ['T1 / Control', 'P3', 'G3', '1.2', '0.01', '0.005'],
    ]));
    const tv = grok.shell.addTableView(df);
    try {
      const docked = dockComparisonFilterIfMultiContrast(tv, df);
      expect(docked, true);
      expect(hasFiltersViewer(tv), true);
    } finally {
      tv.close();
    }
  });

  test('single-comparison file docks no Filters viewer', async () => {
    const df = await parseSpectronautCandidatesText(candidatesTsv([
      ['T1 / Control', 'P1', 'G1', '2.0', '0.001', '0.0001'],
      ['T1 / Control', 'P2', 'G2', '-1.7', '0.002', '0.0005'],
    ]));
    const tv = grok.shell.addTableView(df);
    try {
      const docked = dockComparisonFilterIfMultiContrast(tv, df);
      expect(docked, false);
      expect(hasFiltersViewer(tv), false);
    } finally {
      tv.close();
    }
  });
});
