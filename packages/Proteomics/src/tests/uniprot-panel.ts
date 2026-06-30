import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {
  renderUniProtWidget,
  renderPerGroupBars,
  findHostDataFrameForProtein,
} from '../panels/uniprot-panel';
import {setGroups} from '../analysis/experiment-setup';
import {SEMTYPE} from '../utils/proteomics-types';

/** Builds a fixture DataFrame mimicking a 4-sample, 3-protein post-DE table.
 *  Cells indexed by row 0 (first protein P11111) have known values so the
 *  test can assert mean/SD computation exactly. */
function makeFixtureDf(opts: {
  groupColsExist?: boolean;
  withGroupsTag?: boolean;
  withProteinId?: boolean;
  row0AllNan?: boolean;
} = {}): DG.DataFrame {
  const {
    groupColsExist = true,
    withGroupsTag = true,
    withProteinId = true,
    row0AllNan = false,
  } = opts;

  const cols: DG.Column[] = [];
  if (withProteinId) {
    const pid = DG.Column.fromStrings('Primary Protein ID', ['P11111', 'P22222', 'P33333']);
    pid.semType = SEMTYPE.PROTEIN_ID;
    cols.push(pid);
  }
  if (groupColsExist) {
    // Row 0 — g1: [1, 3] → mean=2, sd=√2; g2: [5, 7] → mean=6, sd=√2
    const s1 = row0AllNan ? [NaN, 1.0, 1.0] : [1.0, 1.5, 2.0];
    const s2 = row0AllNan ? [NaN, 3.0, 1.0] : [3.0, 2.5, 2.5];
    const s3 = row0AllNan ? [NaN, 4.0, 4.0] : [5.0, 4.0, 5.0];
    const s4 = row0AllNan ? [NaN, 6.0, 4.0] : [7.0, 6.0, 5.5];
    const c1 = DG.Column.fromList('double' as DG.ColumnType, 'Sample1', s1);
    const c2 = DG.Column.fromList('double' as DG.ColumnType, 'Sample2', s2);
    const c3 = DG.Column.fromList('double' as DG.ColumnType, 'Sample3', s3);
    const c4 = DG.Column.fromList('double' as DG.ColumnType, 'Sample4', s4);
    if (row0AllNan) {
      c1.set(0, null); c2.set(0, null); c3.set(0, null); c4.set(0, null);
    }
    cols.push(c1, c2, c3, c4);
  }
  const df = DG.DataFrame.fromColumns(cols);
  if (withGroupsTag && groupColsExist) {
    setGroups(df, {
      group1: {name: 'Control', columns: ['Sample1', 'Sample2']},
      group2: {name: 'Treatment', columns: ['Sample3', 'Sample4']},
    });
  }
  df.name = `panel-fixture-${Math.random().toString(36).slice(2, 8)}`;
  return df;
}

/** Mock UniProtEntry compatible with renderUniProtWidget's expected fields. */
function makeMockEntry(accession: string): any {
  return {
    accession,
    proteinDescription: {recommendedName: {fullName: {value: 'Mock Protein'}}},
    genes: [{geneName: {value: 'MOCK1'}}],
    organism: {scientificName: 'Homo sapiens'},
    comments: [],
    uniProtKBCrossReferences: [],
  };
}

category('Proteomics: 14-04', () => {
  test('uniprotPanelPerGroupBarsRender', async () => {
    const df = makeFixtureDf();
    const tv = grok.shell.addTableView(df);
    try {
      const dom = renderUniProtWidget(makeMockEntry('P11111'), 'P11111');
      // Header
      expect(dom.textContent?.includes('Per-Group Quantities'), true);
      // Two SVG rect bars
      const rects = dom.querySelectorAll('rect');
      expect(rects.length, 2);
      // Two `(n=2)` labels
      const texts = Array.from(dom.querySelectorAll('text'))
        .map((t) => t.textContent ?? '');
      const nTwoCount = texts.filter((t) => t.includes('(n=2)')).length;
      expect(nTwoCount, 2);
      // Two mean/SD text rows in 0.85em — one per non-empty group. Match on
      // exact-leaf divs (no nested div children) so the parent container's
      // concatenated textContent doesn't double-count.
      const meanLines = Array.from(dom.querySelectorAll('div'))
        .filter((d) => d.querySelector(':scope > div') === null)
        .map((d) => d.textContent ?? '')
        .filter((t) => /mean=.*SD=/.test(t));
      expect(meanLines.length, 2);
    } finally {
      tv.close();
    }
  });

  test('uniprotPanelPerGroupBarsEmptyState', async () => {
    // No groups tag — findHostDataFrameForProtein returns null, so the section is omitted entirely.
    const df = makeFixtureDf({withGroupsTag: false});
    const tv = grok.shell.addTableView(df);
    try {
      const dom = renderUniProtWidget(makeMockEntry('P11111'), 'P11111');
      // The section is not even rendered when no host DataFrame is found, so the
      // SVG rect count is zero and the header is absent.
      const rects = dom.querySelectorAll('rect');
      expect(rects.length, 0);
      expect(dom.textContent?.includes('Per-Group Quantities'), false);
    } finally {
      tv.close();
    }
  });

  test('uniprotPanelPerGroupBarsColors', async () => {
    const df = makeFixtureDf();
    const tv = grok.shell.addTableView(df);
    try {
      const dom = renderUniProtWidget(makeMockEntry('P11111'), 'P11111');
      const rects = dom.querySelectorAll('rect');
      expect(rects.length, 2);
      expect(rects[0].getAttribute('fill'), '#FF00FF'); // group1 magenta
      expect(rects[1].getAttribute('fill'), '#00FFFF'); // group2 cyan
    } finally {
      tv.close();
    }
  });

  test('uniprotPanelPerGroupMeanSDComputation', async () => {
    // Custom fixture: row 0 group1=[1,3], group2=[5,7] → mean=2/6, SD=√2≈1.414.
    const cols: DG.Column[] = [];
    const pid = DG.Column.fromStrings('Primary Protein ID', ['P11111']);
    pid.semType = SEMTYPE.PROTEIN_ID;
    cols.push(pid);
    cols.push(DG.Column.fromList('double' as DG.ColumnType, 'S1', [1.0]));
    cols.push(DG.Column.fromList('double' as DG.ColumnType, 'S2', [3.0]));
    cols.push(DG.Column.fromList('double' as DG.ColumnType, 'S3', [5.0]));
    cols.push(DG.Column.fromList('double' as DG.ColumnType, 'S4', [7.0]));
    const df = DG.DataFrame.fromColumns(cols);
    setGroups(df, {
      group1: {name: 'Control', columns: ['S1', 'S2']},
      group2: {name: 'Treatment', columns: ['S3', 'S4']},
    });
    df.name = 'mean-sd-fixture';
    const tv = grok.shell.addTableView(df);
    try {
      const dom = renderPerGroupBars(df, 'P11111');
      expect(dom !== null, true);
      const meanLines = Array.from(dom!.querySelectorAll('div'))
        .map((d) => d.textContent ?? '')
        .filter((t) => /mean=.*SD=/.test(t));
      // Two non-empty groups produce two mean/SD lines.
      expect(meanLines.length, 2);
      // Sample SD (n-1 denominator): variance = ((1-2)^2 + (3-2)^2)/(2-1) = 2 → SD = √2 ≈ 1.41
      const g1Match = meanLines.find((l) => l.startsWith('Control:'));
      const g2Match = meanLines.find((l) => l.startsWith('Treatment:'));
      expect(g1Match !== undefined, true);
      expect(g2Match !== undefined, true);
      expect(g1Match!.includes('mean=2.00'), true);
      expect(g1Match!.includes('SD=1.41'), true);
      expect(g2Match!.includes('mean=6.00'), true);
      expect(g2Match!.includes('SD=1.41'), true);
    } finally {
      tv.close();
    }
  });

  test('uniprotPanelPerGroupBarsAllNanProtein', async () => {
    const df = makeFixtureDf({row0AllNan: true});
    const tv = grok.shell.addTableView(df);
    try {
      const dom = renderPerGroupBars(df, 'P11111');
      // Returns the empty-state element (ui.divText), no SVG.
      expect(dom !== null, true);
      expect(dom!.tagName.toLowerCase(), 'div');
      expect((dom!.textContent ?? '').includes('No per-group quantities available'), true);
      expect(dom!.querySelectorAll('rect').length, 0);
    } finally {
      tv.close();
    }
  });

  test('findHostDataFrameForProteinPrefersActiveTable', async () => {
    const df1 = makeFixtureDf();
    const df2 = makeFixtureDf();
    const tv1 = grok.shell.addTableView(df1);
    const tv2 = grok.shell.addTableView(df2);
    try {
      // tv2 is the most recently added, so it is active.
      const host = findHostDataFrameForProtein('P11111');
      expect(host !== null, true);
      // Both DFs match — preference is the active one (tv.dataFrame).
      const activeDf = grok.shell.tv?.dataFrame;
      expect(host === activeDf, true);
    } finally {
      tv1.close();
      tv2.close();
    }
  });
});
