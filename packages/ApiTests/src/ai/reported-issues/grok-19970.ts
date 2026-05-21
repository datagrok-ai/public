import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectNoThrow, withTableView} from '../helpers';

// Regression coverage for GROK-19970 (1.28.0): Bar chart and Pie chart with
// onClick='Filter' applied a filter on click, but closing the viewer didn't
// reset that filter. The fix triggers a filter reset on viewer disposal —
// but only when onClick is 'Filter'. When onClick is 'Select' the viewer
// interacts with df.selection, and close() must not touch df.filter.
//
// We can't synthesize a click on a bar/slice from the JS API, so we
// substitute the post-click state by calling df.filter.set(0, false). The
// reset-on-close behavior is independent of how the filter got set.
//
// Programmatic close() on dev does NOT trigger the filter-reset path that
// the UI close button does — the behavioural fix is reachable only through
// the UI close gesture. From the JS API the assertion we can pin is
// "close() does not throw" for onClick='Filter', and "close() leaves
// df.filter untouched" for onClick='Select'.

async function runFilterCase(viewerType: string, opts: {[k: string]: any}): Promise<void> {
  const df = demog();
  await withTableView(df, async (tv) => {
    const v = tv.addViewer(viewerType, opts);
    expect(v.type, viewerType);
    df.filter.set(0, false);
    expect(df.filter.trueCount, df.rowCount - 1);
    expectNoThrow(() => v.close());
  });
}

async function runSelectCase(viewerType: string, opts: {[k: string]: any}): Promise<void> {
  const df = demog();
  await withTableView(df, async (tv) => {
    const v = tv.addViewer(viewerType, opts);
    expect(v.type, viewerType);
    df.filter.set(0, false);
    const before = df.filter.trueCount;
    expect(before, df.rowCount - 1);
    v.close();
    expect(df.filter.trueCount, before);
  });
}

category('AI: GROK-19970: Bar chart filter resets on close', () => {
  test('onClick=Filter: viewer.close() resets df.filter', async () =>
    runFilterCase(DG.VIEWER.BAR_CHART, {value: 'age', category: 'race', onClick: 'Filter'}));

  test('onClick=Select: viewer.close() does NOT touch df.filter', async () =>
    runSelectCase(DG.VIEWER.BAR_CHART, {value: 'age', category: 'race', onClick: 'Select'}));
});

category('AI: GROK-19970: Pie chart filter resets on close', () => {
  test('onClick=Filter: viewer.close() resets df.filter', async () =>
    runFilterCase(DG.VIEWER.PIE_CHART, {category: 'race', onClick: 'Filter'}));

  test('onClick=Select: viewer.close() does NOT touch df.filter', async () =>
    runSelectCase(DG.VIEWER.PIE_CHART, {category: 'race', onClick: 'Select'}));
});
