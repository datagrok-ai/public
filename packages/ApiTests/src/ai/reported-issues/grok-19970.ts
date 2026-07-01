import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectNoThrow, withTableView} from '../helpers';

// Regression coverage for GROK-19970: Bar/Pie chart filter resets on close (onClick='Filter' only).

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
