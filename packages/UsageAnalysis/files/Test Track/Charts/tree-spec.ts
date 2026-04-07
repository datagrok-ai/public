import { test, expect, Page } from '@playwright/test';

const BASE_URL = 'https://public.datagrok.ai/';

async function closeAll(page: Page) {
  await page.evaluate(() => (window as any).grok.shell.closeAll());
}

test.describe('Tree Viewer — Collaborative Filtering', () => {
  let page: Page;

  test.beforeEach(async ({ page: p }) => {
    page = p;
    await p.goto(BASE_URL);
    await p.waitForFunction(() => typeof (window as any).grok !== 'undefined', { timeout: 30000 });
    await closeAll(p);
  });

  test('Full collaborative filtering scenario', async ({ page }) => {
    // Setup: open demog.csv
    await page.evaluate(() => (window as any).grok.data.getDemoTable('demog.csv'));
    await page.waitForTimeout(1500);

    // Add Tree viewer
    await page.evaluate(() => (window as any).grok.shell.tv.addViewer('Tree'));
    await page.waitForTimeout(1000);

    // Set hierarchy to CONTROL, SEX, RACE
    const setupResult = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const v = tv.viewers.find((v: any) => v.type === 'Tree');
      if (!v) return { error: 'Tree viewer not found' };
      v.props.hierarchyColumnNames = ['CONTROL', 'SEX', 'RACE'];
      return { hierarchySet: true };
    });
    expect(setupResult).not.toHaveProperty('error');
    await page.waitForTimeout(1000);

    await page.screenshot({ path: 'scenario3-tree-setup.png' });

    // Step 1: Select branches false→F→Asian, false→F→Black, false→M→Asian
    const step1 = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const df = tv.dataFrame;
      if (!df) return { error: 'No dataframe' };

      const controlCol = df.col('CONTROL');
      const sexCol = df.col('SEX');
      const raceCol = df.col('RACE');
      if (!controlCol || !sexCol || !raceCol) return { error: 'Columns not found' };

      df.selection.setAll(false);
      const rowCount = df.rowCount;
      for (let i = 0; i < rowCount; i++) {
        const ctrl = controlCol.get(i);
        const sex = sexCol.get(i);
        const race = raceCol.get(i);
        // false→F→Asian OR false→F→Black OR false→M→Asian
        if (ctrl === false && sex === 'F' && race === 'Asian') df.selection.set(i, true);
        else if (ctrl === false && sex === 'F' && race === 'Black') df.selection.set(i, true);
        else if (ctrl === false && sex === 'M' && race === 'Asian') df.selection.set(i, true);
      }
      return { selectedCount: df.selection.trueCount };
    });
    expect(step1).not.toHaveProperty('error');
    expect((step1 as any).selectedCount).toBeGreaterThan(0);

    // Step 2: Apply CONTROL=true filter → selected ∩ filtered should be 0
    const step2 = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const df = tv.dataFrame;
      const controlCol = df.col('CONTROL');

      // Apply filter: CONTROL = true
      df.filter.setAll(false);
      for (let i = 0; i < df.rowCount; i++) {
        if (controlCol.get(i) === true) df.filter.set(i, true);
      }
      df.filter.fireChanged();

      // Count rows that are both selected and filtered
      let overlap = 0;
      for (let i = 0; i < df.rowCount; i++) {
        if (df.selection.get(i) && df.filter.get(i)) overlap++;
      }
      return { filteredCount: df.filter.trueCount, overlapCount: overlap };
    });
    // Expected: 0 overlap (selected rows all have CONTROL=false)
    expect((step2 as any).overlapCount).toBe(0);

    // Step 3: Add All→true→F→Black selection → overlap should be 2
    const step3 = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const df = tv.dataFrame;
      const controlCol = df.col('CONTROL');
      const sexCol = df.col('SEX');
      const raceCol = df.col('RACE');

      // Add true→F→Black to selection
      for (let i = 0; i < df.rowCount; i++) {
        if (controlCol.get(i) === true && sexCol.get(i) === 'F' && raceCol.get(i) === 'Black')
          df.selection.set(i, true);
      }

      let overlap = 0;
      for (let i = 0; i < df.rowCount; i++) {
        if (df.selection.get(i) && df.filter.get(i)) overlap++;
      }
      return { totalSelected: df.selection.trueCount, overlapCount: overlap };
    });
    // Expected: 2 overlap
    expect((step3 as any).overlapCount).toBe(2);

    // Step 4: Clear CONTROL=true filter → selected count = 176
    const step4 = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const df = tv.dataFrame;

      // Clear filter
      df.filter.setAll(true);
      df.filter.fireChanged();

      return { totalSelected: df.selection.trueCount, totalFiltered: df.filter.trueCount };
    });
    // Expected: 176 selected
    expect((step4 as any).totalSelected).toBe(176);
    expect((step4 as any).totalFiltered).toBe(5850);

    await page.screenshot({ path: 'scenario3-tree-final.png' });
  });
});
