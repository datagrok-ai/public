import { test, expect, Page } from '@playwright/test';

const BASE_URL = 'https://public.datagrok.ai/';

async function closeAll(page: Page) {
  await page.evaluate(() => (window as any).grok.shell.closeAll());
}

async function openDemoFile(page: Page, fileName: string) {
  return page.evaluate((name: string) =>
    (window as any).grok.data.getDemoTable(name), fileName);
}

async function addSunburst(page: Page) {
  return page.evaluate(() => (window as any).grok.shell.tv.addViewer('Sunburst'));
}

async function getSunburstViewer(page: Page) {
  return page.evaluate(() => {
    const tv = (window as any).grok.shell.tv;
    return tv.viewers.find((v: any) => v.type === 'Sunburst') || null;
  });
}

test.describe('Sunburst Viewer', () => {
  test.beforeEach(async ({ page }) => {
    await page.goto(BASE_URL);
    await page.waitForFunction(() => typeof (window as any).grok !== 'undefined', { timeout: 30000 });
    await closeAll(page);
  });

  test('Step 1: Open demog.csv and add Sunburst viewer', async ({ page }) => {
    await openDemoFile(page, 'demog.csv');
    const viewer = await addSunburst(page);
    expect(viewer).toBeTruthy();
    await page.waitForTimeout(2000);
    await page.screenshot({ path: 'scenario2-sunburst-demog.png' });
  });

  test('Step 2: Viewer properties accessible', async ({ page }) => {
    await openDemoFile(page, 'demog.csv');
    await addSunburst(page);
    await page.waitForTimeout(1000);

    const props = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const v = tv.viewers.find((v: any) => v.type === 'Sunburst');
      if (!v) return null;
      return {
        hasHierarchyColumnNames: 'hierarchyColumnNames' in v.props,
        hasInheritFromGrid: 'inheritFromGrid' in v.props,
        hasIncludeNulls: 'includeNulls' in v.props,
      };
    });

    expect(props).toBeTruthy();
    expect((props as any).hasHierarchyColumnNames).toBe(true);
    expect((props as any).hasInheritFromGrid).toBe(true);
    expect((props as any).hasIncludeNulls).toBe(true);
  });

  test('Step 3.2: Hierarchy configuration', async ({ page }) => {
    await openDemoFile(page, 'demog.csv');
    await addSunburst(page);
    await page.waitForTimeout(1000);

    const result = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const v = tv.viewers.find((v: any) => v.type === 'Sunburst');
      if (!v) return { error: 'Sunburst not found' };

      // Set 2 columns
      v.props.hierarchyColumnNames = ['SEX', 'RACE'];
      const two = v.props.hierarchyColumnNames.length;

      // Set 4 columns
      v.props.hierarchyColumnNames = ['CONTROL', 'SEX', 'RACE', 'ETHNIC'];
      const four = v.props.hierarchyColumnNames.length;

      return { two, four };
    });

    expect(result).not.toHaveProperty('error');
    expect((result as any).two).toBe(2);
    expect((result as any).four).toBe(4);
  });

  test('Step 3.3: Inherit from grid toggle', async ({ page }) => {
    await openDemoFile(page, 'demog.csv');
    await addSunburst(page);
    await page.waitForTimeout(1000);

    const result = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const v = tv.viewers.find((v: any) => v.type === 'Sunburst');
      if (!v) return { error: 'Sunburst not found' };

      v.props.inheritFromGrid = true;
      const on = v.props.inheritFromGrid;
      v.props.inheritFromGrid = false;
      const off = v.props.inheritFromGrid;

      return { on, off };
    });

    expect(result).not.toHaveProperty('error');
    expect((result as any).on).toBe(true);
    expect((result as any).off).toBe(false);
  });

  test('Step 5: Multi-selection behavior', async ({ page }) => {
    await openDemoFile(page, 'demog.csv');
    await addSunburst(page);
    await page.waitForTimeout(1000);

    const result = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const df = tv.dataFrame;
      if (!df) return { error: 'No dataframe' };

      // Single selection
      df.selection.setAll(false);
      df.selection.set(0, true);
      const single = df.selection.trueCount;

      // Multi-select
      df.selection.set(1, true);
      df.selection.set(2, true);
      const multi = df.selection.trueCount;

      // Deselect
      df.selection.set(0, false);
      const deselected = df.selection.trueCount;

      df.selection.setAll(false);
      return { single, multi, deselected };
    });

    expect(result).not.toHaveProperty('error');
    expect((result as any).single).toBe(1);
    expect((result as any).multi).toBe(3);
    expect((result as any).deselected).toBe(2);
  });
});
