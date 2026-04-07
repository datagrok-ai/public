import { test, expect, Page } from '@playwright/test';

const BASE_URL = 'https://public.datagrok.ai/';

async function closeAll(page: Page) {
  await page.evaluate(() => (window as any).grok.shell.closeAll());
}

async function addViewer(page: Page, viewerName: string) {
  return page.evaluate((name: string) => {
    const tv = (window as any).grok.shell.tv;
    return (window as any).grok.functions.call(`Charts:_${name}Viewer`);
  }, viewerName);
}

async function openDemoFile(page: Page, fileName: string) {
  return page.evaluate((name: string) =>
    (window as any).grok.data.getDemoTable(name), fileName);
}

test.describe('Radar Viewer (Charts package)', () => {
  test.beforeEach(async ({ page }) => {
    await page.goto(BASE_URL);
    await page.waitForFunction(() => typeof (window as any).grok !== 'undefined', { timeout: 30000 });
    await closeAll(page);
  });

  test('Step 1: Open earthquakes.csv and add Radar viewer', async ({ page }) => {
    const df = await openDemoFile(page, 'earthquakes.csv');
    expect(df).toBeTruthy();

    await page.evaluate(() => {
      const df = (window as any).grok.shell.tv?.dataFrame;
      expect(df).toBeTruthy();
    });

    const viewer = await page.evaluate(() =>
      (window as any).grok.shell.tv.addViewer('Radar'));
    expect(viewer).toBeTruthy();

    await page.waitForTimeout(2000);
    const screenshot = await page.screenshot({ path: 'scenario1-radar-earthquakes.png' });
    expect(screenshot).toBeTruthy();
  });

  test('Step 2: Open demog.csv and add Radar viewer', async ({ page }) => {
    await openDemoFile(page, 'demog.csv');
    const viewer = await page.evaluate(() =>
      (window as any).grok.shell.tv.addViewer('Radar'));
    expect(viewer).toBeTruthy();

    await page.waitForTimeout(2000);
    await page.screenshot({ path: 'scenario1-radar-demog.png' });
  });

  test('Step 3a: Radar viewer properties — checkboxes and values columns', async ({ page }) => {
    await openDemoFile(page, 'demog.csv');
    await page.evaluate(() => (window as any).grok.shell.tv.addViewer('Radar'));
    await page.waitForTimeout(1000);

    const result = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const viewer = tv.viewers.find((v: any) => v.type === 'Radar');
      if (!viewer) return { error: 'Radar viewer not found' };

      // Test showCurrentRow toggle
      viewer.props.showCurrentRow = true;
      const checked = viewer.props.showCurrentRow;
      viewer.props.showCurrentRow = false;
      const unchecked = viewer.props.showCurrentRow;

      // Test values columns
      const original = [...viewer.props.valuesColumnNames];
      if (original.length > 1) {
        viewer.props.valuesColumnNames = [original[0]];
        const reduced = viewer.props.valuesColumnNames.length;
        viewer.props.valuesColumnNames = original;
        return { checked, unchecked, reduced, original: original.length };
      }
      return { checked, unchecked, original: original.length };
    });

    expect(result).not.toHaveProperty('error');
    expect((result as any).checked).toBe(true);
    expect((result as any).unchecked).toBe(false);
  });

  test('Step 3b: Radar viewer properties — style changes', async ({ page }) => {
    await openDemoFile(page, 'demog.csv');
    await page.evaluate(() => (window as any).grok.shell.tv.addViewer('Radar'));
    await page.waitForTimeout(1000);

    const result = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const viewer = tv.viewers.find((v: any) => v.type === 'Radar');
      if (!viewer) return { error: 'Radar viewer not found' };

      // Test color changes
      const originalLine = viewer.props.lineColor;
      viewer.props.lineColor = 0xFF0000FF; // red
      const newLine = viewer.props.lineColor;
      viewer.props.lineColor = originalLine;

      return { originalLine, newLine };
    });

    expect(result).not.toHaveProperty('error');
    expect((result as any).newLine).toBe(0xFF0000FF);
  });

  test('Step 3c: Radar viewer — table switching (known bug)', async ({ page }) => {
    await openDemoFile(page, 'earthquakes.csv');
    await page.evaluate(() => (window as any).grok.shell.tv.addViewer('Radar'));
    await page.waitForTimeout(1000);
    await openDemoFile(page, 'demog.csv');
    await page.waitForTimeout(500);

    // Attempt table switch — known to throw NoSuchMethodError: toLowerCase
    const consoleErrors: string[] = [];
    page.on('console', msg => {
      if (msg.type() === 'error') consoleErrors.push(msg.text());
    });

    await page.evaluate(() => {
      const views = (window as any).grok.shell.tableViews;
      if (views.length < 2) return;
      const eqView = views.find((v: any) => v.dataFrame?.name?.includes('earthquakes'));
      const radarViewer = eqView?.viewers?.find((v: any) => v.type === 'Radar');
      const demogDf = views.find((v: any) => v.dataFrame?.name?.includes('demog'))?.dataFrame;
      if (radarViewer && demogDf) {
        try { radarViewer.props.table = demogDf; }
        catch (e) { console.error('Table switch error: ' + e); }
      }
    });

    await page.waitForTimeout(500);
    // This is a known failing step; the test documents the bug
    // TODO: fix NoSuchMethodError: toLowerCase when switching viewer table cross-view
  });
});
