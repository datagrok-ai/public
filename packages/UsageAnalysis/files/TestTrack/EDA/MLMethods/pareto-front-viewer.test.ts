import { test, expect } from './helpers';
import { openDemoCsv, resetShell } from './helpers';

const PARETO = 'Pareto front';

async function addParetoViewerToCurrentTable(page: import('@playwright/test').Page): Promise<void> {

  await page.evaluate((type) => {
    const g = (window as unknown as { grok: any }).grok;
    g.shell.tv.addViewer(type);
  }, PARETO);
  await page.waitForTimeout(2000);
}

async function readParetoLabelColumns(page: import('@playwright/test').Page): Promise<string[]> {
  return page.evaluate(() => {
    const g = (window as unknown as { grok: any }).grok;
    const v = Array.from(g.shell.tv.viewers).find((x: any) => x?.type?.toLowerCase().includes('pareto')) as any;
    const props: any[] = v?.props?.getProperties?.() ?? [];
    const labelProp = props.find((p) => /^labelColumns/i.test(p.name));
    if (labelProp && typeof labelProp.get === 'function') return (labelProp.get(v) ?? []) as string[];
    try { return (v?.props?.labelColumnsColumnNames ?? []) as string[]; }
    catch { return [] as string[]; }
  });
}

async function waitForParetoLabelColumns(
  page: import('@playwright/test').Page, expected: string[], timeoutMs = 15_000,
): Promise<void> {
  await page.waitForFunction((want: string[]) => {
    const g = (window as unknown as { grok: any }).grok;
    const v = Array.from(g.shell.tv.viewers).find((x: any) => x?.type?.toLowerCase().includes('pareto')) as any;
    const props: any[] = v?.props?.getProperties?.() ?? [];
    const labelProp = props.find((p) => /^labelColumns/i.test(p.name));
    let val: string[] = [];
    if (labelProp && typeof labelProp.get === 'function') val = (labelProp.get(v) ?? []) as string[];
    else { try { val = (v?.props?.labelColumnsColumnNames ?? []) as string[]; } catch { val = []; } }
    return val.length === want.length && want.every((n, i) => val[i] === n);
  }, expected, { timeout: timeoutMs });
}

async function openCarsWithMissingFromCarsCsv(page: import('@playwright/test').Page): Promise<void> {
  await openDemoCsv(page, 'cars.csv');
  await page.evaluate(() => {
    const g = (window as unknown as { grok: any }).grok;
    const df = g.shell.tv.dataFrame;
    const turbo = df.col('turbo');
    if (turbo) {
      for (let i = 0; i < df.rowCount; i++) turbo.set(i, null);
    }
    df.name = 'cars-with-missing';
  });
}

test.describe.serial('EDA / Pareto Front Viewer', () => {
  test.afterEach(async ({ page }) => { await resetShell(page); });

  test('cars-with-missing: empty and string columns are excluded from Minimize/Maximize', async ({ page }) => {
    test.setTimeout(180_000);

    await openCarsWithMissingFromCarsCsv(page);
    await addParetoViewerToCurrentTable(page);

    const { eligible, allNames } = await page.evaluate(() => {
      const g = (window as unknown as { grok: any }).grok;
      const v = Array.from(g.shell.tv.viewers).find((x: any) => x?.type?.toLowerCase().includes('pareto'));
      if (!v) return { eligible: [] as string[], allNames: [] as string[] };
      const df = g.shell.tv.dataFrame;
      const all = df.columns.names() as string[];
      const props: any[] = (v as any).props?.getProperties?.() ?? [];
      const minProp = props.find((p) => /minimi[sz]e/i.test(p.name));
      const choices: string[] = (minProp?.choices ?? []).map((c: any) => String(c));
      return { eligible: choices, allNames: all };
    });

    expect(eligible).not.toContain('turbo');
    expect(eligible).not.toContain('model');

    expect(allNames).toEqual(expect.arrayContaining(['turbo', 'model']));

  });

  test('cars.csv auto-selects "model" as Label (unique values)', async ({ page }) => {
    test.setTimeout(120_000);

    await openDemoCsv(page, 'cars.csv');
    await addParetoViewerToCurrentTable(page);

    await waitForParetoLabelColumns(page, ['model']);
    expect(await readParetoLabelColumns(page)).toEqual(['model']);
  });

  test('demog.csv auto-selects "USUBJID" as Label (unique values)', async ({ page }) => {
    test.setTimeout(120_000);

    await openDemoCsv(page, 'demog.csv');
    await addParetoViewerToCurrentTable(page);

    await waitForParetoLabelColumns(page, ['USUBJID']);
    expect(await readParetoLabelColumns(page)).toEqual(['USUBJID']);
  });

  test('Viewer exposes Description/Objectives/Axes/Labels/Legend categories without errors', async ({ page }) => {
    test.setTimeout(120_000);

    await openDemoCsv(page, 'cars.csv');
    await addParetoViewerToCurrentTable(page);

    const expectedCategories = ['Description', 'Objectives', 'Axes', 'Labels', 'Legend'];
    await page.waitForFunction((want: string[]) => {
      const g = (window as unknown as { grok: any }).grok;
      const v = Array.from(g.shell.tv.viewers).find((x: any) => x?.type?.toLowerCase().includes('pareto')) as any;
      const props: any[] = v?.props?.getProperties?.() ?? [];
      const cats = new Set<string>(props.map((p) => String(p.category ?? '')).filter((c) => c.length > 0));
      return want.every((c) => cats.has(c));
    }, expectedCategories, { timeout: 15_000 });

    const categories = await page.evaluate(() => {
      const g = (window as unknown as { grok: any }).grok;
      const v = Array.from(g.shell.tv.viewers).find((x: any) => x?.type?.toLowerCase().includes('pareto')) as any;
      const props: any[] = v?.props?.getProperties?.() ?? [];
      return Array.from(new Set(props.map((p) => p.category ?? '').filter((c: string) => c.length > 0)));
    });
    expect(categories).toEqual(expect.arrayContaining(expectedCategories));
  });
});
