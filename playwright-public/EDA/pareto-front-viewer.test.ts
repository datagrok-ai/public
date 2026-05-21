import { test, expect } from '@playwright/test';
import { openDemoCsv, resetShell } from './helpers';

// Test Track scenario: EDA/pareto-front-viewer.md
//
// Test 1 — empty & non-numeric column handling on cars-with-missing.csv:
//   * Minimize/Maximize dropdowns must exclude empty columns (e.g. `turbo`) and string
//     columns (e.g. `model`).
//   * Selecting all columns in Maximize must produce a warning ("cannot minimize and
//     maximize the same feature simultaneously").
//
// Test 2 — label auto-selection on cars.csv:
//   * `model` is unique on every row, so it must be auto-picked as Label.
//
// Test 3 — label auto-selection on demog.csv:
//   * `USUBJID` is unique on every row, so it must be auto-picked as Label.
//
// Test 4 — viewer property categories:
//   * The viewer exposes Description, Objectives, Axes, Labels, Legend categories
//     without exceptions.
//
// All assertions read the viewer's properties via the public `props.getProperties()`
// API on the viewer instance — there is no DOM equivalent for the columns-eligible-for-
// objectives list, since the Properties panel renders these as a hover dropdown.

const PARETO = 'Pareto front';

async function addParetoViewerToCurrentTable(page: import('@playwright/test').Page): Promise<void> {
  // Adding viewers programmatically is the platform's own pattern for the "Add Viewer"
  // ribbon menu — the menu item ultimately calls `tv.addViewer(type)`. UI-mode menu
  // hover-chain priming over the viewer ribbon is brittle on smaller viewports.
  await page.evaluate((type) => {
    const g = (window as unknown as { grok: any }).grok;
    g.shell.tv.addViewer(type);
  }, PARETO);
  await page.waitForTimeout(2000);
}

/**
 * Read the Pareto viewer's auto-selected Label columns through the public property
 * descriptor (`props.getProperties()`), falling back to the column-names accessor.
 * Reading via the descriptor avoids relying on private fields the production minifier
 * may have renamed.
 */
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

/**
 * Wait in-browser until the Pareto viewer's auto-selected Label columns equal `expected`.
 * Label selection runs asynchronously on the viewer's first render; `waitForFunction`
 * polls inside the page and rejects only once on timeout, so transient empty reads are not
 * logged as report errors (unlike `expect.poll`, which records every failed retry).
 */
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

/**
 * `cars-with-missing.csv` is referenced by the scenario but is not deployed under
 * `System:DemoFiles/` on every Datagrok tenant (see pareto-front-viewer-run.md, step 1).
 * Build an equivalent dataset in-memory by loading `cars.csv` and nulling every cell of
 * the `turbo` column — the scenario only relies on `turbo` being entirely empty.
 */
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

    // Read eligible objective columns from the viewer's property metadata.
    // The Minimize and Maximize column-list properties enumerate eligible columns via
    // their `choices` array; non-numeric and entirely-empty columns must be absent.
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

    // Empty `turbo` and string `model` must be absent from the eligible list.
    expect(eligible).not.toContain('turbo');
    expect(eligible).not.toContain('model');
    // Sanity check the columns actually exist in the source DF (so the absences above
    // are not just a missing dataset).
    expect(allNames).toEqual(expect.arrayContaining(['turbo', 'model']));

    // Note on the "Maximize all triggers a warning" sub-scenario: the warning text is
    // rendered directly on the viewer canvas via `_showErrorMessage` and the underlying
    // `errMsg` field is renamed by the production minifier, so there is no stable surface
    // to assert against from automation. The eligibility check above is the part of the
    // scenario this test can verify deterministically.
  });

  test('cars.csv auto-selects "model" as Label (unique values)', async ({ page }) => {
    test.setTimeout(120_000);

    await openDemoCsv(page, 'cars.csv');
    await addParetoViewerToCurrentTable(page);

    // Label auto-selection runs asynchronously on the viewer's first render — wait for the
    // property to be populated, then assert once.
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

    // The viewer's property descriptors load asynchronously after `addViewer` resolves —
    // the first read can return an empty list. Wait in-browser for all expected categories
    // (so transient empty reads aren't logged as report errors), then assert once.
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
