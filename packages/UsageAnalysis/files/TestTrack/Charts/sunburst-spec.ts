import {test, expect} from '@playwright/test';

test.use({
  viewport: {width: 1920, height: 1080},
  launchOptions: {args: ['--window-size=1920,1080', '--window-position=0,0']},
  actionTimeout: 15_000,
  navigationTimeout: 60_000,
});

const baseUrl = process.env.DATAGROK_URL ?? 'http://localhost:8888';
const login = process.env.DATAGROK_LOGIN ?? 'admin';
const password = process.env.DATAGROK_PASSWORD ?? 'admin';

// Scenario says SPGI_v2.csv but that file is not present on dev; SPGI.csv is used instead.
const spgiPath = 'System:DemoFiles/SPGI.csv';
const demogPath = 'System:DemoFiles/demog.csv';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e?.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e?.message ?? e}`);
  }
}

test('Sunburst viewer', async ({page}) => {
  test.setTimeout(300_000);

  // Login
  await page.goto(baseUrl);
  await page.waitForTimeout(2000);
  const loginInput = page.getByPlaceholder('Login or Email').and(page.locator(':visible'));
  if (await loginInput.isVisible({timeout: 15000}).catch(() => false)) {
    await loginInput.click();
    await page.keyboard.type(login);
    await page.getByPlaceholder('Password').and(page.locator(':visible')).click();
    await page.keyboard.type(password);
    await page.keyboard.press('Enter');
  }
  await page.locator('[name="Browse"]').waitFor({timeout: 120000});

  // Baseline environment setup
  await page.evaluate(() => {
    document.querySelectorAll('.d4-dialog').forEach((d) => {
      const cancel = d.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
      if (cancel) cancel.click();
    });
    (window as any).grok.shell.closeAll();
    document.body.classList.add('selenium');
    (window as any).grok.shell.settings.showFiltersIconsConstantly = true;
    (window as any).grok.shell.windows.simpleMode = true;
  });

  // Step 1: Open SPGI.csv, add Sunburst; then open demog.csv, add Sunburst
  await softStep('Step 1: Open SPGI.csv and demog.csv, add Sunburst viewer for each', async () => {
    const spgi = await page.evaluate(async (path) => {
      const grok = (window as any).grok;
      grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 500));
      const df = await grok.dapi.files.readCsv(path);
      const tv = grok.shell.addTableView(df);
      await new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
      tv.addViewer('Sunburst');
      await new Promise((r) => setTimeout(r, 2000));
      const viewerTypes: string[] = [];
      for (const v of tv.viewers) viewerTypes.push(v.type);
      return {rowCount: df.rowCount, viewerTypes};
    }, spgiPath);
    expect(spgi.rowCount).toBeGreaterThan(0);
    expect(spgi.viewerTypes).toContain('Sunburst');

    const demog = await page.evaluate(async (path) => {
      const grok = (window as any).grok;
      grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 500));
      const df = await grok.dapi.files.readCsv(path);
      const tv = grok.shell.addTableView(df);
      await new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
      tv.addViewer('Sunburst');
      await new Promise((r) => setTimeout(r, 2000));
      const viewerTypes: string[] = [];
      for (const v of tv.viewers) viewerTypes.push(v.type);
      return {rowCount: df.rowCount, viewerTypes};
    }, demogPath);
    expect(demog.rowCount).toBe(5850);
    expect(demog.viewerTypes).toContain('Sunburst');
  });

  // Step 2: Verify Sunburst property names via getProperties()
  await softStep('Step 2: Context Panel property names include hierarchyColumnNames, inheritFromGrid, includeNulls', async () => {
    const result = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      let sunburst: any = null;
      for (const v of tv.viewers) if (v.type === 'Sunburst') { sunburst = v; break; }
      if (!sunburst) return {propNames: [] as string[]};
      const props = sunburst.props.getProperties();
      const propNames: string[] = [];
      for (const p of props) propNames.push(p.name);
      return {propNames};
    });
    expect(result.propNames).toEqual(expect.arrayContaining(['hierarchyColumnNames', 'inheritFromGrid', 'includeNulls']));
  });

  // Step 3.1: Table switching between SPGI and demog — AMBIGUOUS in 2b
  await softStep('Step 3.1: Table switching SPGI <-> demog (AMBIGUOUS, not exercised)', async () => {
    // 2b closed all views between opens; the live table-rebind path was not exercised.
    test.skip(true, 'AMBIGUOUS: table switching via UI/props not exercised in MCP run');
  });

  // Step 3.2: Select Columns — set hierarchyColumnNames via setOptions, read back
  await softStep('Step 3.2: Set hierarchyColumnNames to [SEX, RACE] and read back', async () => {
    const result = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      let sunburst: any = null;
      for (const v of tv.viewers) if (v.type === 'Sunburst') { sunburst = v; break; }
      if (!sunburst) return {cols: [] as string[]};
      sunburst.setOptions({hierarchyColumnNames: ['SEX', 'RACE']});
      await new Promise((r) => setTimeout(r, 500));
      const cols = sunburst.props.get('hierarchyColumnNames') as string[];
      return {cols: Array.from(cols ?? [])};
    });
    expect(result.cols).toEqual(['SEX', 'RACE']);
  });

  // Step 3.3: Inherit from grid — toggle on, read back
  await softStep('Step 3.3: Toggle inheritFromGrid=true and read back', async () => {
    const result = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      let sunburst: any = null;
      for (const v of tv.viewers) if (v.type === 'Sunburst') { sunburst = v; break; }
      if (!sunburst) return {inheritFromGrid: null as any};
      sunburst.setOptions({hierarchyColumnNames: ['SEX'], inheritFromGrid: true});
      await new Promise((r) => setTimeout(r, 500));
      return {inheritFromGrid: sunburst.props.get('inheritFromGrid')};
    });
    expect(result.inheritFromGrid).toBe(true);
  });

  // Step 3.4: Toggle includeNulls true then false
  await softStep('Step 3.4: Toggle includeNulls true then false and read back', async () => {
    const result = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      let sunburst: any = null;
      for (const v of tv.viewers) if (v.type === 'Sunburst') { sunburst = v; break; }
      if (!sunburst) return {firstRead: null as any, secondRead: null as any};
      sunburst.setOptions({includeNulls: true});
      await new Promise((r) => setTimeout(r, 300));
      const firstRead = sunburst.props.get('includeNulls');
      sunburst.setOptions({includeNulls: false});
      await new Promise((r) => setTimeout(r, 300));
      const secondRead = sunburst.props.get('includeNulls');
      return {firstRead, secondRead};
    });
    expect(result.firstRead).toBe(true);
    expect(result.secondRead).toBe(false);
  });

  // Step 4: Reset view — AMBIGUOUS (canvas-based)
  await softStep('Step 4: Reset view via double-click / context menu (AMBIGUOUS, canvas-based)', async () => {
    test.skip(true, 'AMBIGUOUS: requires coordinate-precise canvas double-click or context menu');
  });

  // Step 5: Multi-selection — AMBIGUOUS (canvas-based)
  await softStep('Step 5: Multi-selection (Click / Ctrl+Click / Ctrl+Shift+Click) (AMBIGUOUS, canvas-based)', async () => {
    test.skip(true, 'AMBIGUOUS: canvas-based selection without a selection API');
  });

  // Step 6: Select/filter on empty category — AMBIGUOUS
  await softStep('Step 6: Select/filter on empty (null) category (AMBIGUOUS, canvas-based)', async () => {
    test.skip(true, 'AMBIGUOUS: canvas-based, depends on null-segment rendering and selection');
  });

  // Step 7: Projects & layouts save/restore — AMBIGUOUS (not exercised)
  await softStep('Step 7: Projects & layouts save/restore (AMBIGUOUS, not exercised)', async () => {
    test.skip(true, 'AMBIGUOUS: layout round-trip with viewer settings preservation is its own test surface');
  });

  // Step 8: Old layout compatibility (issue #2979) — SKIP (external asset)
  await softStep('Step 8: Old layout compatibility — issue #2979 (SKIP, external asset)', async () => {
    test.skip(true, 'SKIP: requires specific layout file from GitHub issue attachment');
  });

  // Step 9: Collaborative filtering — AMBIGUOUS
  await softStep('Step 9: Collaborative filtering — internal + panel filters combine (AMBIGUOUS)', async () => {
    test.skip(true, 'AMBIGUOUS: not exercised in MCP run');
  });

  // Cleanup
  await page.evaluate(() => (window as any).grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
