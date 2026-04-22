import {test, expect} from '@playwright/test';

test.use({
  viewport: {width: 1920, height: 1080},
  launchOptions: {args: ['--window-size=1920,1080', '--window-position=0,0']},
  actionTimeout: 15_000,
  navigationTimeout: 60_000,
});

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try { await test.step(name, fn); }
  catch (e: any) { stepErrors.push({step: name, error: e.message ?? String(e)}); }
}

test('Pareto Front viewer scenario', async ({page}) => {
  test.setTimeout(300_000);

  const baseUrl = process.env.DATAGROK_URL ?? 'http://localhost:8888';
  const login = process.env.DATAGROK_LOGIN ?? 'admin';
  const password = process.env.DATAGROK_PASSWORD ?? 'admin';

  await page.goto(baseUrl);
  // Give the Dart app a moment to boot before checking the login input.
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

  // Baseline setup: selenium class, Tabs mode, close everything.
  await page.evaluate(async () => {
    const g: any = (window as any).grok;
    document.body.classList.add('selenium');
    g.shell.settings.showFiltersIconsConstantly = true;
    g.shell.windows.simpleMode = true;
    g.shell.closeAll();
  });

  await softStep('Step 1 — Open cars-with-missing.csv from Demo Files', async () => {
    // 2b observed: readCsv returned a bare 502 and cars-with-missing.csv is not present
    // under System:DemoFiles on dev (only cars.csv exists). Record that verbatim.
    const result = await page.evaluate(async () => {
      const g: any = (window as any).grok;
      const out: {opened: boolean; error?: string; listed?: string[]} = {opened: false};
      try {
        const df = await g.dapi.files.readCsv('System:DemoFiles/cars-with-missing.csv');
        g.shell.addTableView(df);
        out.opened = true;
      } catch (e: any) {
        out.error = e?.message ?? String(e);
      }
      try {
        const items = await g.dapi.files.list('System:DemoFiles/', true);
        out.listed = items
          .map((f: any) => (f?.fileName ?? f?.name ?? String(f)))
          .filter((n: string) => typeof n === 'string' && n.toLowerCase().includes('car'));
      } catch { /* ignore secondary listing failure */ }
      return out;
    });
    if (result.opened)
      throw new Error('cars-with-missing.csv unexpectedly opened — scenario expected a missing file on dev.');
    throw new Error(
      'cars-with-missing.csv not available on dev (readCsv error: ' +
      (result.error ?? 'unknown') + '). System:DemoFiles car-related files: ' +
      JSON.stringify(result.listed ?? []));
  });

  await softStep('Step 2 — Add Pareto front viewer, open Properties panel', async () => {
    throw new Error('SKIP: blocked by Step 1 (cars-with-missing.csv is missing on dev).');
  });

  await softStep('Step 3 — Minimize/Maximize lists exclude empty & string cols', async () => {
    throw new Error('SKIP: blocked by Step 1 (cars-with-missing.csv is missing on dev).');
  });

  await softStep('Step 4 — Select all columns in Maximize; expect warning', async () => {
    throw new Error('SKIP: blocked by Step 1 (cars-with-missing.csv is missing on dev).');
  });

  await softStep('Step 5 — Open cars.csv; add Pareto Front; "model" auto-selected as Label', async () => {
    const info = await page.evaluate(async () => {
      const g: any = (window as any).grok;
      g.shell.closeAll();
      const df = await g.dapi.files.readCsv('System:DemoFiles/cars.csv');
      const tv = g.shell.addTableView(df);
      await new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
      const paretoV = tv.addViewer('Pareto Front');
      // Allow a brief settle for viewer props to populate.
      await new Promise((r) => setTimeout(r, 1500));
      return {
        type: paretoV.type,
        labelColumnsColumnNames: paretoV.props.labelColumnsColumnNames,
        minimizeColumnNames: paretoV.props.minimizeColumnNames,
        maximizeColumnNames: paretoV.props.maximizeColumnNames,
        xAxisColumnName: paretoV.props.xAxisColumnName,
        yAxisColumnName: paretoV.props.yAxisColumnName,
        autoLabelsSelection: paretoV.props.autoLabelsSelection,
      };
    });

    // Scenario expectation: "model" column automatically selected as Label on cars.csv.
    expect(info.labelColumnsColumnNames).toContain('model');
  });

  await softStep('Step 6 — Open demog; Label empty by default OR unique-value category auto-selects', async () => {
    const info = await page.evaluate(async () => {
      const g: any = (window as any).grok;
      g.shell.closeAll();
      const df = await g.dapi.files.readCsv('System:DemoFiles/demog.csv');
      const tv = g.shell.addTableView(df);
      await new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
      const paretoV = tv.addViewer('Pareto Front');
      await new Promise((r) => setTimeout(r, 1500));

      const labelCols: string[] = paretoV.props.labelColumnsColumnNames ?? [];
      const rowCount: number = df.rowCount;
      const uniqueCounts: Record<string, number> = {};
      for (const name of labelCols) {
        const col = df.columns.byName(name);
        if (!col) continue;
        const seen = new Set<any>();
        for (let i = 0; i < col.length; i++) seen.add(col.get(i));
        uniqueCounts[name] = seen.size;
      }
      return {labelCols, rowCount, uniqueCounts};
    });

    // Scenario allows: (a) empty by default, OR (b) every auto-selected column has unique values.
    if (info.labelCols.length === 0)
      return;
    for (const name of info.labelCols) {
      const u = info.uniqueCounts[name] ?? 0;
      expect(u, `Auto-selected label column "${name}" must have unique values (got ${u}/${info.rowCount})`)
        .toBe(info.rowCount);
    }
  });

  await softStep('Step 7 — Review viewer properties (Description, Objectives, Axes, Labels, Legend)', async () => {
    const props = await page.evaluate(async () => {
      const g: any = (window as any).grok;
      g.shell.closeAll();
      const df = await g.dapi.files.readCsv('System:DemoFiles/cars.csv');
      const tv = g.shell.addTableView(df);
      await new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
      const paretoV = tv.addViewer('Pareto Front');
      await new Promise((r) => setTimeout(r, 1500));
      const allProps = paretoV.props.getProperties() as Array<any>;
      const categories = new Set<string>(
        allProps.map((p: any) => p?.category ?? p?._category ?? '').filter((c: string) => !!c));
      return {categories: Array.from(categories), count: allProps.length};
    });

    expect(props.count).toBeGreaterThan(0);
    const expected = ['Description', 'Objectives', 'Axes', 'Labels', 'Legend'];
    for (const cat of expected)
      expect(props.categories, `Missing property category "${cat}" (got ${JSON.stringify(props.categories)})`)
        .toContain(cat);
  });

  if (stepErrors.length > 0)
    throw new Error('Step failures:\n' + stepErrors.map((e) => `- ${e.step}: ${e.error}`).join('\n'));
});
