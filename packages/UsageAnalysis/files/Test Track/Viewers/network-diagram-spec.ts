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

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try { await test.step(name, fn); }
  catch (e: any) { stepErrors.push({step: name, error: e?.message ?? String(e)}); }
}

test('Network diagram', async ({page}) => {
  test.setTimeout(300_000);

  await page.goto(baseUrl);
  const loginInput = page.getByPlaceholder('Login or Email').and(page.locator(':visible'));
  if (await loginInput.isVisible({timeout: 15000}).catch(() => false)) {
    await loginInput.click();
    await page.keyboard.type(login);
    await page.getByPlaceholder('Password').and(page.locator(':visible')).click();
    await page.keyboard.type(password);
    await page.keyboard.press('Enter');
  }
  await page.locator('[name="Browse"]').waitFor({timeout: 120000});

  await page.evaluate(async () => {
    const w = window as any;
    document.body.classList.add('selenium');
    w.grok.shell.settings.showFiltersIconsConstantly = true;
    w.grok.shell.windows.simpleMode = true;
    w.grok.shell.closeAll();
    const df = await w.grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
    const tv = w.grok.shell.addTableView(df);
    await new Promise((resolve: any) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
    tv.getFiltersGroup();
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  await softStep('Step 2: Add Network diagram via Toolbox icon', async () => {
    await page.locator('[name="icon-network-diagram"]').click();
    await page.locator('[name="viewer-Network-diagram"]').waitFor({timeout: 15000});
    await page.waitForTimeout(1500);
    const info = await page.evaluate(() => {
      const w = window as any;
      const nd = w.grok.shell.tv.viewers.find((v: any) => v.type === 'Network diagram');
      return {added: !!nd, node1: nd?.props.node1ColumnName, node2: nd?.props.node2ColumnName};
    });
    expect(info.added).toBe(true);
    expect(info.node1).toBeTruthy();
    expect(info.node2).toBeTruthy();
  });

  await softStep('Step 3: Switch Node 1 to RACE, Node 2 to DEMOG', async () => {
    const info = await page.evaluate(async () => {
      const w = window as any;
      const nd = w.grok.shell.tv.viewers.find((v: any) => v.type === 'Network diagram');
      nd.props.node1ColumnName = 'RACE';
      nd.props.node2ColumnName = 'DEMOG';
      await new Promise((r: any) => setTimeout(r, 800));
      return {node1: nd.props.node1ColumnName, node2: nd.props.node2ColumnName};
    });
    expect(info.node1).toBe('RACE');
    expect(info.node2).toBe('DEMOG');
  });

  // Steps 4-7: canvas click/shift+click/ctrl+click/dblclick on vis.js-rendered nodes — omitted
  // because the viewer renders into canvas (no DOM handles for nodes/edges).

  await softStep('Step 8: Open Property Pane via Gear icon', async () => {
    await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Network-diagram"]')!;
      const panel = root.closest('.panel-base')!;
      (panel.querySelector('[name="icon-font-icon-settings"]') as HTMLElement).click();
    });
    await page.waitForTimeout(800);
    const hasData = await page.evaluate(() => {
      return !!Array.from(document.querySelectorAll('.property-grid-category'))
        .find(e => e.textContent?.trim() === 'Data');
    });
    expect(hasData).toBe(true);
  });

  await softStep('Step 9: Set edge/node color/size/width columns', async () => {
    const after = await page.evaluate(async () => {
      const w = window as any;
      const nd = w.grok.shell.tv.viewers.find((v: any) => v.type === 'Network diagram');
      nd.props.edgeColorColumnName = 'AGE';
      nd.props.edgeColorAggrType = 'avg';
      nd.props.edgeWidthColumnName = 'WEIGHT';
      nd.props.edgeWidthAggrType = 'avg';
      nd.props.node1SizeColumnName = 'AGE';
      nd.props.node1ColorColumnName = 'SEX';
      await new Promise((r: any) => setTimeout(r, 1000));
      return {
        edgeColor: nd.props.edgeColorColumnName,
        edgeWidth: nd.props.edgeWidthColumnName,
        node1Size: nd.props.node1SizeColumnName,
        node1Color: nd.props.node1ColorColumnName,
      };
    });
    expect(after.edgeColor).toBe('AGE');
    expect(after.edgeWidth).toBe('WEIGHT');
    expect(after.node1Size).toBe('AGE');
    expect(after.node1Color).toBe('SEX');
  });

  await softStep('Step 10: Style toggles (selectors, arrows, simulation)', async () => {
    const res = await page.evaluate(async () => {
      const w = window as any;
      const nd = w.grok.shell.tv.viewers.find((v: any) => v.type === 'Network diagram');
      nd.props.showColumnSelectors = false;
      await new Promise((r: any) => setTimeout(r, 400));
      const hiddenCount = Array.from(document.querySelectorAll('[name="viewer-Network-diagram"] [name^="div-column-combobox"]'))
        .filter(e => !!(e as HTMLElement).offsetParent).length;
      nd.props.showColumnSelectors = true;
      await new Promise((r: any) => setTimeout(r, 400));
      const shownCount = Array.from(document.querySelectorAll('[name="viewer-Network-diagram"] [name^="div-column-combobox"]'))
        .filter(e => !!(e as HTMLElement).offsetParent).length;
      nd.props.showArrows = 'to';
      await new Promise((r: any) => setTimeout(r, 400));
      nd.props.suspendSimulation = true;
      await new Promise((r: any) => setTimeout(r, 400));
      const suspT = nd.props.suspendSimulation;
      nd.props.suspendSimulation = false;
      return {hiddenCount, shownCount, arrows: nd.props.showArrows, suspT, suspF: nd.props.suspendSimulation};
    });
    expect(res.hiddenCount).toBe(0);
    expect(res.shownCount).toBe(2);
    expect(res.arrows).toBe('to');
    expect(res.suspT).toBe(true);
    expect(res.suspF).toBe(false);
  });

  await softStep('Step 11: Filter AGE > 40 and toggle Show Filtered Out Nodes', async () => {
    const res = await page.evaluate(async () => {
      const w = window as any;
      const tv = w.grok.shell.tv;
      const df = tv.dataFrame;
      const bs = df.filter;
      const age = df.col('AGE');
      for (let i = 0; i < df.rowCount; i++) bs.set(i, age.get(i) > 40, false);
      bs.fireChanged();
      await new Promise((r: any) => setTimeout(r, 800));
      const filteredCount = bs.trueCount;
      const nd = tv.viewers.find((v: any) => v.type === 'Network diagram');
      nd.props.showFilteredOutNodes = true;
      await new Promise((r: any) => setTimeout(r, 500));
      return {filteredCount, total: df.rowCount, showFilteredOut: nd.props.showFilteredOutNodes};
    });
    expect(res.filteredCount).toBeGreaterThan(0);
    expect(res.filteredCount).toBeLessThan(res.total);
    expect(res.showFilteredOut).toBe(true);
  });

  await softStep('Step 12: Close viewer via × icon', async () => {
    await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Network-diagram"]')!;
      const panel = root.closest('.panel-base')!;
      (panel.querySelector('[name="Close"]') as HTMLElement).click();
    });
    await page.waitForTimeout(800);
    const gone = await page.evaluate(() => {
      const w = window as any;
      return !w.grok.shell.tv.viewers.find((v: any) => v.type === 'Network diagram');
    });
    expect(gone).toBe(true);
  });

  if (stepErrors.length > 0)
    throw new Error('Step failures:\n' + stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n'));
});
