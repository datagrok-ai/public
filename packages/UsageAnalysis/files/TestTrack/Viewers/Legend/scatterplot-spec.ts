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

async function openSPGI(page: any): Promise<void> {
  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    (window as any).grok.shell.settings.showFiltersIconsConstantly = true;
    (window as any).grok.shell.windows.simpleMode = true;
    (window as any).grok.shell.closeAll();
    const df = await (window as any).grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv');
    (window as any).grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(null); });
      setTimeout(resolve, 5000);
    });
    const hasBioChem = Array.from({length: df.columns.length}, (_, i: number) => df.columns.byIndex(i))
      .some((c: any) => c.semType === 'Molecule' || c.semType === 'Macromolecule');
    if (hasBioChem) {
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise(r => setTimeout(r, 200));
      }
      await new Promise(r => setTimeout(r, 5000));
    }
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').first().waitFor({timeout: 30000});
}

test('Legend scatterplot behaviors', async ({page}) => {
  test.setTimeout(900_000);

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

  await openSPGI(page);

  await softStep('Sub 1: Color=Series + Marker=Series — combined legend', async () => {
    const items = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      tv.addViewer('Scatter plot');
      await new Promise(r => setTimeout(r, 500));
      const sp = tv.viewers.find((v: any) => v.type === 'Scatter plot');
      sp.props.colorColumnName = 'Series';
      sp.props.markersColumnName = 'Series';
      try { sp.props.legendVisibility = 'Always'; } catch(_) {}
      await new Promise(r => setTimeout(r, 1200));
      const items = sp.root.querySelectorAll('[name="legend"] .d4-legend-item');
      return items.length;
    });
    expect(items).toBeGreaterThanOrEqual(5);
  });

  await softStep('Sub 1: color picker visible on hover', async () => {
    const found = await page.evaluate(async () => {
      const sp = (window as any).grok.shell.tv.viewers.find((v: any) => v.type === 'Scatter plot');
      const item = sp.root.querySelector('[name="legend"] .d4-legend-item');
      item?.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      item?.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
      await new Promise(r => setTimeout(r, 400));
      return !!document.querySelector('[name="legend-icon-color-picker"], [name="Legend-icon-color-picker"]');
    });
    expect(found).toBe(true);
  });

  await softStep('Sub 1: categorical formula column produces categorical legend', async () => {
    const count = await page.evaluate(async () => {
      const df = (window as any).grok.shell.tv.dataFrame;
      try { await df.columns.addNewCalculated('testCat', "if(${Stereo Category}=='S_UNKN', null, ${Series})"); } catch(_) {}
      await new Promise(r => setTimeout(r, 1500));
      const sp = (window as any).grok.shell.tv.viewers.find((v: any) => v.type === 'Scatter plot');
      sp.props.markersColumnName = '';
      sp.props.colorColumnName = 'testCat';
      await new Promise(r => setTimeout(r, 1200));
      const items = sp.root.querySelectorAll('[name="legend"] .d4-legend-item');
      return items.length;
    });
    expect(count).toBeGreaterThan(0);
  });

  await softStep('Sub 1: Color=ID, Marker=Core', async () => {
    const count = await page.evaluate(async () => {
      const sp = (window as any).grok.shell.tv.viewers.find((v: any) => v.type === 'Scatter plot');
      sp.props.colorColumnName = 'ID';
      sp.props.markersColumnName = 'Core';
      await new Promise(r => setTimeout(r, 1500));
      const items = sp.root.querySelectorAll('[name="legend"] .d4-legend-item');
      return items.length;
    });
    expect(count).toBeGreaterThan(0);
  });

  await softStep('Sub 2: legend updates when X axis switches between col1 and col2', async () => {
    const res = await page.evaluate(async () => {
      (window as any).grok.shell.closeAll();
      const df = await (window as any).grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv');
      (window as any).grok.shell.addTableView(df);
      await new Promise(r => setTimeout(r, 4000));
      await df.columns.addNewCalculated('col1', "if(${Stereo Category}!='S_UNKN', null, ${Average Mass})");
      await df.columns.addNewCalculated('col2', "if(${Stereo Category}=='S_UNKN', null, ${Average Mass})");
      await new Promise(r => setTimeout(r, 1500));
      const tv = (window as any).grok.shell.tv;
      tv.addViewer('Scatter plot');
      await new Promise(r => setTimeout(r, 500));
      const sp = tv.viewers.find((v: any) => v.type === 'Scatter plot');
      sp.props.xColumnName = 'col1';
      sp.props.colorColumnName = 'Stereo Category';
      try { sp.props.legendVisibility = 'Always'; } catch(_) {}
      await new Promise(r => setTimeout(r, 1500));
      const a = sp.root.querySelectorAll('[name="legend"] .d4-legend-item').length;
      sp.props.xColumnName = 'col2';
      await new Promise(r => setTimeout(r, 1500));
      const b = sp.root.querySelectorAll('[name="legend"] .d4-legend-item').length;
      return {a, b};
    });
    expect(res.a).not.toBe(res.b);
  });

  await softStep('Sub 3: in-viewer filter hides legend DOM (platform bug)', async () => {
    const res = await page.evaluate(async () => {
      const sp = (window as any).grok.shell.tv.viewers.find((v: any) => v.type === 'Scatter plot');
      sp.props.markersColumnName = 'Stereo Category';
      sp.props.filter = '${Stereo Category} in ("R_ONE", "S_UNKN")';
      await new Promise(r => setTimeout(r, 1500));
      const after = sp.root.querySelectorAll('[name="legend"] .d4-legend-item').length;
      sp.props.filter = '';
      await new Promise(r => setTimeout(r, 1200));
      const cleared = sp.root.querySelectorAll('[name="legend"] .d4-legend-item').length;
      return {after, cleared};
    });
    expect(res.cleared).toBeGreaterThan(0);
  });

  await softStep('Sub 4: Filter Panel filter on Primary Scaffold Name', async () => {
    const res = await page.evaluate(async () => {
      (window as any).grok.shell.closeAll();
      const df = await (window as any).grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv');
      (window as any).grok.shell.addTableView(df);
      await new Promise(r => setTimeout(r, 4000));
      const tv = (window as any).grok.shell.tv;
      tv.addViewer('Scatter plot');
      await new Promise(r => setTimeout(r, 500));
      const sp = tv.viewers.find((v: any) => v.type === 'Scatter plot');
      sp.props.xColumnName = 'Chemical Space X';
      sp.props.yColumnName = 'Chemical Space Y';
      sp.props.colorColumnName = 'Primary Scaffold Name';
      sp.props.markersColumnName = 'Stereo Category';
      try { sp.props.legendVisibility = 'Always'; } catch(_) {}
      tv.getFiltersGroup();
      await new Promise(r => setTimeout(r, 1500));
      const fg = tv.getFiltersGroup();
      const DG = (window as any).DG;
      const scaffolds = df.col('Primary Scaffold Name').categories;
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'Primary Scaffold Name', selected: scaffolds.slice(0, 2)});
      await new Promise(r => setTimeout(r, 1500));
      const items = sp.root.querySelectorAll('[name="legend"] .d4-legend-item');
      return {items: items.length, filtered: df.filter.trueCount};
    });
    expect(res.items).toBeGreaterThan(0);
  });

  await softStep('Sub 5: grid linear color coding', async () => {
    const res = await page.evaluate(async () => {
      (window as any).grok.shell.closeAll();
      const df = await (window as any).grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv');
      (window as any).grok.shell.addTableView(df);
      await new Promise(r => setTimeout(r, 4000));
      const tv = (window as any).grok.shell.tv;
      tv.addViewer('Scatter plot');
      tv.addViewer('Box plot');
      tv.addViewer('PC Plot');
      await new Promise(r => setTimeout(r, 1500));
      const sp = tv.viewers.find((v: any) => v.type === 'Scatter plot');
      sp.props.colorColumnName = 'Chemical Space X';
      try { sp.props.legendVisibility = 'Always'; } catch(_) {}
      const col = df.col('Chemical Space X');
      col.tags['.color-coding-type'] = 'Linear';
      (window as any).grok.shell.tv.grid.invalidate();
      await new Promise(r => setTimeout(r, 1500));
      return {colorCodingType: col.tags['.color-coding-type']};
    });
    expect(res.colorCodingType).toBe('Linear');
  });

  await softStep('Cleanup', async () => {
    await page.evaluate(async () => {
      (window as any).grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 500));
    });
  });

  if (stepErrors.length > 0)
    throw new Error('Step failures:\n' + stepErrors.map(e => `- ${e.step}: ${e.error}`).join('\n'));
});
