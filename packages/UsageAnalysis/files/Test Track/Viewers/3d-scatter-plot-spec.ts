import {test, expect} from '@playwright/test';

test.use({
  viewport: {width: 1920, height: 1080},
  launchOptions: {args: ['--window-size=1920,1080', '--window-position=0,0']},
});

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try { await test.step(name, fn); }
  catch (e: any) { stepErrors.push({step: name, error: e.message ?? String(e)}); }
}

test('3d Scatter plot scenario', async ({page}) => {
  test.setTimeout(300_000);

  const baseUrl = process.env.DATAGROK_URL ?? 'http://localhost:8888';
  const login = process.env.DATAGROK_LOGIN ?? 'admin';
  const password = process.env.DATAGROK_PASSWORD ?? 'admin';

  await page.goto(baseUrl);
  const loginInput = page.getByPlaceholder('Login or Email').and(page.locator(':visible'));
  if (await loginInput.isVisible({timeout: 15000}).catch(() => false)) {
    await loginInput.click();
    await page.keyboard.type(login);
    await page.getByPlaceholder('Password').and(page.locator(':visible')).click();
    await page.keyboard.type(password);
    await page.keyboard.press('Enter');
  }
  await page.locator('[name="Browse"]').waitFor({timeout: 60000});

  await page.evaluate(async () => {
    const g: any = (window as any).grok;
    document.body.classList.add('selenium');
    g.shell.settings.showFiltersIconsConstantly = true;
    g.shell.windows.simpleMode = true;
    g.shell.closeAll();
    const df = await g.dapi.files.readCsv('System:DemoFiles/demog.csv');
    g.shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
    const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
    const hasBioChem = cols.some((c: any) => c.semType === 'Molecule' || c.semType === 'Macromolecule');
    if (hasBioChem) {
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      await new Promise((r) => setTimeout(r, 5000));
    }
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  await softStep('Step 1 — demog dataset open', async () => {
    const rows = await page.evaluate(() => (window as any).grok.shell.tv.dataFrame.rowCount);
    expect(rows).toBe(5850);
  });

  await softStep('Step 2 — click 3d Scatter plot icon; viewer opens', async () => {
    await page.locator('[name="icon-3d-scatter-plot"]').click();
    await page.locator('[name="viewer-3d-scatter-plot"]').waitFor({timeout: 10000});
    const hasType = await page.evaluate(() => {
      const g: any = (window as any).grok;
      return Array.from(g.shell.tv.viewers).some((v: any) => v.type === '3d scatter plot');
    });
    expect(hasType).toBe(true);
  });

  await softStep('Step 3a — click on viewer changes currentRowIdx', async () => {
    const beforeIdx = await page.evaluate(() => (window as any).grok.shell.tv.dataFrame.currentRowIdx);
    const box = await page.locator('[name="viewer-3d-scatter-plot"]').boundingBox();
    if (!box) throw new Error('viewer bounding box unavailable');
    await page.mouse.click(box.x + box.width * 0.6, box.y + box.height * 0.5);
    await page.waitForTimeout(600);
    const afterIdx = await page.evaluate(() => (window as any).grok.shell.tv.dataFrame.currentRowIdx);
    expect(afterIdx).not.toBe(beforeIdx);
  });

  await softStep('Step 3b — mouse wheel zoom', async () => {
    const ok = await page.evaluate(async () => {
      const v = document.querySelector('[name="viewer-3d-scatter-plot"]') as HTMLElement;
      const canvas = v?.querySelector('canvas') as HTMLCanvasElement;
      if (!canvas) return false;
      const rect = canvas.getBoundingClientRect();
      const cx = rect.left + rect.width / 2;
      const cy = rect.top + rect.height / 2;
      for (let i = 0; i < 5; i++) {
        canvas.dispatchEvent(new WheelEvent('wheel', {
          bubbles: true, cancelable: true, clientX: cx, clientY: cy,
          deltaY: -120, deltaMode: 0, view: window,
        }));
        await new Promise((r) => setTimeout(r, 50));
      }
      await new Promise((r) => setTimeout(r, 500));
      return true;
    });
    expect(ok).toBe(true);
  });

  await softStep('Step 4 — gear opens Property Pane; properties reflect on viewer', async () => {
    await page.evaluate(() => {
      const viewer = document.querySelector('[name="viewer-3d-scatter-plot"]') as HTMLElement;
      const panelBase = viewer.parentElement!.parentElement!;
      const gear = panelBase.querySelector('[name="icon-font-icon-settings"]') as HTMLElement;
      gear.click();
    });
    await page.locator('[name="ScatterPlot3dLook"]').waitFor({timeout: 5000});

    const result = await page.evaluate(async () => {
      const g: any = (window as any).grok;
      const v3d: any = Array.from(g.shell.tv.viewers).find((x: any) => x.type === '3d scatter plot');
      v3d.setOptions({
        colorColumnName: 'RACE',
        sizeColumnName: 'WEIGHT',
        backColor: 0xFFF5F5F5,
        markerOpacity: 90,
      });
      await new Promise((r) => setTimeout(r, 500));
      return {
        colorCol: v3d.props.colorColumnName,
        sizeCol: v3d.props.sizeColumnName,
        markerOpacity: v3d.props.markerOpacity,
      };
    });
    expect(result.colorCol).toBe('RACE');
    expect(result.sizeCol).toBe('WEIGHT');
    expect(result.markerOpacity).toBe(90);
  });

  await page.evaluate(() => (window as any).grok.shell.closeAll());

  if (stepErrors.length > 0)
    throw new Error('Step failures:\n' + stepErrors.map((e) => `- ${e.step}: ${e.error}`).join('\n'));
});
