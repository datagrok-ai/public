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

test('3D Scatter Plot — Test Track scenarios', async ({page}) => {
  test.setTimeout(600_000);

  await page.goto(baseUrl);
  const submitLogin = async () => {
    const loginInput = page.getByPlaceholder('Login or Email').and(page.locator(':visible'));
    if (!(await loginInput.isVisible({timeout: 15000}).catch(() => false))) return;
    await loginInput.click();
    await page.keyboard.type(login);
    await page.getByPlaceholder('Password').and(page.locator(':visible')).click();
    await page.keyboard.type(password);
    await page.keyboard.press('Enter');
  };
  await submitLogin();
  const browse = page.locator('[name="Browse"]');
  try {
    await browse.waitFor({timeout: 30000});
  } catch {
    // Login form may have silently rejected the first submission; retry once.
    await submitLogin();
    await browse.waitFor({timeout: 120000});
  }

  // Phase 2: open dataset
  await page.evaluate(async () => {
    const g: any = (window as any).grok;
    document.body.classList.add('selenium');
    g.shell.settings.showFiltersIconsConstantly = true;
    g.shell.windows.simpleMode = true;
    g.shell.closeAll();
    const df = await g.dapi.files.readCsv('System:DemoFiles/demog.csv');
    g.shell.addTableView(df);
    await new Promise<void>(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').first().waitFor({timeout: 30000});

  // Phase 3: add 3D Scatter Plot
  await page.evaluate(() => {
    const icon = document.querySelector('[name="icon-3d-scatter-plot"]') as HTMLElement;
    icon.click();
  });
  await page.locator('[name="viewer-3d-scatter-plot"]').first().waitFor({timeout: 30000});
  await page.waitForTimeout(500);

  const sp = () => page.evaluate(() => {
    const g: any = (window as any).grok;
    return !!g.shell.tv.viewers.find((v: any) => v.type === '3d scatter plot');
  });
  expect(await sp()).toBe(true);

  // JS API helpers — scatter-plot-spec.ts pattern
  const setProps = async (patch: Record<string, any>) => {
    await page.evaluate((p) => {
      const g: any = (window as any).grok;
      const v: any = g.shell.tv.viewers.find((x: any) => x.type === '3d scatter plot');
      for (const k of Object.keys(p)) v.props[k] = p[k];
    }, patch);
    await page.waitForTimeout(100);
  };
  const getProps = async (keys: string[]) => {
    return await page.evaluate((ks) => {
      const g: any = (window as any).grok;
      const v: any = g.shell.tv.viewers.find((x: any) => x.type === '3d scatter plot');
      const o: any = {};
      for (const k of ks) o[k] = v.props[k];
      return o;
    }, keys);
  };

  await softStep('Axis column assignment', async () => {
    await setProps({xColumnName: 'AGE', yColumnName: 'HEIGHT', zColumnName: 'WEIGHT'});
    let s = await getProps(['xColumnName', 'yColumnName', 'zColumnName']);
    expect(s).toEqual({xColumnName: 'AGE', yColumnName: 'HEIGHT', zColumnName: 'WEIGHT'});

    await setProps({xColumnName: 'WEIGHT', yColumnName: 'AGE', zColumnName: 'HEIGHT'});
    s = await getProps(['xColumnName', 'yColumnName', 'zColumnName']);
    expect(s).toEqual({xColumnName: 'WEIGHT', yColumnName: 'AGE', zColumnName: 'HEIGHT'});

    await setProps({xColumnName: 'AGE', yColumnName: 'HEIGHT', zColumnName: 'WEIGHT'});
    s = await getProps(['xColumnName', 'yColumnName', 'zColumnName']);
    expect(s).toEqual({xColumnName: 'AGE', yColumnName: 'HEIGHT', zColumnName: 'WEIGHT'});
  });

  await softStep('Axis types', async () => {
    await setProps({xAxisType: 'logarithmic'});
    expect((await getProps(['xAxisType'])).xAxisType).toBe('logarithmic');
    await setProps({yAxisType: 'logarithmic'});
    expect((await getProps(['yAxisType'])).yAxisType).toBe('logarithmic');
    await setProps({zAxisType: 'logarithmic'});
    expect((await getProps(['zAxisType'])).zAxisType).toBe('logarithmic');
    await setProps({xAxisType: 'linear', yAxisType: 'linear', zAxisType: 'linear'});
    const s = await getProps(['xAxisType', 'yAxisType', 'zAxisType']);
    expect(s).toEqual({xAxisType: 'linear', yAxisType: 'linear', zAxisType: 'linear'});
  });

  await softStep('Color coding - categorical', async () => {
    await setProps({colorColumnName: 'SEX'});
    expect((await getProps(['colorColumnName'])).colorColumnName).toBe('SEX');
    await setProps({colorColumnName: 'RACE'});
    expect((await getProps(['colorColumnName'])).colorColumnName).toBe('RACE');
    await setProps({colorColumnName: ''});
    expect((await getProps(['colorColumnName'])).colorColumnName).toBeFalsy();
  });

  await softStep('Color coding - numerical', async () => {
    await setProps({colorColumnName: 'AGE'});
    expect((await getProps(['colorColumnName'])).colorColumnName).toBe('AGE');
    await setProps({colorColumnName: ''});
    expect((await getProps(['colorColumnName'])).colorColumnName).toBeFalsy();
  });

  await softStep('Size coding', async () => {
    await setProps({sizeColumnName: 'WEIGHT'});
    expect((await getProps(['sizeColumnName'])).sizeColumnName).toBe('WEIGHT');
    await setProps({sizeColumnName: 'AGE'});
    expect((await getProps(['sizeColumnName'])).sizeColumnName).toBe('AGE');
    await setProps({sizeColumnName: ''});
    expect((await getProps(['sizeColumnName'])).sizeColumnName).toBeFalsy();
  });

  await softStep('Labels', async () => {
    await setProps({labelColumnName: 'SEX'});
    expect((await getProps(['labelColumnName'])).labelColumnName).toBe('SEX');
    await setProps({labelColumnName: ''});
    expect((await getProps(['labelColumnName'])).labelColumnName).toBeFalsy();
  });

  await softStep('Marker type', async () => {
    for (const m of ['sphere', 'box', 'cylinder', 'tetrahedron', 'dodecahedron', 'octahedron']) {
      await setProps({markerType: m});
      expect((await getProps(['markerType'])).markerType).toBe(m);
    }
  });

  await softStep('Marker opacity and rotation', async () => {
    for (const v of [20, 100, 69]) {
      await setProps({markerOpacity: v});
      expect((await getProps(['markerOpacity'])).markerOpacity).toBe(v);
    }
    await setProps({markerRandomRotation: true});
    expect((await getProps(['markerRandomRotation'])).markerRandomRotation).toBe(true);
    await setProps({markerRandomRotation: false});
    expect((await getProps(['markerRandomRotation'])).markerRandomRotation).toBe(false);
  });

  await softStep('Filtered out points', async () => {
    await page.evaluate(async () => {
      const g: any = (window as any).grok;
      g.shell.tv.getFiltersGroup().updateOrAdd({type: 'histogram', column: 'AGE', min: 20, max: 40});
      await new Promise(r => setTimeout(r, 300));
    });
    await setProps({showFilteredOutPoints: true});
    expect((await getProps(['showFilteredOutPoints'])).showFilteredOutPoints).toBe(true);
    await setProps({showFilteredOutPoints: false});
    expect((await getProps(['showFilteredOutPoints'])).showFilteredOutPoints).toBe(false);
    const res = await page.evaluate(async () => {
      const g: any = (window as any).grok;
      const col = g.shell.tv.dataFrame.col('AGE');
      const s = col.stats;
      g.shell.tv.getFiltersGroup().updateOrAdd({type: 'histogram', column: 'AGE', min: s.min, max: s.max});
      await new Promise(r => setTimeout(r, 300));
      const df = g.shell.tv.dataFrame;
      return {f: df.filter.trueCount, t: df.rowCount};
    });
    expect(res.f).toBe(res.t);
  });

  await softStep('Axes visibility and grid lines', async () => {
    await setProps({showAxes: false});
    expect((await getProps(['showAxes'])).showAxes).toBe(false);
    await setProps({showAxes: true});
    expect((await getProps(['showAxes'])).showAxes).toBe(true);
    await setProps({showVerticalGridLines: false});
    expect((await getProps(['showVerticalGridLines'])).showVerticalGridLines).toBe(false);
    await setProps({showHorizontalGridLines: false});
    expect((await getProps(['showHorizontalGridLines'])).showHorizontalGridLines).toBe(false);
    await setProps({showVerticalGridLines: true, showHorizontalGridLines: true});
    const s = await getProps(['showVerticalGridLines', 'showHorizontalGridLines']);
    expect(s).toEqual({showVerticalGridLines: true, showHorizontalGridLines: true});
  });

  await softStep('Background and colors (JS API)', async () => {
    const origAxis = (await getProps(['axisLineColor'])).axisLineColor;
    await setProps({backColor: 0xFF000000});
    expect((await getProps(['backColor'])).backColor).toBe(0xFF000000);
    await setProps({axisLineColor: 0xFFFFFFFF});
    expect((await getProps(['axisLineColor'])).axisLineColor).toBe(0xFFFFFFFF);
    await setProps({backColor: 0xFFFFFFFF, axisLineColor: origAxis});
    const s = await getProps(['backColor', 'axisLineColor']);
    expect(s.backColor).toBe(0xFFFFFFFF);
    expect(s.axisLineColor).toBe(origAxis);
  });

  await softStep('Dynamic camera movement', async () => {
    await setProps({dynamicCameraMovement: true});
    expect((await getProps(['dynamicCameraMovement'])).dynamicCameraMovement).toBe(true);
    await setProps({dynamicCameraMovement: false});
    expect((await getProps(['dynamicCameraMovement'])).dynamicCameraMovement).toBe(false);
  });

  await softStep('Zoom and navigation', async () => {
    await page.evaluate(async () => {
      const viewer = document.querySelector('[name="viewer-3d-scatter-plot"]') as HTMLElement;
      const canvas = viewer.querySelector('canvas') as HTMLCanvasElement;
      const rect = canvas.getBoundingClientRect();
      const cx = rect.left + rect.width / 2;
      const cy = rect.top + rect.height / 2;
      for (let i = 0; i < 5; i++) {
        canvas.dispatchEvent(new WheelEvent('wheel', {bubbles: true, cancelable: true, clientX: cx, clientY: cy, deltaY: -120, deltaMode: 0}));
        await new Promise(r => setTimeout(r, 50));
      }
      for (let i = 0; i < 5; i++) {
        canvas.dispatchEvent(new WheelEvent('wheel', {bubbles: true, cancelable: true, clientX: cx, clientY: cy, deltaY: 120, deltaMode: 0}));
        await new Promise(r => setTimeout(r, 50));
      }
      canvas.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, button: 2, clientX: cx, clientY: cy}));
      await new Promise(r => setTimeout(r, 300));
      const label = Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label'))
        .find(e => (e as HTMLElement).textContent?.trim() === 'Reset View') as HTMLElement;
      if (label) (label.closest('.d4-menu-item') as HTMLElement).click();
      await new Promise(r => setTimeout(r, 300));
    });
    const menuOpen = await page.evaluate(() => !!document.querySelector('.d4-menu-popup'));
    expect(menuOpen).toBe(false);
  });

  await softStep('Mouse-over row group highlight', async () => {
    await page.evaluate(async () => {
      const icon = document.querySelector('[name="icon-bar-chart"]') as HTMLElement;
      icon.click();
      for (let i = 0; i < 30; i++) {
        if (document.querySelector('[name="viewer-Bar-chart"]')) break;
        await new Promise(r => setTimeout(r, 100));
      }
    });
    expect(await page.evaluate(() => !!document.querySelector('[name="viewer-Bar-chart"]'))).toBe(true);
    expect((await getProps(['showMouseOverRowGroup'])).showMouseOverRowGroup).toBe(true);

    await page.evaluate(async () => {
      const bc = document.querySelector('[name="viewer-Bar-chart"]') as HTMLElement;
      const canvas = bc.querySelector('canvas') as HTMLCanvasElement;
      const rect = canvas.getBoundingClientRect();
      canvas.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, clientX: rect.left + rect.width * 0.3, clientY: rect.top + rect.height * 0.6}));
      await new Promise(r => setTimeout(r, 200));
    });

    await setProps({showMouseOverRowGroup: false});
    expect((await getProps(['showMouseOverRowGroup'])).showMouseOverRowGroup).toBe(false);

    await page.evaluate(async () => {
      const bc = document.querySelector('[name="viewer-Bar-chart"]') as HTMLElement;
      const canvas = bc.querySelector('canvas') as HTMLCanvasElement;
      const rect = canvas.getBoundingClientRect();
      canvas.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, clientX: rect.left + rect.width * 0.5, clientY: rect.top + rect.height * 0.6}));
      await new Promise(r => setTimeout(r, 200));
    });

    await setProps({showMouseOverRowGroup: true});
    expect((await getProps(['showMouseOverRowGroup'])).showMouseOverRowGroup).toBe(true);

    await page.evaluate(async () => {
      const bc = document.querySelector('[name="viewer-Bar-chart"]') as HTMLElement;
      const pb = bc.closest('.panel-base') as HTMLElement;
      const closeBtn = pb?.querySelector('.panel-titlebar-button-close') as HTMLElement;
      closeBtn?.click();
      await new Promise(r => setTimeout(r, 300));
    });
    expect(await page.evaluate(() => !!document.querySelector('[name="viewer-Bar-chart"]'))).toBe(false);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
