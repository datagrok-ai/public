import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Predictive models: Train / Apply / Apply on new dataset / Delete', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);

  // ── Scenario 1: Train ────────────────────────────────────────────────────
  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    const g: any = (window as any).grok;
    g.shell.settings.showFiltersIconsConstantly = true;
    g.shell.windows.simpleMode = true;
    g.shell.closeAll();
    const df = await g.dapi.files.readCsv('System:DemoFiles/sensors/accelerometer.csv');
    g.shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});

  await softStep('1.1 Open Demo/Sensors/accelerometer.csv', async () => {
    const info = await page.evaluate(() => {
      const df = (window as any).grok.shell.tv?.dataFrame;
      return {rows: df?.rowCount ?? 0, cols: df?.columns?.length ?? 0};
    });
    expect(info.rows).toBeGreaterThan(0);
    expect(info.cols).toBe(4);
  });

  await softStep('1.2 Open ML > Models > Train Model...', async () => {
    await page.evaluate(async () => {
      const ml = document.querySelector('[name="div-ML"]') as HTMLElement;
      ml.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      ml.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
      await new Promise(r => setTimeout(r, 400));
      (document.querySelector('[name="div-ML---Models---Train-Model..."]') as HTMLElement).click();
    });
    await page.waitForFunction(() => {
      const g: any = (window as any).grok;
      return g.shell.v?.type === 'PredictiveModel';
    }, {timeout: 15_000});
  });

  await softStep('1.3 Set Features to accel_y, accel_z, time_offset', async () => {
    await page.locator('[name="div-Features"]').click();
    await page.locator('.d4-dialog[name="dialog-Select-columns..."]').waitFor({timeout: 10_000});
    await page.locator('.d4-dialog [name="label-All"]').click();
    const counter = await page.evaluate(async () => {
      const overlay = document.querySelector('.d4-dialog [name="viewer-Grid"] [name="overlay"]') as HTMLElement;
      const canvas = document.querySelector('.d4-dialog [name="viewer-Grid"] [name="canvas"]') as HTMLElement;
      const rect = canvas.getBoundingClientRect();
      const clientX = rect.right - 20;
      const clientY = rect.top + 34;
      ['mousedown', 'mouseup', 'click'].forEach(type => {
        overlay.dispatchEvent(new MouseEvent(type, {bubbles: true, cancelable: true, clientX, clientY, button: 0}));
      });
      await new Promise(r => setTimeout(r, 300));
      return document.querySelector('.d4-dialog label[style*="margin-left"]')?.textContent;
    });
    expect(counter).toMatch(/3\s*checked/);
    await page.locator('.d4-dialog [name="button-OK"]').click();
    await page.waitForFunction(() =>
      document.querySelector('[name="input-host-Features"] .ui-input-column-names')?.textContent?.includes('(3)'),
      {timeout: 10_000});
  });

  await softStep('1.4 Set Model Engine to Eda: PLS Regression', async () => {
    await page.evaluate(async () => {
      const sel = document.querySelector('[name="input-Model-Engine"]') as HTMLSelectElement;
      sel.focus();
      const setter = Object.getOwnPropertyDescriptor(HTMLSelectElement.prototype, 'value')!.set!;
      setter.call(sel, 'Eda: PLS Regression');
      sel.dispatchEvent(new Event('input', {bubbles: true}));
      sel.dispatchEvent(new Event('change', {bubbles: true}));
    });
    await page.waitForFunction(() => {
      const text = Array.from(document.querySelectorAll('h3, h4, [class*="card-header"]'))
        .map(e => e.textContent?.trim() ?? '').find(t => /Eda:\s*PLS/i.test(t));
      return !!text;
    }, {timeout: 30_000});
  });

  await softStep('1.5 Components defaults to 3', async () => {
    const v = await page.locator('[name="input-Components"]').inputValue();
    expect(v).toBe('3');
  });

  await softStep('1.6 Save as Accelerometer_model_PLS', async () => {
    await page.locator('[name="button-Save"]').click();
    await page.locator('.d4-dialog [name="input-Name"]').waitFor({timeout: 10_000});
    await page.evaluate(() => {
      const input = document.querySelector('.d4-dialog [name="input-Name"]') as HTMLInputElement;
      input.focus();
      const setter = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value')!.set!;
      setter.call(input, 'Accelerometer_model_PLS');
      input.dispatchEvent(new Event('input', {bubbles: true}));
      input.dispatchEvent(new Event('change', {bubbles: true}));
    });
    await page.locator('.d4-dialog [name="button-OK"]').click();
    await page.waitForFunction(() => !document.querySelector('.d4-dialog'), {timeout: 30_000});
    const exists = await page.evaluate(async () => {
      const g: any = (window as any).grok;
      const list = await g.dapi.models.filter('friendlyName = "Accelerometer_model_PLS"').list();
      return list.length > 0;
    });
    expect(exists).toBe(true);
  });

  await softStep('1.7 Switch Model Engine to Eda: Linear Regression', async () => {
    await page.evaluate(async () => {
      const sel = document.querySelector('[name="input-Model-Engine"]') as HTMLSelectElement;
      sel.focus();
      const setter = Object.getOwnPropertyDescriptor(HTMLSelectElement.prototype, 'value')!.set!;
      setter.call(sel, 'Eda: Linear Regression');
      sel.dispatchEvent(new Event('input', {bubbles: true}));
      sel.dispatchEvent(new Event('change', {bubbles: true}));
    });
    await page.waitForFunction(() => {
      const text = Array.from(document.querySelectorAll('h3, h4, [class*="card-header"]'))
        .map(e => e.textContent?.trim() ?? '').find(t => /Eda:\s*Linear/i.test(t) && !/PLS/i.test(t));
      return !!text;
    }, {timeout: 30_000});
  });

  await softStep('1.8 Save as Accelerometer_model_LR', async () => {
    await page.locator('[name="button-Save"]').click();
    await page.locator('.d4-dialog [name="input-Name"]').waitFor({timeout: 10_000});
    await page.evaluate(() => {
      const input = document.querySelector('.d4-dialog [name="input-Name"]') as HTMLInputElement;
      input.focus();
      const setter = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value')!.set!;
      setter.call(input, 'Accelerometer_model_LR');
      input.dispatchEvent(new Event('input', {bubbles: true}));
      input.dispatchEvent(new Event('change', {bubbles: true}));
    });
    await page.locator('.d4-dialog [name="button-OK"]').click();
    await page.waitForFunction(() => !document.querySelector('.d4-dialog'), {timeout: 30_000});
    const exists = await page.evaluate(async () => {
      const g: any = (window as any).grok;
      const list = await g.dapi.models.filter('friendlyName = "Accelerometer_model_LR"').list();
      return list.length > 0;
    });
    expect(exists).toBe(true);
  });

  // ── Scenario 2: Apply ────────────────────────────────────────────────────
  await page.evaluate(async () => {
    const g: any = (window as any).grok;
    g.shell.closeAll();
    await new Promise(r => setTimeout(r, 500));
    const df = await g.dapi.files.readCsv('System:DemoFiles/sensors/accelerometer.csv');
    g.shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});

  await softStep('2.1 Re-open accelerometer.csv', async () => {
    const info = await page.evaluate(() => {
      const df = (window as any).grok.shell.tv?.dataFrame;
      return {rows: df?.rowCount ?? 0};
    });
    expect(info.rows).toBeGreaterThan(0);
  });

  await softStep('2.2 ML > Models > Apply Model... → PLS model, inputs (3/3)', async () => {
    await page.evaluate(async () => {
      const ml = document.querySelector('[name="div-ML"]') as HTMLElement;
      ml.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      ml.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
      await new Promise(r => setTimeout(r, 400));
      (document.querySelector('[name="div-ML---Models---Apply-Model..."]') as HTMLElement).click();
    });
    await page.locator('.d4-dialog').waitFor({timeout: 10_000});
    await page.evaluate(async () => {
      const sel = document.querySelector('[name="input-host-Model"] select') as HTMLSelectElement;
      const options = Array.from(sel.options);
      // PLS was saved first → has older timestamp → shown at the bottom of dropdown
      const plsOption = options[options.length - 1];
      const setter = Object.getOwnPropertyDescriptor(HTMLSelectElement.prototype, 'value')!.set!;
      setter.call(sel, plsOption.value);
      sel.dispatchEvent(new Event('input', {bubbles: true}));
      sel.dispatchEvent(new Event('change', {bubbles: true}));
      await new Promise(r => setTimeout(r, 500));
    });
    const inputsText = await page.locator('[name="input-host-Inputs"]').textContent();
    expect(inputsText).toMatch(/3\/3/);
  });

  await softStep('2.3 Apply PLS → prediction column added', async () => {
    await page.locator('.d4-dialog [name="button-OK"]').click();
    await page.waitForFunction(() => !document.querySelector('.d4-dialog'), {timeout: 60_000});
    const cols = await page.evaluate(() => {
      const df = (window as any).grok.shell.tv.dataFrame;
      return Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i).name);
    });
    expect(cols.length).toBeGreaterThan(4);
  });

  await softStep('2.4 Apply LR → second prediction column added', async () => {
    await page.evaluate(async () => {
      const ml = document.querySelector('[name="div-ML"]') as HTMLElement;
      ml.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      ml.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
      await new Promise(r => setTimeout(r, 400));
      (document.querySelector('[name="div-ML---Models---Apply-Model..."]') as HTMLElement).click();
    });
    await page.locator('.d4-dialog').waitFor({timeout: 10_000});
    await page.evaluate(async () => {
      const sel = document.querySelector('[name="input-host-Model"] select') as HTMLSelectElement;
      // LR saved second → newer timestamp → first in dropdown
      const lrOption = Array.from(sel.options)[0];
      const setter = Object.getOwnPropertyDescriptor(HTMLSelectElement.prototype, 'value')!.set!;
      setter.call(sel, lrOption.value);
      sel.dispatchEvent(new Event('input', {bubbles: true}));
      sel.dispatchEvent(new Event('change', {bubbles: true}));
      await new Promise(r => setTimeout(r, 500));
    });
    await page.locator('.d4-dialog [name="button-OK"]').click();
    await page.waitForFunction(() => !document.querySelector('.d4-dialog'), {timeout: 60_000});
    const cols = await page.evaluate(() => {
      const df = (window as any).grok.shell.tv.dataFrame;
      return Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i).name);
    });
    expect(cols.length).toBeGreaterThanOrEqual(6);
  });

  // ── Scenario 3: Apply on a new dataset ───────────────────────────────────
  await softStep('3.1 Tools > Dev > Open Test Dataset (1000 rows, 10 cols, random walk)', async () => {
    await page.evaluate(() => {
      const g: any = (window as any).grok;
      const toolsMenu = g.shell.topMenu.find('Tools');
      const item = toolsMenu.root.querySelector('[name="div-Tools---Dev---Open-Test-Dataset"]');
      item.click();
    });
    await page.locator('.d4-dialog [name="input-rows"]').waitFor({timeout: 10_000});
    await page.evaluate(() => {
      const setInput = (sel: string, v: string) => {
        const el = document.querySelector(sel) as HTMLInputElement;
        el.focus();
        const setter = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value')!.set!;
        setter.call(el, v);
        el.dispatchEvent(new Event('input', {bubbles: true}));
        el.dispatchEvent(new Event('change', {bubbles: true}));
      };
      setInput('[name="input-rows"]', '1000');
      setInput('[name="input-cols"]', '10');
      const radios = Array.from(document.querySelectorAll('.d4-dialog [type="radio"]')) as HTMLInputElement[];
      // Radio order: demog(0), biosensor(1), plates(2), random walk(3), ...
      radios[3].click();
    });
    await page.locator('.d4-dialog [name="button-OK"]').click();
    await page.waitForFunction(() => (window as any).grok.shell.v?.name?.includes('randomWalk'), {timeout: 15_000});
  });

  await softStep('3.2 ML > Models > Apply LR on random walk', async () => {
    await page.evaluate(async () => {
      const ml = document.querySelector('[name="div-ML"]') as HTMLElement;
      ml.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      ml.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
      await new Promise(r => setTimeout(r, 400));
      (document.querySelector('[name="div-ML---Models---Apply-Model..."]') as HTMLElement).click();
    });
    await page.locator('.d4-dialog').waitFor({timeout: 10_000});
    await page.evaluate(async () => {
      const sel = document.querySelector('[name="input-host-Model"] select') as HTMLSelectElement;
      const lrOption = Array.from(sel.options)[0];
      const setter = Object.getOwnPropertyDescriptor(HTMLSelectElement.prototype, 'value')!.set!;
      setter.call(sel, lrOption.value);
      sel.dispatchEvent(new Event('input', {bubbles: true}));
      sel.dispatchEvent(new Event('change', {bubbles: true}));
      await new Promise(r => setTimeout(r, 500));
    });
    await page.locator('.d4-dialog [name="button-OK"]').click();
    await page.waitForFunction(() => !document.querySelector('.d4-dialog'), {timeout: 60_000});
    const cols = await page.evaluate(() => {
      const df = (window as any).grok.shell.tv.dataFrame;
      return Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i).name);
    });
    expect(cols).toContain('accel_x');
  });

  // ── Scenario 4: Delete ───────────────────────────────────────────────────
  await softStep('4.1 Browse > Platform > Predictive models', async () => {
    await page.evaluate(async () => {
      const g: any = (window as any).grok;
      g.shell.windows.showBrowse = true;
      await new Promise(r => setTimeout(r, 600));
      const labels = Array.from(document.querySelectorAll('.d4-tree-view-group-label')) as HTMLElement[];
      const pm = labels.find(l => l.textContent?.trim() === 'Predictive models')!;
      pm.click();
    });
    await page.waitForFunction(() =>
      Array.from(document.querySelectorAll('.grok-gallery-grid-item-title'))
        .some(l => l.textContent?.trim() === 'Accelerometer_model_PLS'),
      {timeout: 15_000});
  });

  await softStep('4.2 Context panel shows Details / Performance / Activity / Sharing', async () => {
    await page.evaluate(() => {
      const label = Array.from(document.querySelectorAll('.grok-gallery-grid-item-title'))
        .find(l => l.textContent?.trim() === 'Accelerometer_model_PLS') as HTMLElement;
      (label.closest('.grok-gallery-grid-item') as HTMLElement).click();
    });
    await page.waitForTimeout(1000);
    const panes = await page.evaluate(() =>
      Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
        .map(h => h.textContent?.trim() ?? ''));
    expect(panes).toContain('Details');
    expect(panes).toContain('Performance');
    expect(panes.some((p: string) => p.startsWith('Activity'))).toBe(true);
    expect(panes).toContain('Sharing');
  });

  const deleteCard = async (name: string) => {
    await page.evaluate((targetName) => {
      const label = Array.from(document.querySelectorAll('.grok-gallery-grid-item-title'))
        .find(l => l.textContent?.trim() === targetName) as HTMLElement;
      const card = label.closest('.grok-gallery-grid-item') as HTMLElement;
      card.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true, button: 2}));
    }, name);
    await page.locator('.d4-menu-popup').waitFor({timeout: 5000});
    await page.evaluate(() => {
      const popup = document.querySelector('.d4-menu-popup')!;
      const del = Array.from(popup.querySelectorAll('.d4-menu-item'))
        .find(i => i.textContent?.trim() === 'Delete') as HTMLElement;
      del.click();
    });
    await page.locator('.d4-dialog').waitFor({timeout: 5000});
    await page.evaluate(() => {
      const dlg = document.querySelector('.d4-dialog')!;
      const btn = Array.from(dlg.querySelectorAll('button'))
        .find(b => b.textContent?.trim() === 'DELETE') as HTMLElement;
      btn.click();
    });
    await page.waitForFunction(() => !document.querySelector('.d4-dialog'), {timeout: 15_000});
  };

  await softStep('4.3 Delete Accelerometer_model_LR', async () => {
    await deleteCard('Accelerometer_model_LR');
    await page.waitForFunction(async () => {
      const g: any = (window as any).grok;
      const list = await g.dapi.models.filter('friendlyName = "Accelerometer_model_LR"').list();
      return list.length === 0;
    }, {timeout: 15_000});
  });

  await softStep('4.4 Delete Accelerometer_model_PLS', async () => {
    await deleteCard('Accelerometer_model_PLS');
    await page.waitForFunction(async () => {
      const g: any = (window as any).grok;
      const list = await g.dapi.models.filter('friendlyName = "Accelerometer_model_PLS"').list();
      return list.length === 0;
    }, {timeout: 15_000});
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
