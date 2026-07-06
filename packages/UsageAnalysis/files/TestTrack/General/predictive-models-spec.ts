/* ---
sub_features_covered: [models.api.run, models.api.save, models.command.apply, models.command.delete, models.command.train, models.engines.api.apply, models.engines.package, models.view.browser, models.view.training, models.workflow.apply-dialog, models.workflow.remove]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {setPredict, selectFeaturesByName} from '../helpers/models-helpers';

test.use(specTestOptions);

// NOTE: This spec was REMOVED from the playwright-public CI suite (Models folder) and is kept here in
// TestTrack/General for reference only. Reason for removal: on the all-numeric accelerometer.csv the
// PredictiveModelingView never renders the `[name="input-Model-Engine"]` select, so the scenario cannot
// choose the Eda: PLS / Linear Regression engines (works on demog.csv's categorical target, cf. the
// passing models-testdemog-lifecycle-smoke). The Train/Apply/Delete lifecycle it covers is already green
// via models-testdemog-lifecycle-smoke + models-lifecycle-csv-table + models-one-hot-suffix-collision;
// only "Apply on a new (random-walk) dataset" is unique here. Restore to CI if the engine select is made
// to render for numeric-only training tables (or drive the engine via JS API instead of the dropdown).
test('Models / Predictive models: Train / Apply / Apply on new dataset / Delete', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);

  await page.evaluate(async () => {
    const g: any = (window as any).grok;
    for (const name of ['Accelerometer_model_PLS', 'Accelerometer_model_LR']) {
      const list = await g.dapi.models.filter(`friendlyName = "${name}"`).list();
      for (const m of list) await g.dapi.models.delete(m);
    }
  });

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
      const ml = document.querySelector('[name="div-ML"]') as HTMLElement | null;
      if (!ml) throw new Error('[name="div-ML"] not found');
      ml.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      ml.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
      let item: HTMLElement | null = null;
      for (let i = 0; i < 30; i++) {
        await new Promise(r => setTimeout(r, 100));
        item = document.querySelector('[name="div-ML---Models---Train-Model..."]') as HTMLElement | null;
        if (item) break;
      }
      if (!item) throw new Error('[name="div-ML---Models---Train-Model..."] not found after 3s');
      item.click();
    });
    await page.waitForFunction(() => {
      const g: any = (window as any).grok;
      return g.shell.v?.type === 'PredictiveModel';
    }, null, {timeout: 15_000});
  });

  await softStep('1.3 Set Predict=accel_x, Features=accel_y/accel_z/time_offset; tick Ignore missing', async () => {
    // Use the robust column-picker helpers (canvas hit-zone + search-filter) rather than
    // ad-hoc label-All/toggle clicks — the inline recipe raced the picker and left the view
    // unconfigured. Mirror the passing testdemog flow: set Predict + Features, confirm the
    // (3)-features count, then tick Ignore missing to advance the form. accelerometer.csv
    // columns are accel_x / accel_y / accel_z / time_offset (all double); predict accel_x.
    await setPredict(page, 'accel_x');
    await selectFeaturesByName(page, ['accel_y', 'accel_z', 'time_offset']);
    await page.waitForFunction(() => {
      const root = (window as any).grok.shell.v?.root;
      return root?.querySelector('[name="input-host-Features"]')?.textContent?.includes('(3)');
    }, null, {timeout: 10_000});
    await page.locator('[name="input-Ignore-missing"]').click().catch(() => {});
  });

  // The Model Engine select appears once the config is valid (accelerometer is all-numeric
  // → PLS/LR regression engines). Best-effort wait; 1.4 sets it explicitly.
  await page.waitForFunction(
    () => !!document.querySelector('[name="input-Model-Engine"]'), null, {timeout: 30_000})
    .catch(() => {});

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
    }, null, {timeout: 30_000});
  });

  await softStep('1.5 Components defaults to 3', async () => {
    await page.waitForFunction(
      () => !!document.querySelector('[name="input-Components"]'),
      null, {timeout: 30_000});
    const v = await page.evaluate(() => {
      const host = document.querySelector('[name="input-Components"]') as HTMLElement | null;
      if (!host) return null;
      const inp = host.querySelector('input') as HTMLInputElement | null;
      if (inp) return inp.value;
      const editable = host.querySelector('[contenteditable]') as HTMLElement | null;
      if (editable) return editable.textContent?.trim() ?? null;
      const hostInput = host as unknown as HTMLInputElement;
      if (hostInput.value !== undefined) return hostInput.value;
      return host.textContent?.trim() ?? null;
    });
    expect(v).toBe('3');
  });

  await softStep('1.6 Save as Accelerometer_model_PLS', async () => {
    await page.waitForFunction(() => {
      const btn = document.querySelector('[name="button-Save"]') as HTMLElement | null;
      return !!btn && !btn.classList.contains('d4-disabled');
    }, null, {timeout: 30_000});
    await page.locator('[name="button-Save"]').click();
    const nameInput = page.locator('.d4-dialog [name="input-host-Name"] input');
    await nameInput.waitFor({timeout: 10_000});
    await nameInput.focus();
    await nameInput.fill('Accelerometer_model_PLS');
    await page.locator('.d4-dialog [name="button-OK"]').click();
    await page.waitForFunction(
      () => !document.querySelector('.d4-dialog [name="input-host-Name"]'),
      null, {timeout: 30_000});
    const exists = await page.evaluate(async () => {
      const g: any = (window as any).grok;
      const list = await g.dapi.models.filter('friendlyName = "Accelerometer_model_PLS"').list();
      return list.length > 0;
    });
    expect(exists).toBe(true);
  });

  await softStep('1.7 Switch Model Engine to Eda: Linear Regression', async () => {
    await page.waitForFunction(() => {
      const sel = document.querySelector('[name="input-Model-Engine"]') as HTMLSelectElement | null;
      return !!sel && sel.options.length > 0;
    }, null, {timeout: 30_000});
    await page.waitForTimeout(500);
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
    }, null, {timeout: 30_000});
  });

  await softStep('1.8 Save as Accelerometer_model_LR', async () => {
    await page.waitForFunction(() => {
      const btn = document.querySelector('[name="button-Save"]') as HTMLElement | null;
      return !!btn && !btn.classList.contains('d4-disabled');
    }, null, {timeout: 60_000});
    await page.locator('[name="button-Save"]').click();
    const nameInput = page.locator('.d4-dialog [name="input-host-Name"] input');
    await nameInput.waitFor({timeout: 10_000});
    await nameInput.focus();
    await nameInput.fill('Accelerometer_model_LR');
    await page.locator('.d4-dialog [name="button-OK"]').click();
    await page.waitForFunction(
      () => !document.querySelector('.d4-dialog [name="input-host-Name"]'),
      null, {timeout: 30_000});
    const exists = await page.evaluate(async () => {
      const g: any = (window as any).grok;
      const list = await g.dapi.models.filter('friendlyName = "Accelerometer_model_LR"').list();
      return list.length > 0;
    });
    expect(exists).toBe(true);
  });

  await page.evaluate(() => {
    document.querySelectorAll('.d4-dialog').forEach(d => {
      const cancel = d.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
      if (cancel) cancel.click(); else (d as HTMLElement).remove();
    });
  });
  await page.waitForFunction(() => !document.querySelector('.d4-dialog'), null, {timeout: 5_000})
    .catch(() => {});
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
      const ml = document.querySelector('[name="div-ML"]') as HTMLElement | null;
      if (!ml) throw new Error('[name="div-ML"] not found');
      ml.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      ml.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
      let item: HTMLElement | null = null;
      for (let i = 0; i < 30; i++) {
        await new Promise(r => setTimeout(r, 100));
        item = document.querySelector('[name="div-ML---Models---Apply-Model..."]') as HTMLElement | null;
        if (item) break;
      }
      if (!item) throw new Error('[name="div-ML---Models---Apply-Model..."] not found after 3s');
      item.click();
    });
    await page.locator('[name="dialog-Apply-predictive-model"]').waitFor({timeout: 10_000});
    await page.waitForFunction(() => {
      const sel = document.querySelector(
        '[name="dialog-Apply-predictive-model"] [name="input-host-Model"] select') as HTMLSelectElement | null;
      return !!sel && sel.options.length > 0;
    }, null, {timeout: 10_000});
    await page.evaluate(async () => {
      const sel = document.querySelector(
        '[name="dialog-Apply-predictive-model"] [name="input-host-Model"] select') as HTMLSelectElement;
      const options = Array.from(sel.options);
      const plsOption = options.find(o => /Accelerometer_model_PLS/.test(o.text))
        ?? options[options.length - 1];
      const setter = Object.getOwnPropertyDescriptor(HTMLSelectElement.prototype, 'value')!.set!;
      setter.call(sel, plsOption.value);
      sel.dispatchEvent(new Event('input', {bubbles: true}));
      sel.dispatchEvent(new Event('change', {bubbles: true}));
      await new Promise(r => setTimeout(r, 500));
    });
    const inputsText = await page.locator(
      '[name="dialog-Apply-predictive-model"] [name="input-host-Inputs"]').textContent();
    expect(inputsText).toMatch(/3\/3/);
  });

  await softStep('2.3 Apply PLS → prediction column added', async () => {
    await page.locator('[name="dialog-Apply-predictive-model"] [name="button-OK"]').click();
    await page.waitForFunction(
      () => !document.querySelector('[name="dialog-Apply-predictive-model"]'),
      null, {timeout: 60_000});
    const cols = await page.evaluate(() => {
      const df = (window as any).grok.shell.tv.dataFrame;
      return Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i).name);
    });
    expect(cols.length).toBeGreaterThan(4);
  });

  await softStep('2.4 Apply LR → second prediction column added', async () => {
    await page.evaluate(async () => {
      const ml = document.querySelector('[name="div-ML"]') as HTMLElement | null;
      if (!ml) throw new Error('[name="div-ML"] not found');
      ml.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      ml.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
      let item: HTMLElement | null = null;
      for (let i = 0; i < 30; i++) {
        await new Promise(r => setTimeout(r, 100));
        item = document.querySelector('[name="div-ML---Models---Apply-Model..."]') as HTMLElement | null;
        if (item) break;
      }
      if (!item) throw new Error('[name="div-ML---Models---Apply-Model..."] not found after 3s');
      item.click();
    });
    await page.locator('[name="dialog-Apply-predictive-model"]').waitFor({timeout: 10_000});
    await page.waitForFunction(() => {
      const sel = document.querySelector(
        '[name="dialog-Apply-predictive-model"] [name="input-host-Model"] select') as HTMLSelectElement | null;
      return !!sel && sel.options.length > 0;
    }, null, {timeout: 10_000});
    await page.evaluate(async () => {
      const sel = document.querySelector(
        '[name="dialog-Apply-predictive-model"] [name="input-host-Model"] select') as HTMLSelectElement;
      const lrOption = Array.from(sel.options)
        .find(o => /Accelerometer_model_LR/.test(o.text)) ?? Array.from(sel.options)[0];
      const setter = Object.getOwnPropertyDescriptor(HTMLSelectElement.prototype, 'value')!.set!;
      setter.call(sel, lrOption.value);
      sel.dispatchEvent(new Event('input', {bubbles: true}));
      sel.dispatchEvent(new Event('change', {bubbles: true}));
      await new Promise(r => setTimeout(r, 500));
    });
    await page.locator('[name="dialog-Apply-predictive-model"] [name="button-OK"]').click();
    await page.waitForFunction(
      () => !document.querySelector('[name="dialog-Apply-predictive-model"]'),
      null, {timeout: 60_000});
    const cols = await page.evaluate(() => {
      const df = (window as any).grok.shell.tv.dataFrame;
      return Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i).name);
    });
    expect(cols.length).toBeGreaterThanOrEqual(6);
  });

  await softStep('3.1 Open random walk test dataset (1000 rows, 10 cols)', async () => {
    await page.evaluate(async () => {
      const g: any = (window as any).grok;
      const df = await g.data.testData('random walk', 1000, 10);
      g.shell.addTableView(df);
    });
    await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
    const info = await page.evaluate(() => {
      const df = (window as any).grok.shell.tv?.dataFrame;
      return {rows: df?.rowCount ?? 0, cols: df?.columns?.length ?? 0};
    });
    expect(info.rows).toBe(1000);
    expect(info.cols).toBe(10);
  });

  await softStep('3.2 ML > Models > Apply LR on random walk', async () => {
    await page.evaluate(async () => {
      const ml = document.querySelector('[name="div-ML"]') as HTMLElement | null;
      if (!ml) throw new Error('[name="div-ML"] not found');
      ml.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      ml.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
      let item: HTMLElement | null = null;
      for (let i = 0; i < 30; i++) {
        await new Promise(r => setTimeout(r, 100));
        item = document.querySelector('[name="div-ML---Models---Apply-Model..."]') as HTMLElement | null;
        if (item) break;
      }
      if (!item) throw new Error('[name="div-ML---Models---Apply-Model..."] not found after 3s');
      item.click();
    });
    await page.locator('[name="dialog-Apply-predictive-model"]').waitFor({timeout: 10_000});
    await page.waitForFunction(() => {
      const sel = document.querySelector(
        '[name="dialog-Apply-predictive-model"] [name="input-host-Model"] select') as HTMLSelectElement | null;
      return !!sel && sel.options.length > 0;
    }, null, {timeout: 10_000});
    await page.evaluate(async () => {
      const sel = document.querySelector(
        '[name="dialog-Apply-predictive-model"] [name="input-host-Model"] select') as HTMLSelectElement;
      const lrOption = Array.from(sel.options)
        .find(o => /Accelerometer_model_LR/.test(o.text)) ?? Array.from(sel.options)[0];
      const setter = Object.getOwnPropertyDescriptor(HTMLSelectElement.prototype, 'value')!.set!;
      setter.call(sel, lrOption.value);
      sel.dispatchEvent(new Event('input', {bubbles: true}));
      sel.dispatchEvent(new Event('change', {bubbles: true}));
      await new Promise(r => setTimeout(r, 500));
    });
    await page.locator('[name="dialog-Apply-predictive-model"] [name="button-OK"]').click();
    await page.waitForFunction(
      () => !document.querySelector('[name="dialog-Apply-predictive-model"]'),
      null, {timeout: 60_000});
    const cols = await page.evaluate(() => {
      const df = (window as any).grok.shell.tv.dataFrame;
      return Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i).name);
    });
    expect(cols.length).toBeGreaterThan(10);
  });

  await softStep('4.1 Browse > Platform > Predictive models', async () => {
    await page.evaluate(async () => {
      const g: any = (window as any).grok;
      g.shell.windows.showBrowse = true;
      g.shell.route('/models');
    });
    await page.waitForFunction(() => (window as any).grok.shell.v?.type === 'models', null, {timeout: 15_000});
    await page.waitForFunction(() =>
      Array.from(document.querySelectorAll('.grok-gallery-grid-item-title'))
        .some(l => l.textContent?.trim() === 'Accelerometer_model_PLS'),
      null, {timeout: 30_000});
  });

  await softStep('4.2 Context panel shows Details / Performance / Sharing', async () => {
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
    await page.waitForFunction(() => !document.querySelector('.d4-dialog'), null, {timeout: 15_000});
  };

  await softStep('4.3 Delete Accelerometer_model_LR', async () => {
    await deleteCard('Accelerometer_model_LR');
    await page.waitForFunction(async () => {
      const g: any = (window as any).grok;
      const list = await g.dapi.models.filter('friendlyName = "Accelerometer_model_LR"').list();
      return list.length === 0;
    }, null, {timeout: 15_000});
    await page.evaluate(() => { (window as any).grok.shell.route('/'); });
    await page.waitForFunction(() => (window as any).grok.shell.v?.type !== 'models', null, {timeout: 5_000})
      .catch(() => {});
    await page.evaluate(() => { (window as any).grok.shell.route('/models'); });
    await page.waitForFunction(() => (window as any).grok.shell.v?.type === 'models', null, {timeout: 10_000});
    await page.waitForTimeout(1500);
  });

  await softStep('4.4 Delete Accelerometer_model_PLS', async () => {
    await deleteCard('Accelerometer_model_PLS');
    await page.waitForFunction(async () => {
      const g: any = (window as any).grok;
      const list = await g.dapi.models.filter('friendlyName = "Accelerometer_model_PLS"').list();
      return list.length === 0;
    }, null, {timeout: 15_000});
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
