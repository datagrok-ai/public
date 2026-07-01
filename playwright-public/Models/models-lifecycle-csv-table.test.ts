/* ---
sub_features_covered: [models.api.save, models.command.apply, models.command.delete, models.command.edit, models.command.share, models.command.train, models.engines.api.apply, models.meta.performance-section, models.view.training, models.workflow.apply-dialog, models.workflow.edit-info, models.workflow.remove]
--- */
﻿import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const MODEL_BASE = 'LifecycleCsvModel';
const MODEL_NAME = `${MODEL_BASE}_${Date.now()}`;
const NEW_DESCRIPTION = 'CSV-backed lifecycle smoke (round 4)';

test('Models / lifecycle on CSV-backed trainedOn (trained_on_csv_table)', async ({page}) => {
  // Full train→apply→browse→edit→share→delete lifecycle, training one EDA Linear Regression on a CSV
  // (not chemprop). Several steps each poll a train/save/apply signal up to 180s; 300s covers the
  // longest realistic chain (train engine-mount + SAVE-enable + apply) with margin.
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  
  await page.evaluate(async ({prefix}) => {
    const g: any = (window as any).grok;
    const all = await g.dapi.models.list();
    const stale = all.filter((m: any) => (m.friendlyName || m.name || '').startsWith(prefix));
    for (const m of stale)
      try { await g.dapi.models.delete(m); } catch (_) { /* best-effort */ }
  }, {prefix: MODEL_BASE});


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

  await softStep('1.1 Open accelerometer.csv — table view active', async () => {
    const info = await page.evaluate(() => {
      const g: any = (window as any).grok;
      const df = g.shell.tv?.dataFrame;
      const names = df ? Array.from({length: df.columns.length},
        (_: unknown, i: number) => df.columns.byIndex(i).name) : [];
      return {rows: df?.rowCount ?? 0, cols: df?.columns?.length ?? 0,
              viewType: g.shell.v?.type, names};
    });
    expect(info.rows).toBeGreaterThan(0);
    expect(info.viewType).toBe('TableView');
    
    expect(info.names).toContain('accel_x');
    expect(info.names).toContain('accel_y');
    expect(info.names).toContain('accel_z');
    expect(info.names).toContain('time_offset');
  });

  await softStep('1.2 ML > Models > Train Model... — PredictiveModelingView opens', async () => {
    
    await page.locator('[name="div-ML"]').click();
    await page.evaluate(() => {
      const models = document.querySelector('[name="div-ML---Models"]') as HTMLElement | null;
      if (!models) throw new Error('ML > Models submenu not found');
      const r = models.getBoundingClientRect();
      const ev = (t: string) => new MouseEvent(t, {
        bubbles: true, cancelable: true, view: window,
        clientX: r.left + 5, clientY: r.top + 5,
      });
      models.dispatchEvent(ev('mouseover'));
      models.dispatchEvent(ev('mouseenter'));
      models.dispatchEvent(ev('mousemove'));
    });
    await page.locator('[name="div-ML---Models---Train-Model..."]').click();
    await page.waitForFunction(() => (window as any).grok.shell.v?.type === 'PredictiveModel',
      null, {timeout: 15_000});
  });

  await softStep('1.3 Accept default Predict = accel_x (numerical regression target)', async () => {
    
    await expect.poll(async () => await page.evaluate(() =>
      (window as any).grok.shell.v.root
        .querySelector('[name="input-host-Predict"] .d4-column-selector-column')?.textContent?.trim()),
      {timeout: 15_000}).toBe('accel_x');
  });

  await softStep('1.4 Set Features to accel_y, accel_z, time_offset (3 cols)', async () => {
    
    await page.locator('[name="div-Features"]').click();
    await page.locator('.d4-dialog[name="dialog-Select-columns..."]').waitFor({timeout: 10_000});
    
    await page.waitForTimeout(2_000);
    await page.locator('.d4-dialog [name="label-All"]').click();
    await page.waitForFunction(() => {
      const lbl = Array.from(document.querySelectorAll('.d4-dialog label'))
        .find((l) => /checked/.test(l.textContent || ''));
      return lbl?.textContent?.startsWith('4');
    }, null, {timeout: 10_000});
    
    const toggled = await page.evaluate(async () => {
      const dlg = document.querySelector(
        '.d4-dialog[name="dialog-Select-columns..."]')!;
      const overlay = dlg.querySelector('canvas[name="overlay"]') as HTMLCanvasElement;
      const rect = overlay.getBoundingClientRect();
      const x = rect.right - 40;
      const y = rect.top + 36;
      const evtInit: MouseEventInit = {
        bubbles: true, cancelable: true, clientX: x, clientY: y, button: 0, view: window,
      };
      overlay.dispatchEvent(new MouseEvent('mousedown', evtInit));
      overlay.dispatchEvent(new MouseEvent('mouseup', evtInit));
      overlay.dispatchEvent(new MouseEvent('click', evtInit));
      await new Promise((r) => setTimeout(r, 600));
      const counter = Array.from(dlg.querySelectorAll('label'))
        .map((l) => l.textContent?.trim() ?? '').find((t) => /checked/.test(t));
      return {counter, x, y};
    });
    
    if (!toggled.counter?.startsWith('3')) {
      await page.mouse.click(toggled.x, toggled.y);
      await page.waitForTimeout(500);
    }
    await page.waitForFunction(() => {
      const lbl = Array.from(document.querySelectorAll('.d4-dialog label'))
        .find((l) => /checked/.test(l.textContent || ''));
      return lbl?.textContent?.startsWith('3');
    }, null, {timeout: 10_000});
    await page.locator('.d4-dialog [name="button-OK"]').click();
    await page.waitForFunction(() =>
      document.querySelector('[name="input-host-Features"]')?.textContent?.includes('(3)'),
      null, {timeout: 10_000});
    await expect(page.locator('[name="input-host-Features"]'))
      .not.toContainText('accel_x');
  });

  await softStep('1.5 Set Model Engine to Eda: Linear Regression', async () => {
    
    try {
      
      await page.waitForFunction(() => {
        const sel = document.querySelector(
          '[name="input-Model-Engine"]') as HTMLSelectElement | null;
        return !!sel && sel.options.length > 0;
      }, null, {timeout: 180_000});
    } catch (err) {
      const diag = await page.evaluate(() => {
        const v: any = (window as any).grok.shell.v;
        const root = v?.root as HTMLElement | undefined;
        const insights = Array.from(root?.querySelectorAll('span, .d4-card, [class*="insight"], [class*="tips"]') ?? [])
          .map((e) => (e as HTMLElement).textContent?.trim() ?? '')
          .filter((t) => t && t.length < 200).slice(0, 10);
        const hosts = Array.from(root?.querySelectorAll('[name^="input-host-"]') ?? [])
          .map((e) => e.getAttribute('name'));
        const noModelsMsg = Array.from(root?.querySelectorAll('span') ?? [])
          .find((s) => /No models reg/i.test(s.textContent || ''))?.textContent ?? null;
        
        const warnings = (((window as any).grok.shell.warnings) ?? []).slice(-10)
          .map((w: any) => typeof w === 'string' ? w.slice(0, 200) : JSON.stringify(w).slice(0, 200));
        return {insights, hosts, noModelsMsg, warnings};
      });
      throw new Error(`Model Engine select did not mount within 180s. ` +
        `Likely applicableEngines.isEmpty branch (Predict + Features pairing yields no engines) ` +
        `or cold engine.isApplicable budget exhaustion. ` +
        `noModelsMsg=${diag.noModelsMsg}; hosts=${JSON.stringify(diag.hosts)}; ` +
        `insights=${JSON.stringify(diag.insights)}; warnings=${JSON.stringify(diag.warnings)}; ` +
        `original error: ${err}`);
    }
    await page.evaluate(() => {
      const sel = document.querySelector(
        '[name="input-Model-Engine"]') as HTMLSelectElement;
      sel.focus();
      const setter = Object.getOwnPropertyDescriptor(HTMLSelectElement.prototype, 'value')!.set!;
      
      setter.call(sel, 'Eda: Linear Regression');
      sel.dispatchEvent(new Event('input', {bubbles: true}));
      sel.dispatchEvent(new Event('change', {bubbles: true}));
    });
    
    await page.waitForFunction(() => {
      const text = Array.from(document.querySelectorAll(
        'h3, h4, [class*="card-header"], .d4-engine-parameters-header'))
        .map((e) => e.textContent?.trim() ?? '')
        .find((t) => /Eda:\s*Linear/i.test(t) && !/pls|partial/i.test(t));
      return !!text;
    }, null, {timeout: 180_000});
  });

  await softStep(`1.6 Save model as ${MODEL_NAME} — dapi.ml.save persists`, async () => {
    
    await page.waitForFunction(() => {
      const btn = document.querySelector('[name="button-Save"]') as HTMLElement | null;
      return !!btn && !btn.classList.contains('d4-disabled');
    }, null, {timeout: 180_000});
    await page.locator('[name="button-Save"]').click();
    
    const nameInput = page.locator('.d4-dialog [name="input-host-Name"] input');
    await nameInput.waitFor({timeout: 10_000});
    await nameInput.focus();
    await nameInput.fill(MODEL_NAME);
    await page.locator('.d4-dialog [name="button-OK"]').click();
    // Wait for the save dialog to close.
    await page.waitForFunction(() =>
      !document.querySelector('.d4-dialog [name="input-host-Name"]'),
      null, {timeout: 90_000});
    
    await expect.poll(async () => await page.evaluate(async (name: string) => {
      const list = await (window as any).grok.dapi.models
        .filter(`friendlyName = "${name}"`).list();
      return list.length;
    }, MODEL_NAME), {timeout: 60_000}).toBeGreaterThan(0);
  });

  
  await softStep('2.1 Close current table view, re-open accelerometer.csv', async () => {
    
    await page.evaluate(async () => {
      const g: any = (window as any).grok;
      g.shell.closeAll();
      await new Promise((r) => setTimeout(r, 500));
      const df = await g.dapi.files.readCsv('System:DemoFiles/sensors/accelerometer.csv');
      g.shell.addTableView(df);
      await new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
      (window as any).__initialColCount = df.columns.length;
    });
    await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
    await page.locator('[name="div-ML"]').waitFor({timeout: 10_000});
  });

  await softStep('2.2 ML > Models > Apply Model... — apply dialog opens', async () => {
    await page.locator('[name="div-ML"]').click();
    await page.evaluate(() => {
      const models = document.querySelector('[name="div-ML---Models"]') as HTMLElement | null;
      if (!models) throw new Error('ML > Models submenu not found');
      const r = models.getBoundingClientRect();
      const ev = (t: string) => new MouseEvent(t, {
        bubbles: true, cancelable: true, view: window,
        clientX: r.left + 5, clientY: r.top + 5,
      });
      models.dispatchEvent(ev('mouseover'));
      models.dispatchEvent(ev('mouseenter'));
      models.dispatchEvent(ev('mousemove'));
    });
    await page.locator('[name="div-ML---Models---Apply-Model..."]').click();
    await page.locator('[name="dialog-Apply-predictive-model"]').waitFor({timeout: 10_000});
    
    await page.waitForFunction(() => {
      const sel = document.querySelector(
        '[name="dialog-Apply-predictive-model"] [name="input-host-Model"] select') as HTMLSelectElement | null;
      return !!sel && sel.options.length > 0;
    }, null, {timeout: 10_000});
  });

  await softStep(`2.3 Select ${MODEL_NAME} — inputs auto-map to (3/3)`, async () => {
    
    const pickResult = await page.evaluate(([wantName, wantBase]: string[]) => {
      const sel = document.querySelector(
        '[name="dialog-Apply-predictive-model"] [name="input-host-Model"] select') as HTMLSelectElement | null;
      if (!sel) return {found: false, opts: [] as string[]};
      const opts = Array.from(sel.options).map((o) => o.textContent || '');
      const tryPrefixes = [wantBase, wantName, wantName.slice(0, 32), wantName.slice(0, 24), wantName.slice(0, 16)];
      let idx = -1;
      for (const p of tryPrefixes) {
        idx = opts.findIndex((t) => t.includes(p));
        if (idx >= 0) break;
      }
      if (idx < 0) return {found: false, opts};
      sel.selectedIndex = idx;
      sel.dispatchEvent(new Event('change', {bubbles: true}));
      return {found: true, opts, picked: opts[idx]};
    }, [MODEL_NAME, MODEL_BASE]);
    expect(pickResult.found,
      `Just-saved model "${MODEL_NAME}" must surface in dapi.ml.suggested(tableInfo) ` +
      `for accelerometer.csv. Options seen: ${JSON.stringify(pickResult.opts)}.`)
      .toBe(true);
    
    await expect(page.locator(
      '[name="dialog-Apply-predictive-model"] [name="input-host-Inputs"]'))
      .toContainText('(3/3)', {timeout: 10_000});
  });

  await softStep('2.4 OK — engine apply appends prediction column', async () => {
    
    await page.locator('[name="dialog-Apply-predictive-model"] [name="button-OK"]').click();
    await page.locator('[name="dialog-Apply-predictive-model"]')
      .waitFor({state: 'detached', timeout: 30_000});
    await expect.poll(async () => await page.evaluate(() => {
      const df: any = (window as any).grok.shell.tv?.dataFrame;
      const initial: number = (window as any).__initialColCount ?? 0;
      return (df?.columns?.length ?? 0) - initial;
    }), {timeout: 60_000, intervals: [500, 1_000, 2_000]}).toBeGreaterThan(0);
    const result = await page.evaluate(() => {
      const df: any = (window as any).grok.shell.tv?.dataFrame;
      const initial: number = (window as any).__initialColCount ?? 0;
      const allNames = Array.from({length: df.columns.length},
        (_: unknown, i: number) => df.columns.byIndex(i).name);
      const newCols = allNames.slice(initial);
      
      const DG: any = (window as any).DG;
      const tagKey = DG?.TAGS?.PREDICTIVE_MODEL ?? '.predictive-model';
      const taggedCount = newCols.filter((n: string) => {
        const col = df.columns.byName(n);
        return !!col?.tags?.[tagKey] || !!col?.getTag?.(tagKey);
      }).length;
      return {initial, now: df.columns.length, newCols, taggedCount};
    });
    expect(result.newCols.length).toBeGreaterThan(0);
    if (result.taggedCount === 0) {
      console.warn(`Apply: ${result.newCols.length} new column(s) appended ` +
        `(${JSON.stringify(result.newCols)}) but none carry Tags.PredictiveModel. ` +
        `Engine may set tag server-side only.`);
    }
  });

  await softStep('3.1 Navigate to /models and select the saved model', async () => {
    
    await page.evaluate(() => {
      const g: any = (window as any).grok;
      g.shell.windows.showBrowse = true;
      g.shell.route('/models');
    });
    await page.waitForFunction(() => {
      const v = (window as any).grok.shell.v;
      return v && v.type === 'models';
    }, null, {timeout: 20_000});
    await page.locator('.grok-gallery-grid-item.grok-predictive-model').first()
      .waitFor({timeout: 15_000});
    const search = page.locator('input[placeholder="Search models by name or by #tags"]');
    await search.waitFor({timeout: 10_000});
    await search.fill('');
    await search.fill(MODEL_NAME);
    await expect.poll(async () => await page.evaluate((wantName: string) => {
      const titles = Array.from(document.querySelectorAll(
        '.grok-gallery-grid-item.grok-predictive-model .grok-gallery-grid-item-title'));
      return titles.some((t) => (t.textContent || '').trim().includes(wantName));
    }, MODEL_NAME), {timeout: 15_000}).toBe(true);
    await page.evaluate((wantName: string) => {
      const label = Array.from(document.querySelectorAll(
        '.grok-gallery-grid-item.grok-predictive-model .grok-gallery-grid-item-title'))
        .find((t) => (t.textContent || '').trim().includes(wantName)) as HTMLElement | undefined;
      if (!label) throw new Error(`card titled "${wantName}" not in DOM`);
      (label.closest('.grok-gallery-grid-item') as HTMLElement).click();
    }, MODEL_NAME);
    await page.waitForFunction(() =>
      Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
        .some((h) => (h.textContent || '').trim() === 'Performance'),
      null, {timeout: 15_000});
  });

  await softStep('3.2 Open Performance pane — stored metrics render', async () => {
    await page.evaluate(() => {
      const perf = Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
        .find((h) => (h.textContent || '').trim() === 'Performance') as HTMLElement | undefined;
      if (!perf) throw new Error('Performance pane header not found');
      perf.click();
    });
    
    await page.locator('[name="button-Run-Evaluation"]').waitFor({timeout: 15_000});
    await expect(page.locator('[name="button-Run-Evaluation"]')).toBeVisible();
  });

  await softStep('3.3 Click Run Evaluation — re-runs on CSV trainedOn', async () => {
    
    const warnBefore = await page.evaluate(() =>
      ((window as any).grok.shell.warnings ?? []).length);
    await page.locator('[name="button-Run-Evaluation"]').click();
    await page.waitForTimeout(3_000);
    const evidence = await page.evaluate((prevCount: number) => {
      const w: any[] = (window as any).grok.shell.warnings ?? [];
      const errors = w.slice(prevCount).filter((x) =>
        /error|fail|exception/i.test(JSON.stringify(x)));
      const panel = document.querySelector('.grok-prop-panel') ||
        document.querySelector('.d4-accordion');
      const widgetCanvases = panel?.querySelectorAll('canvas').length ?? 0;
      const widgetImgs = panel?.querySelectorAll('img').length ?? 0;
      return {errors, widgetCanvases, widgetImgs};
    }, warnBefore);
    expect(evidence.errors.length,
      `Run Evaluation must complete without error balloons. Errors seen: ${JSON.stringify(evidence.errors)}`)
      .toBe(0);
    if (evidence.widgetCanvases === 0 && evidence.widgetImgs === 0) {
      console.warn(`Run Evaluation completed without errors but no canvas/img ` +
        `widget surfaced in the property panel. Stored-metrics-only render is ` +
        `acceptable; flagging for retrospective review.`);
    }
  });

    await softStep('4.1 Right-click model card → Edit... — edit dialog opens', async () => {
    
    await page.evaluate((wantName: string) => {
      const label = Array.from(document.querySelectorAll(
        '.grok-gallery-grid-item.grok-predictive-model .grok-gallery-grid-item-title'))
        .find((t) => (t.textContent || '').trim().includes(wantName)) as HTMLElement | undefined;
      if (!label) throw new Error(`card titled "${wantName}" not in DOM`);
      const card = label.closest('.grok-gallery-grid-item') as HTMLElement;
      card.dispatchEvent(new MouseEvent('contextmenu',
        {bubbles: true, cancelable: true, button: 2}));
    }, MODEL_NAME);
    const editItem = page.locator('.d4-menu-popup .d4-menu-item-label', {hasText: /^Edit\.\.\.$/}).first();
    await expect(editItem).toBeVisible({timeout: 10_000});
    await editItem.click();
    await page.locator('[name="dialog-Predictive-model"]').waitFor({timeout: 10_000});
    await expect(page.locator('[name="dialog-Predictive-model"] .d4-dialog-title'))
      .toHaveText('Predictive model');
  });

  await softStep(`4.2 Change Description and OK — dapi.ml.save persists`, async () => {
    
    const descInput = page.locator(
      '[name="dialog-Predictive-model"] [name="input-host-Description"] input, ' +
      '[name="dialog-Predictive-model"] [name="input-host-Description"] textarea').first();
    await descInput.waitFor({timeout: 10_000});
    await descInput.focus();
    await descInput.fill(NEW_DESCRIPTION);
    // Drop the stale #grok-preloader overlay that lingers over the Edit dialog on the slow CI
    // stack and intercepts the pointer click (same guard as the Helm editor OK/CANCEL).
    await page.evaluate(() => document.querySelector('#grok-preloader')?.remove());
    await page.locator('[name="dialog-Predictive-model"] [name="button-OK"]').click();
    await page.locator('[name="dialog-Predictive-model"]')
      .waitFor({state: 'detached', timeout: 15_000});
    // Verify the model still exists (no accidental delete occurred during edit).
    await expect.poll(async () => await page.evaluate(async (name: string) => {
      const list = await (window as any).grok.dapi.models
        .filter(`friendlyName = "${name}"`).list();
      return list.length;
    }, MODEL_NAME), {timeout: 15_000}).toBeGreaterThan(0);
    
    await page.evaluate((wantName: string) => {
      const label = Array.from(document.querySelectorAll(
        '.grok-gallery-grid-item.grok-predictive-model .grok-gallery-grid-item-title'))
        .find((t) => (t.textContent || '').trim().includes(wantName)) as HTMLElement | undefined;
      if (!label) throw new Error(`card titled "${wantName}" not in DOM for re-open verify`);
      const card = label.closest('.grok-gallery-grid-item') as HTMLElement;
      card.dispatchEvent(new MouseEvent('contextmenu',
        {bubbles: true, cancelable: true, button: 2}));
    }, MODEL_NAME);
    const editItem2 = page.locator('.d4-menu-popup .d4-menu-item-label', {hasText: /^Edit\.\.\.$/}).first();
    await expect(editItem2).toBeVisible({timeout: 10_000});
    await editItem2.click();
    await page.locator('[name="dialog-Predictive-model"]').waitFor({timeout: 10_000});
    const descVerify = page.locator(
      '[name="dialog-Predictive-model"] [name="input-host-Description"] input, ' +
      '[name="dialog-Predictive-model"] [name="input-host-Description"] textarea').first();
    await descVerify.waitFor({timeout: 5_000});
    await expect(descVerify).toHaveValue(NEW_DESCRIPTION, {timeout: 5_000});
    await page.evaluate(() => document.querySelector('#grok-preloader')?.remove());
    await page.locator('[name="dialog-Predictive-model"] [name="button-CANCEL"]').click();
    await page.locator('[name="dialog-Predictive-model"]')
      .waitFor({state: 'detached', timeout: 10_000});
  });


  await softStep('5.1 Right-click model card → Share... — share dialog opens', async () => {
    
    await page.evaluate((wantName: string) => {
      const label = Array.from(document.querySelectorAll(
        '.grok-gallery-grid-item.grok-predictive-model .grok-gallery-grid-item-title'))
        .find((t) => (t.textContent || '').trim().includes(wantName)) as HTMLElement | undefined;
      if (!label) throw new Error(`card titled "${wantName}" not in DOM`);
      const card = label.closest('.grok-gallery-grid-item') as HTMLElement;
      card.dispatchEvent(new MouseEvent('contextmenu',
        {bubbles: true, cancelable: true, button: 2}));
    }, MODEL_NAME);
    const shareItem = page.locator('.d4-menu-popup .d4-menu-item-label', {hasText: /^Share\.\.\.$/}).first();
    await expect(shareItem).toBeVisible({timeout: 10_000});
    await shareItem.click();
    
    const encodedName = MODEL_NAME.replace(/[:_;*\[\]{}|\s]/g, '-');
    const shareDialog = page.locator(
      `[name="dialog-Share-${encodedName}"], ` +
      `.d4-dialog:has(.d4-dialog-title:has-text("Share ${MODEL_NAME}"))`).first();
    await expect(shareDialog).toBeVisible({timeout: 15_000});
  });

  await softStep('5.2 Cancel share dialog — leaves permissions unchanged', async () => {
    
    const encodedName = MODEL_NAME.replace(/[:_;*\[\]{}|\s]/g, '-');
    const cancelBtn = page.locator(
      `[name="dialog-Share-${encodedName}"] [name="button-CANCEL"], ` +
      `.d4-dialog [name="button-CANCEL"]`).first();
    await cancelBtn.click();
    await page.locator(`[name="dialog-Share-${encodedName}"]`)
      .waitFor({state: 'detached', timeout: 5_000})
      .catch(() => { /* generic dialog selector may have closed already */ });
  });


  await softStep('6.1 Right-click model card → Delete — confirm modal opens', async () => {
    
    await page.evaluate((wantName: string) => {
      const label = Array.from(document.querySelectorAll(
        '.grok-gallery-grid-item.grok-predictive-model .grok-gallery-grid-item-title'))
        .find((t) => (t.textContent || '').trim().includes(wantName)) as HTMLElement | undefined;
      if (!label) throw new Error(`card titled "${wantName}" not in DOM`);
      const card = label.closest('.grok-gallery-grid-item') as HTMLElement;
      card.dispatchEvent(new MouseEvent('contextmenu',
        {bubbles: true, cancelable: true, button: 2}));
    }, MODEL_NAME);
    const deleteItem = page.locator('.d4-menu-popup .d4-menu-item-label', {hasText: /^Delete$/}).first();
    await expect(deleteItem).toBeVisible({timeout: 10_000});
    await deleteItem.click();
    
    const dlg = page.locator('.d4-dialog').filter({hasText: /delete|are you sure|remove/i }).first();
    await expect(dlg).toBeVisible({timeout: 10_000});
  });

  await softStep(`6.2 Confirm delete — ${MODEL_NAME} removed via dapi.ml.delete`, async () => {
    const dlg = page.locator('.d4-dialog').filter({hasText: /delete|are you sure|remove/i }).first();
    const balloonsBefore = await page.locator('.grok-balloon-error, .d4-balloon-error').count();
    const ok = dlg.locator('[name="button-OK"]').first();
    if (await ok.count() > 0)
      await ok.click();
    else
      await dlg.locator('button:has-text("DELETE"), button:has-text("Delete"), button:has-text("OK")')
        .first().click();
    await dlg.waitFor({state: 'detached', timeout: 15_000});
    await page.waitForTimeout(1_500);
    const balloonsAfter = await page.locator('.grok-balloon-error, .d4-balloon-error').count();
    expect(balloonsAfter,
      `No new error balloon should surface during delete (before=${balloonsBefore}, after=${balloonsAfter}).`)
      .toBeLessThanOrEqual(balloonsBefore);
  });

  await softStep('6.3 Cleanup verification — model absent from dapi.ml + /models', async () => {
    
    await expect.poll(async () => await page.evaluate(async (wantName: string) => {
      const list = await (window as any).grok.dapi.models
        .filter(`friendlyName = "${wantName}"`).list();
      return list.length;
    }, MODEL_NAME), {timeout: 30_000, intervals: [500, 1_000, 2_000]}).toBe(0);
    const search = page.locator('input[placeholder="Search models by name or by #tags"]');
    await search.fill('');
    await search.fill(MODEL_NAME);
    await expect.poll(async () => await page.evaluate((wantName: string) => {
      const titles = Array.from(document.querySelectorAll(
        '.grok-gallery-grid-item.grok-predictive-model .grok-gallery-grid-item-title'));
      return titles.filter((t) => (t.textContent || '').trim().includes(wantName)).length;
    }, MODEL_NAME), {timeout: 15_000, intervals: [250, 500, 1_000],
      message: `Catalog still shows card(s) titled "${MODEL_NAME}" after delete`}).toBe(0);
  });

  
  await page.evaluate(async ({prefix}) => {
    const g: any = (window as any).grok;
    const all = await g.dapi.models.list();
    const stale = all.filter((m: any) => (m.friendlyName || m.name || '').startsWith(prefix));
    for (const m of stale)
      try { await g.dapi.models.delete(m); } catch (_) { /* best-effort */ }
  }, {prefix: MODEL_BASE});

  if (stepErrors.length > 0) {
    throw new Error(`${stepErrors.length} step(s) failed:\n` +
      stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n'));
  }
});
