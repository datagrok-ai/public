/* ---
sub_features_covered: [models.command.train, models.command.apply, models.command.delete, models.view.training, models.view.browser, models.workflow.apply-dialog, models.workflow.remove, models.engines.api.apply, models.engines.package, models.api.save, models.api.run]
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: absent (scenario coverage_type: smoke; section ui-smoke role per Notes)
//   sub_features_covered: see block above
//   ui_coverage_responsibility: section ui-smoke owner — Train/Apply/Apply-on-new-dataset/Delete UI flows DOM-driven
//   related_bugs: [GROK-18612, GROK-846, GROK-19177, GROK-2381, GROK-873, GROK-19550]
// Selectors are class-1 (all [name=...] strings appear verbatim in
// .claude/skills/grok-browser/references/models.md: div-ML, div-ML---Models---Train-Model...,
// div-ML---Models---Apply-Model..., dialog-Select-columns..., button-Save, button-OK,
// input-host-Features, input-host-Model, input-host-Inputs, input-Table, input-Model).
// MCP recon 2026-06-09 (chrome-devtools list_pages + evaluate_script on dev.datagrok.ai)
// confirmed div-ML, div-ML---Models---Train-Model..., div-ML---Models---Apply-Model... present
// after addTableView; Tools > Dev > Open Test Dataset addressed via grok.shell.topMenu.find('Tools').root.
//
// Retry-round-1 hypothesis (Gate B FAIL B-RUN-PASS,B-STAB-01 cascade from step 1.7):
//   Round-1 fixes were tactical waits + dialog scoping. They did NOT address
//   the real failure root at step 1.3.
//
// Retry-round-2 hypothesis (cycle 2026-06-09-models-automate-01; Gate B FAIL
// again, attempts 1-3 all failed at step 1.3 with counter still "4 checked"
// instead of "3 checked", cascading downstream):
//   Round-2 prior author claimed "synthetic dispatchEvent doesn't toggle the
//   column-grid overlay; pivot to page.mouse.click(x,y)". This was an
//   under-investigated hypothesis: the picked coordinates were wrong, not the
//   event mechanism. The prior spec computed
//     stride = canvas.height / 4 = 35.5
//     y = top + stride * 0.5 = top + ~17.75 (HEADER strip, not row 0)
//     x = left + 15 (left edge — checkboxes sit at right edge)
//   Header strip top..top+22 is a select-all toggle (not row deselect), AND
//   the checkbox column sits at right-40, NOT left+15.
//
// Retry-round-4 hypothesis (cycle 2026-06-09-models-automate-02 hypothesis_retry,
// this dispatch), backed by fresh live MCP recon on dev.datagrok.ai 2026-06-10
// on a clean accelerometer.csv Train view:
//
//   Recon disambiguated two failure surfaces from the Gate B FAIL cascade.
//
//   Step 1.3 (column picker): the round-3 spec's coordinate-correction
//   (synthetic dispatch at right - 40, top + 36) was reconfirmed to work
//   under live MCP — label-All yields "4 checked" → synthetic dispatch on
//   overlay canvas at the right-edge checkbox row 0 → "3 checked" → OK →
//   Features reads "(3) accel_y, accel_z, time_offset". The picker DOES
//   present all 4 columns; accel_x is NOT auto-excluded from it (an
//   intermediate hypothesis that the picker showed only 3 rows was
//   falsified by recon). Round-3's tactical fix is correct as authored.
//
//   Step 4.1 (Browse > Models gallery): the round-1/2 spec's Browse-tree
//   path FAILED because on this build the Platform group is expanded by
//   default — pre-emptively clicking it COLLAPSES Platform and hides
//   "Predictive models", which was the attempt-3 failure mode. Fresh recon
//   confirms grok.shell.route('/models') opens PredictiveModelsView with
//   grok.shell.v.type === 'models' (also confirmed in
//   grok-browser/references/models.md under "/models — PredictiveModelsView
//   browser"). This dispatch keeps the route-based step 4.1.
//
//   Net change vs round-3: step 4.1 only. Step 1.3 keeps round-3's
//   coordinate-correction synthetic dispatch.
//
//   Steps 4.2/4.3/4.4 selector recon: gallery card structure verified
//   2026-06-10 — .grok-gallery-grid-item-title is a <label> child whose
//   closest('.grok-gallery-grid-item') resolves to the card div as
//   expected. Their prior failures were entirely cascade-driven (no cards
//   existed because no models were saved). Once 1.3 + 4.1 are fixed, the
//   gallery handlers work as authored.
//
// Round-5 tactical fixes (cycle 2026-06-10-models-automate-02, initial dispatch,
// MCP unavailable — session expired; cheap-checks only):
//
//   Fix A — Step 4.2 Activity pane assertion removed.
//     models.md recon note (DOM 2026-06-10) states the context panel accordion
//     shows "Details / Performance / Sharing / Chats / Sticky meta" with NO
//     Activity pane [DOM models.md §"Right-side context panel"]. The prior spec
//     asserted Activity — this is a confirmed wrong assertion (atlas drift, not
//     a spec design choice). Removed; now asserts only Details + Performance +
//     Sharing which are all confirmed by models.md DOM evidence.
//
//   Fix B — Step 3.1 Tools > Dev menu replaced with grok.data.testData().
//     Scenario Notes section explicitly documents: "For deterministic CI, the
//     equivalent JS-API path is grok.data.testData('random walk', 1000, 10); the
//     migrated spec should prefer that path if dev-mode is not guaranteed." The
//     prior spec used the Tools > Dev menu (uncertain in CI). B-STAB-01 is
//     consistent with non-deterministic menu availability. Replaced with the
//     sanctioned JS-API path per scenario Notes authority.
//
// Selector-recon notes (class-2: live MCP-observed 2026-06-10, not yet in
// grok-browser/references/models.md):
//   grok.shell.route('/models') — opens PredictiveModelsView with
//     grok.shell.v.type === 'models' (refdoc-confirmed under "/models —
//     PredictiveModelsView browser"); used here in lieu of Browse tree
//     navigation because the tree's Platform group expansion state is
//     not stable across cycles.
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Predictive models: Train / Apply / Apply on new dataset / Delete', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);

  // Pre-test cleanup: remove stale models from previous runs so gallery assertions are deterministic.
  await page.evaluate(async () => {
    const g: any = (window as any).grok;
    for (const name of ['Accelerometer_model_PLS', 'Accelerometer_model_LR']) {
      const list = await g.dapi.models.filter(`friendlyName = "${name}"`).list();
      for (const m of list) await g.dapi.models.delete(m);
    }
  });

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

  await softStep('1.3 Set Features to accel_y, accel_z, time_offset', async () => {
    // Open the column-picker via the Features editor (class-1: [name="div-Features"]).
    await page.locator('[name="div-Features"]').click();
    await page.locator('.d4-dialog[name="dialog-Select-columns..."]').waitFor({timeout: 10_000});
    // Settle wait — column-grid canvas paints async (~1.2s on dev).
    await page.waitForTimeout(1200);
    // Click "All" → 4 columns checked, then deselect accel_x (row 0).
    // The column-grid is canvas-only; the checkbox column sits at the RIGHT
    // edge of the overlay canvas, not the left. Synthetic MouseEvent dispatch
    // at (right - 40, top + 36) toggles row 0 — verified live MCP 2026-06-10.
    await page.locator('.d4-dialog [name="label-All"]').click();
    await page.waitForFunction(() => {
      const lbl = Array.from(document.querySelectorAll('.d4-dialog label'))
        .find((l) => /checked/.test(l.textContent || ''));
      return lbl?.textContent?.startsWith('4');
    }, null, {timeout: 10_000});
    // Synthetic dispatch on overlay canvas at right-edge row-0 checkbox.
    const toggled = await page.evaluate(async () => {
      const dlg = document.querySelector('.d4-dialog[name="dialog-Select-columns..."]')!;
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
    // Defence in depth: if synthetic dispatch left counter at "4 checked"
    // (e.g. future isTrusted gating), retry once with page.mouse.click which
    // dispatches isTrusted events.
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
  });

  // Engine select only appears after a missing-value strategy is chosen (or after initial auto-train
  // on a complete dataset). Wait up to 5s; if still absent, tick Ignore missing to unlock it.
  const engineVisible = await page.waitForFunction(
    () => !!document.querySelector('[name="input-Model-Engine"]'), null, {timeout: 5_000})
    .then(() => true).catch(() => false);
  if (!engineVisible) {
    await page.locator('[name="input-Ignore-missing"]').click();
    await page.waitForFunction(
      () => !!document.querySelector('[name="input-Model-Engine"]'), null, {timeout: 15_000});
  }

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
    // Wait for the PLS-specific Components parameter to render — it appears after the engine
    // preview card fully loads, which is async after the select change.
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
    // Wait for Save button to be enabled (model preview must be ready, d4-disabled cleared).
    await page.waitForFunction(() => {
      const btn = document.querySelector('[name="button-Save"]') as HTMLElement | null;
      return !!btn && !btn.classList.contains('d4-disabled');
    }, null, {timeout: 30_000});
    await page.locator('[name="button-Save"]').click();
    // Sibling-spec convention: target the Name input via its host wrapper +
    // trailing input (the host is the wrapper div, the actual <input> is nested).
    const nameInput = page.locator('.d4-dialog [name="input-host-Name"] input');
    await nameInput.waitFor({timeout: 10_000});
    await nameInput.focus();
    await nameInput.fill('Accelerometer_model_PLS');
    await page.locator('.d4-dialog [name="button-OK"]').click();
    // Use the save-dialog-specific selector for close detection (not generic .d4-dialog)
    // to avoid false closures if any unrelated dialog is open.
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
    // After Save dialog closes, the Predictive model view re-renders and the engine
    // <select> is briefly detached. Wait for the element to be present AND populated
    // before manipulating it (this was the cascade origin in the prior Gate B run).
    await page.waitForFunction(() => {
      const sel = document.querySelector('[name="input-Model-Engine"]') as HTMLSelectElement | null;
      return !!sel && sel.options.length > 0;
    }, null, {timeout: 30_000});
    // Settle wait to let any re-render queue drain before we set the value.
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
    // Same enable-gate as 1.6 — d4-disabled must clear once the LR preview is ready.
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

  // ── Scenario 2: Apply ────────────────────────────────────────────────────
  // Boundary cleanup: close any stale dialogs that may be lingering from scenario 1
  // before opening a fresh view (prevents strict-mode .d4-dialog violations downstream).
  await page.evaluate(() => {
    document.querySelectorAll('.d4-dialog').forEach(d => {
      const cancel = d.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
      if (cancel) cancel.click(); else (d as HTMLElement).remove();
    });
  });
  await page.waitForFunction(() => !document.querySelector('.d4-dialog'), null, {timeout: 5_000})
    .catch(() => {/* best-effort */});
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
    // Scope to the dialog by its name (sibling apply-spec.ts pattern, class-1) so
    // any stragglers from prior steps don't trigger strict-mode violations.
    await page.locator('[name="dialog-Apply-predictive-model"]').waitFor({timeout: 10_000});
    // Wait for the model <select> to be populated before reading options.
    await page.waitForFunction(() => {
      const sel = document.querySelector(
        '[name="dialog-Apply-predictive-model"] [name="input-host-Model"] select') as HTMLSelectElement | null;
      return !!sel && sel.options.length > 0;
    }, null, {timeout: 10_000});
    await page.evaluate(async () => {
      const sel = document.querySelector(
        '[name="dialog-Apply-predictive-model"] [name="input-host-Model"] select') as HTMLSelectElement;
      const options = Array.from(sel.options);
      // Select PLS model by name (not by position — timestamp ordering is fragile)
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
      // Select LR model by name (not by position — timestamp ordering is fragile)
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

  // ── Scenario 3: Apply on a new dataset ───────────────────────────────────
  await softStep('3.1 Open random walk test dataset (1000 rows, 10 cols)', async () => {
    // Scenario Notes: "For deterministic CI, the equivalent JS-API path is
    // grok.data.testData('random walk', 1000, 10); the migrated spec should
    // prefer that path if dev-mode is not guaranteed." Using the sanctioned
    // JS-API path instead of Tools > Dev > Open Test Dataset which requires
    // dev-mode enabled (non-deterministic in CI — identified as B-STAB-01 root).
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
      // Select LR model by name (not by index position which is timestamp-order-dependent)
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
    // Random walk dataset starts with 10 columns (#0..#9). After applying
    // the LR model via isApplicable Levenshtein/JaroWinkler input matching,
    // at least one prediction column is appended — scenario expects "the
    // inputs and the resulting prediction column appear".
    expect(cols.length).toBeGreaterThan(10);
  });

  // ── Scenario 4: Delete ───────────────────────────────────────────────────
  await softStep('4.1 Browse > Platform > Predictive models', async () => {
    // Use the route helper to open the PredictiveModelsView directly. The
    // Browse tree's "Platform" group expansion state is not stable across
    // cycles (it was expanded by default in the 2026-06-10 recon, so a
    // pre-emptive .click() COLLAPSES it and hides Predictive models — the
    // attempt-3 failure mode). grok.shell.route('/models') is the canonical
    // entry per models.md `/models — PredictiveModelsView browser`.
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
    // Activity pane is absent from the models browser accordion in this build
    // (models.md DOM 2026-06-10: "Details / Performance / Sharing / Chats / Sticky meta").
    // Prior spec asserted Activity — confirmed incorrect (atlas drift).
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
    // Navigate away and back to force gallery re-query so PLS card is visible for step 4.4.
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
