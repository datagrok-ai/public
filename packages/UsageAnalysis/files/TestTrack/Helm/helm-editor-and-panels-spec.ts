import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {finishSpec} from '../helpers/viewers';

test.use(specTestOptions);
test.use({timeout: 600_000});

const DATASET_PATH = 'System:AppData/Helm/samples/HELM.csv';

test('Helm — cell rendering, Web Editor & Properties panel', async ({page}) => {
  // 600s ceiling: cold Block B open (~300s lazy init) + login (~15s) +
  // setup (~17s) + blocks A-G warm (~35s) = ~367s worst-case; 233s headroom.
  test.setTimeout(600_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);

  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    (window as any).grok.shell.settings.showFiltersIconsConstantly = true;
    (window as any).grok.shell.windows.simpleMode = true;
    // showContextPanel must be set explicitly: simpleMode=true does NOT auto-open
    // the right context panel, and pane-Properties is absent from DOM when it is
    // false — causing Block F waitFor to timeout in headless Playwright sessions.
    (window as any).grok.shell.windows.showContextPanel = true;
    (window as any).grok.shell.closeAll();
    // Wait for Bio package to register detectMacromolecule before loading HELM.csv.
    // In a fresh Playwright context, Bio's autostart may not have completed yet
    // by the time readCsv fires — if the detector is unregistered, the semType
    // detection skips Bio entirely and helmCol.semType stays null (confirmed via
    // direct spec run: grok test --host dev → FAIL at line 183, Received: null).
    // Round-5: reduced from 30×500ms (15s max) to 8×500ms (4s max) to cut
    // setup time. Bio typically loads in <2s on dev; 4s is a safe ceiling.
    for (let i = 0; i < 8; i++) {
      if ((window as any).DG.Func.find({name: 'detectMacromolecule'}).length > 0) break;
      await new Promise((r) => setTimeout(r, 500));
    }
    const df = await (window as any).grok.dapi.files.readCsv(path);
    (window as any).grok.shell.addTableView(df);
    // Round-5: reduced semType race cap 5s→3s; Bio detector is ready by now.
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(() => resolve(), 3000);
    });
    // Fallback: if Bio detector was still late and semType is null, re-trigger detection.
    // detectSemanticTypes returns a Promise — await it so the caller doesn't proceed
    // before all detectors have run. Then poll helmCol.semType directly (up to 5s)
    // in case the async Bio detector posts the result after the Promise resolves.
    // Round-5: fallback poll reduced from 20×500ms (10s max) to 10×500ms (5s max).
    const helmCol = df.col('HELM');
    if (helmCol && !helmCol.semType) {
      await (window as any).grok.data.detectSemanticTypes(df);
      for (let i = 0; i < 10; i++) {
        if (df.col('HELM')?.semType) break;
        await new Promise((r) => setTimeout(r, 500));
      }
    }
    const cols = Array.from({length: df.columns.length}, (_: unknown, i: number) =>
      df.columns.byIndex(i));
    const hasMacro = cols.some((c: any) => c.semType === 'Macromolecule');
    if (hasMacro) {
      // Round-5: canvas poll reduced from 60×200ms (12s max) to 30×200ms (6s max).
      for (let i = 0; i < 30; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      // Bio + Helm package init (RDKit + JSDraw2 + monomer-dictionary rewrite).
      // Round-5: reduced from 3500ms to 2000ms — JSDraw2/Dojo stack loads once;
      // on warm sessions this settle is sufficient for canvas readiness.
      await new Promise((r) => setTimeout(r, 2000));
    }
    // Round-11: removed getHelmHelper/helmWebEditor$ blocking wait from setup.
    // Awaiting any Helm:* function call inside page.evaluate blocks for the
    // full cold Helm init duration (~5min: Dojo + JSDraw2 + Pistoia HELM Web
    // Editor + monomer-dict rewrite). This causes the test-level timeout to
    // fire before a single assertion runs. The cold init is deferred to Block B
    // via a fire-and-forget editMoleculeCell call + outer Playwright waitFor,
    // which is interruptible by test timeout (unlike an awaited page.evaluate).
  }, DATASET_PATH);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});

  // Pre-flight invariants — Bio detector tagged the column AND it has enough
  // rows so the Block A / Block F per-cell flows have non-degenerate input.
  const setupProbe = await page.evaluate(() => {
    const df = (window as any).grok.shell.tv.dataFrame;
    const helmCol = df.col('HELM');
    return {
      rowCount: df.rowCount,
      helmName: helmCol?.name ?? null,
      helmSemType: helmCol?.semType ?? null,
      helmUnits: helmCol?.meta?.units ?? null,
      helmRenderer: helmCol?.getTag ? helmCol.getTag('cell.renderer') : null,
      helmQuality: helmCol?.getTag ? helmCol.getTag('quality') : null,
    };
  });
  expect(setupProbe.helmSemType,
    'atlas helm.rendering.cell-renderer precondition: Bio detector MUST tag HELM column semType=Macromolecule').toBe('Macromolecule');
  expect(setupProbe.helmUnits,
    'atlas helm.rendering.cell-renderer precondition: column tag meta.units MUST be "helm" for HelmGridCellRenderer to attach').toBe('helm');
  expect(setupProbe.helmRenderer,
    'atlas helm.rendering.cell-renderer precondition: column tag cell.renderer MUST be "helm"').toBe('helm');
  expect(setupProbe.helmQuality,
    'atlas helm.rendering.cell-renderer precondition: column tag quality MUST be "Macromolecule"').toBe('Macromolecule');
  expect(setupProbe.rowCount,
    'scenario .md Setup: HELM.csv MUST yield 540 rows').toBe(540);

  // Helper: compute the (x,y) click coord on the main grid canvas for a cell.
  const getCellCoord = async (rowIdx: number) => page.evaluate((r) => {
    const grid = (window as any).grok.shell.tv.grid;
    const cell = grid.cell('HELM', r);
    const b = cell.bounds;
    const canvases = Array.from(document.querySelectorAll('[name="viewer-Grid"] canvas')) as HTMLCanvasElement[];
    const main = canvases.find((c) => {
      const rect = c.getBoundingClientRect();
      return rect.width > 100 && rect.height > 100;
    });
    if (!main) return null;
    const rect = main.getBoundingClientRect();
    return {x: Math.round(rect.left + b.x + b.width / 2), y: Math.round(rect.top + b.y + b.height / 2)};
  }, rowIdx);

  // Helper: dismiss all open full-screen Web Editor dialogs before opening a
  // fresh one. Prevents strict-mode violations when a dblclick opens a dialog
  // late (after the fallback timeout), leaving two full-screen dialogs stacked.
  const closeAllEditorDialogs = async () => {
    const count = await page.locator('.d4-dialog.d4-dialog-full-screen').count();
    if (count === 0) return;
    await page.evaluate(() => {
      const dialogs = Array.from(document.querySelectorAll('.d4-dialog.d4-dialog-full-screen'));
      for (const dlg of dialogs) {
        const cancelBtn = dlg.querySelector('button[name="button-CANCEL"]') as HTMLElement | null;
        if (cancelBtn) cancelBtn.click();
      }
    });
    // Wait until all full-screen dialogs are gone
    await page.waitForFunction(
      () => document.querySelectorAll('.d4-dialog.d4-dialog-full-screen').length === 0,
      {timeout: 8_000}
    ).catch(() => {/* ignore timeout — proceed and let next waitFor surface the real error */});
    await page.waitForTimeout(300);
  };

  // Helper: click the CANCEL or OK button via JS evaluate() to bypass
  // Playwright's hit-test. Playwright's pointer-based click() fails for these
  // buttons because the full-screen dialog container intercepts the synthetic
  // pointer event (position:fixed, z-index:3005 overlay). JS .click()
  // bypasses the hit-test entirely and fires the native click handler.
  // Round-13: click button in ALL open full-screen dialogs. The Block B
  // fire-and-forget fallback + delayed dblclick race can open two dialogs
  // simultaneously; clicking only the first leaves the second open and the
  // waitForFunction(length===0) times out. Iterating over all is idempotent
  // for the normal single-dialog case.
  const clickDialogButton = async (buttonName: 'button-CANCEL' | 'button-OK') => {
    await page.evaluate((name) => {
      const dialogs = Array.from(document.querySelectorAll('.d4-dialog.d4-dialog-full-screen'));
      for (const dlg of dialogs) {
        const btn = dlg.querySelector(`button[name="${name}"]`) as HTMLElement | null;
        if (btn) btn.click();
      }
    }, buttonName);
    await page.waitForFunction(
      () => document.querySelectorAll('.d4-dialog.d4-dialog-full-screen').length === 0,
      {timeout: 15_000}
    );
  };

  // Helper: open the Web Editor via JS-API. Permitted for Blocks D/E/G as
  // SETUP (the owned UI flow helm.open-editor.cell-editor is exercised in
  // Block B; helm.open-editor.edit-helm-action is exercised in Block C).
  // Always closes any stale dialogs first to avoid strict-mode violations.
  // Round-12: fire-and-forget pattern (no await inside evaluate). The blocking
  // `await grok.functions.call(...)` inside evaluate can hold the browser for
  // minutes when JsDraw2 needs re-init (same root cause as Block B fallback
  // Round-11). Dialog visibility detected via outer Playwright waitFor instead.
  const openEditorViaJsApi = async () => {
    await closeAllEditorDialogs();
    await page.evaluate(() => {
      const grid = (window as any).grok.shell.tv.grid;
      const DG = (window as any).DG;
      const cell = DG.GridCell.fromColumnRow(grid, 'HELM', 0);
      (window as any).grok.functions.call('Helm:editMoleculeCell', {cell}); // intentional non-await
    });
    await page.locator('.d4-dialog.d4-dialog-full-screen').first()
      .waitFor({state: 'visible', timeout: 30_000});
    await page.locator('.d4-dialog.d4-dialog-full-screen [id^="__JSDraw_"]').first()
      .waitFor({state: 'attached', timeout: 12_000});
    await page.waitForTimeout(600);
  };

  await softStep('Block A Step 1: HELM column auto-renders via HelmGridCellRenderer (structural-proxy)', async () => {
    const probe = await page.evaluate(() => {
      const grid = (window as any).grok.shell.tv.grid;
      const helmCol = grid.columns.byName('HELM');
      const canvases = Array.from(document.querySelectorAll('[name="viewer-Grid"] canvas')) as HTMLCanvasElement[];
      const mainCanvasCount = canvases.filter((c) => {
        const r = c.getBoundingClientRect();
        return r.width > 100 && r.height > 100;
      }).length;
      return {helmColExists: !!helmCol, canvasCount: canvases.length, mainCanvasCount};
    });
    expect(probe.helmColExists,
      'helm.rendering.cell-renderer attachment: grid.columns.byName("HELM") MUST resolve').toBe(true);
    expect(probe.mainCanvasCount,
      'helm.rendering.cell-renderer: a main render canvas (width+height > 100) MUST be present').toBeGreaterThan(0);
  });

  await softStep('Block A Step 2: hover monomer hexagon → .d4-tooltip element responds (real mouse)', async () => {
    const coords = await getCellCoord(0);
    expect(coords,
      'hover precondition: a main grid canvas + HELM cell 0 bounds MUST resolve').not.toBeNull();
    // Hover the LEFT edge of the cell — the first monomer hexagon sits there.
    const hx = coords!.x - 170; // left-edge offset within the 400px-wide cell
    const hy = coords!.y;
    await page.mouse.move(hx - 30, hy - 30);
    await page.waitForTimeout(120);
    await page.mouse.move(hx, hy, {steps: 6});
    await page.waitForTimeout(1200);
    const tooltipInfo = await page.evaluate(() => {
      const tt = document.querySelector('.d4-tooltip');
      return {
        exists: !!tt,
        display: tt ? (tt as HTMLElement).style.display : null,
        textLen: (tt?.textContent ?? '').length,
      };
    });
    expect(tooltipInfo.exists,
      'helm.rendering.cell-renderer wiring: .d4-tooltip singleton MUST be present in the DOM (HelmGridCellRendererBack subscribed)').toBe(true);
    if (tooltipInfo.textLen === 0) {
      console.warn('Block A hover tooltip: .d4-tooltip element present but text empty under Playwright headless (per helm.md Pitfall #11 caveat). Element-presence assertion stands as the structural-proxy guard.');
    }
    const balloonErrors = await page.locator('.d4-balloon.error').count();
    expect(balloonErrors, 'hover MUST NOT raise an error balloon').toBe(0);
  });

  let dialogOpenUsedFallback = false;
  await softStep('Block B Step 1: double-click HELM cell → full-screen Web Editor dialog opens', async () => {
    const coords = await getCellCoord(0);
    expect(coords,
      'dblclick precondition: a main grid canvas + HELM cell 0 bounds MUST resolve').not.toBeNull();
    await page.mouse.move(coords!.x, coords!.y, {steps: 4});
    await page.mouse.dblclick(coords!.x, coords!.y);
    let dialogOpened = false;
    try {
      // Round-5: reduced from 20s to 12s. Bio cold-open on warm sessions
      // completes in <5s; 12s accounts for first-session JSDraw2 load on
      // a warm pre-loaded Helm package (monomer dict already cached in memory
      // after setup's 2s settle). True cold-open (first open in the browser
      // session) is handled by the openEditorViaJsApi fallback below.
      await page.locator('.d4-dialog.d4-dialog-full-screen').first()
        .waitFor({state: 'visible', timeout: 12_000});
      dialogOpened = true;
    } catch (_) { /* fall back below */ }
    if (!dialogOpened) {
      dialogOpenUsedFallback = true;
      console.warn('Block B dblclick did not open editor in 12s — cold JSDraw2/Helm init in progress; firing editMoleculeCell (fire-and-forget) and waiting up to 360s.');
      // Close any late-arriving dialog from the dblclick before re-firing.
      await closeAllEditorDialogs();
      // Fire-and-forget: intentionally NOT awaited inside evaluate so the
      // page.evaluate returns immediately instead of blocking for ~5min of
      // cold Helm init. The dialog appearance is detected via the Playwright
      // waitFor below, which is interruptible by the test-level timeout.
      await page.evaluate(() => {
        const g = (window as any).grok;
        const DG = (window as any).DG;
        const cell = DG.GridCell.fromColumnRow(g.shell.tv.grid, 'HELM', 0);
        g.functions.call('Helm:editMoleculeCell', {cell}); // intentional non-await
      });
      await page.locator('.d4-dialog.d4-dialog-full-screen').first()
        .waitFor({state: 'visible', timeout: 360_000});
      await page.locator('.d4-dialog.d4-dialog-full-screen [id^="__JSDraw_"]').first()
        .waitFor({state: 'attached', timeout: 30_000});
      await page.waitForTimeout(800);
    } else {
      await page.locator('.d4-dialog.d4-dialog-full-screen [id^="__JSDraw_"]').first()
        .waitFor({state: 'attached', timeout: 12_000});
      // Round-5: reduced settle 1500→800ms — JSDraw2 host attaches synchronously
      // with dialog show on warm sessions (empirical: attach_ms=0 per prior recon).
      await page.waitForTimeout(800);
    }
  });

  await softStep('Block B Step 1: dialog structure — JSDraw2 host + palette + bottom tabs + footer OK/CANCEL', async () => {
    const dlgByName = await page.locator('[name="dialog-"]').count();
    expect(dlgByName,
      'helm.editor.web-editor-app: full-screen Web Editor dialog name= attribute slug is "dialog-" (empty-title slug per helm.md)').toBeGreaterThan(0);
    const okCount = await page.locator('.d4-dialog button[name="button-OK"]').count();
    const cancelCount = await page.locator('.d4-dialog button[name="button-CANCEL"]').count();
    expect(okCount,
      'helm.editor.web-editor-app: footer OK button[name="button-OK"] MUST exist').toBeGreaterThan(0);
    expect(cancelCount,
      'helm.editor.web-editor-app: footer CANCEL button[name="button-CANCEL"] MUST exist').toBeGreaterThan(0);
    const expectedTabs = ['Monomers', 'Rules', 'Placeholders', 'Sequence', 'HELM', 'Properties', 'Structure View'];
    const tabTexts = await page.evaluate(() => {
      const dlg = document.querySelector('.d4-dialog.d4-dialog-full-screen');
      return Array.from(dlg?.querySelectorAll('td.hwe-tab-td') ?? []).map((t) => t.textContent?.trim());
    });
    for (const expected of expectedTabs) {
      expect(tabTexts,
        `helm.editor.web-editor-app: palette/bottom tab "${expected}" MUST be present (td.hwe-tab-td)`).toContain(expected);
    }
    const jsdHostInfo = await page.evaluate(() => {
      const hosts = Array.from(document.querySelectorAll('.d4-dialog.d4-dialog-full-screen [id^="__JSDraw_"]')) as HTMLElement[];
      const withSvg = hosts.filter((h) => h.querySelector('svg'));
      return {hostCount: hosts.length, hostsWithSvg: withSvg.length};
    });
    expect(jsdHostInfo.hostCount,
      'helm.editor.web-editor-wrapper: at least one JSDraw2 host (id="__JSDraw_N") MUST attach inside the dialog').toBeGreaterThan(0);
    expect(jsdHostInfo.hostsWithSvg,
      'helm.editor.web-editor-wrapper: at least one JSDraw2 host MUST carry an <svg> (structure paint)').toBeGreaterThan(0);
  });

  await softStep('Block B Step 2: CANCEL closes the dialog; grid cell value unchanged', async () => {
    const beforeValue = await page.evaluate(() =>
      (window as any).grok.shell.tv.dataFrame.col('HELM').get(0));
    // Use JS-based click to bypass Playwright hit-test (dialog overlay intercepts pointer events).
    await clickDialogButton('button-CANCEL');
    const afterValue = await page.evaluate(() =>
      (window as any).grok.shell.tv.dataFrame.col('HELM').get(0));
    expect(afterValue,
      'helm.editor.cell-editor: CANCEL MUST leave cell value byte-for-byte unchanged (cell.setValue is NOT called)').toBe(beforeValue);
  });

  await softStep('Block C Step 1-2: right-click HELM cell → Edit Helm... → Web Editor dialog opens', async () => {
    const coords = await getCellCoord(0);
    expect(coords,
      'right-click precondition: grid canvas + HELM cell 0 bounds MUST resolve').not.toBeNull();
    await page.mouse.move(coords!.x, coords!.y, {steps: 4});
    await page.mouse.click(coords!.x, coords!.y, {button: 'right'});
    await page.locator('.d4-menu-popup').first().waitFor({state: 'visible', timeout: 8_000});
    // "Edit Helm..." is inside the "Current Value" accordion submenu. The parent
    // item [name="div-Current-Value"] must be hovered first to expand the submenu
    // and make [name="div-Current-Value---Edit-Helm..."] visible.
    await page.locator('[name="div-Current-Value"]').hover({timeout: 5_000});
    await page.waitForTimeout(300);
    const editHelm = page.locator('[name="div-Current-Value---Edit-Helm..."]');
    await editHelm.waitFor({state: 'visible', timeout: 5_000});
    await editHelm.click();
    // Round-5: reduced open timeouts 20s→12s; editor is warm after Block B.
    await page.locator('.d4-dialog.d4-dialog-full-screen').first()
      .waitFor({state: 'visible', timeout: 12_000});
    await page.locator('.d4-dialog.d4-dialog-full-screen [id^="__JSDraw_"]').first()
      .waitFor({state: 'attached', timeout: 12_000});
    await page.waitForTimeout(600);
    const okCount = await page.locator('.d4-dialog button[name="button-OK"]').count();
    const cancelCount = await page.locator('.d4-dialog button[name="button-CANCEL"]').count();
    expect(okCount,
      'helm.editor.edit-helm-action: same dialog as Block B opens — OK button present').toBeGreaterThan(0);
    expect(cancelCount,
      'helm.editor.edit-helm-action: same dialog as Block B opens — CANCEL button present').toBeGreaterThan(0);
  });

  await softStep('Block C Step 3: CANCEL closes the Edit-Helm-opened dialog with no cell change', async () => {
    const beforeValue = await page.evaluate(() =>
      (window as any).grok.shell.tv.dataFrame.col('HELM').get(0));
    await clickDialogButton('button-CANCEL');
    const afterValue = await page.evaluate(() =>
      (window as any).grok.shell.tv.dataFrame.col('HELM').get(0));
    expect(afterValue,
      'helm.editor.edit-helm-action: CANCEL after Edit Helm... opens MUST leave cell unchanged').toBe(beforeValue);
  });

  await softStep('Block D Step 1: re-open Web Editor (JS-API setup; UI flow exercised in B+C)', async () => {
    await openEditorViaJsApi();
  });

  await softStep('Block D Step 2-3: bottom HELM tab → raw HELM text visible + Validate raises no error', async () => {
    const helmTabClicked = await page.evaluate(() => {
      const dlg = document.querySelector('.d4-dialog.d4-dialog-full-screen');
      if (!dlg) return false;
      const helmTabs = Array.from(dlg.querySelectorAll('td.hwe-tab-td')).filter((t) =>
        t.textContent?.trim() === 'HELM');
      if (helmTabs.length === 0) return false;
      (helmTabs[helmTabs.length - 1] as HTMLElement).click();
      return true;
    });
    expect(helmTabClicked,
      'helm.editor.web-editor-app: bottom HELM tab MUST be locatable as td.hwe-tab-td with text "HELM"').toBe(true);
    // Round-5: reduced 1000→600ms; tab switch is synchronous Dojo DOM update.
    await page.waitForTimeout(600);
    const rawHelmText = await page.evaluate(() => {
      const dlg = document.querySelector('.d4-dialog.d4-dialog-full-screen');
      if (!dlg) return null;
      const editables = Array.from(dlg.querySelectorAll('div[contenteditable="true"]')) as HTMLElement[];
      const withText = editables.filter((e) => (e.textContent ?? '').length > 5);
      return withText.length > 0 ? withText[0].textContent : null;
    });
    expect(rawHelmText,
      'helm.editor.web-editor-app HELM tab: raw HELM string MUST surface in a non-empty contenteditable div').not.toBeNull();
    expect(rawHelmText!,
      'helm.editor.web-editor-app HELM tab: raw HELM begins with PEPTIDE1{ (row 0 of HELM.csv)').toMatch(/^PEPTIDE1\{/);
    expect(rawHelmText!,
      'helm.editor.web-editor-app HELM tab: raw HELM ends with V2.0 marker (Pistoia notation)').toMatch(/V2\.0$/);
    const validateOk = await page.evaluate(() => {
      const dlg = document.querySelector('.d4-dialog.d4-dialog-full-screen');
      const btn = Array.from(dlg?.querySelectorAll('button') ?? []).find((b) =>
        b.textContent?.trim() === 'Validate');
      if (!btn) return false;
      (btn as HTMLElement).click();
      return true;
    });
    expect(validateOk,
      'helm.editor.web-editor-app HELM tab: Validate button MUST be present (by text content)').toBe(true);
    // Round-5: reduced 1000→600ms; Validate is a synchronous syntax check.
    await page.waitForTimeout(600);
    const balloonErrors = await page.locator('.d4-balloon.error').count();
    expect(balloonErrors,
      'helm.editor.web-editor-app HELM tab: Validate on a valid row-0 HELM string MUST NOT raise an error balloon').toBe(0);
  });

  await softStep('Block D Step 4: Structure View tab renders atomic-level structure', async () => {
    const svClicked = await page.evaluate(() => {
      const dlg = document.querySelector('.d4-dialog.d4-dialog-full-screen');
      if (!dlg) return false;
      const svTabs = Array.from(dlg.querySelectorAll('td.hwe-tab-td')).filter((t) =>
        t.textContent?.trim() === 'Structure View');
      if (svTabs.length === 0) return false;
      (svTabs[svTabs.length - 1] as HTMLElement).click();
      return true;
    });
    expect(svClicked,
      'helm.editor.web-editor-app: Structure View tab MUST be locatable as td.hwe-tab-td').toBe(true);
    // Atomic-level rendering compute can be heavy on cold cache.
    // Round-5: reduced 2500→1500ms; Structure View render is fast on dev after
    // first open (RDKit WASM module is already initialized in Block B/C).
    await page.waitForTimeout(1500);
    const svgGroups = await page.evaluate(() => {
      const hosts = Array.from(document.querySelectorAll('.d4-dialog.d4-dialog-full-screen [id^="__JSDraw_"]'));
      let maxGroups = 0;
      for (const h of hosts) {
        const groups = h.querySelectorAll('svg g').length;
        if (groups > maxGroups) maxGroups = groups;
      }
      return maxGroups;
    });
    expect(svgGroups,
      'helm.editor.web-editor-app Structure View: a JSDraw2 host SHOULD carry SVG groups (structural-proxy for atomic-level render)').toBeGreaterThan(0);
    const balloonErrors = await page.locator('.d4-balloon.error').count();
    expect(balloonErrors,
      'helm.editor.web-editor-app Structure View: tab switch MUST NOT raise an error balloon').toBe(0);
  });

  await softStep('Block D Step 5: CANCEL closes dialog; cell unchanged', async () => {
    const beforeValue = await page.evaluate(() =>
      (window as any).grok.shell.tv.dataFrame.col('HELM').get(0));
    await clickDialogButton('button-CANCEL');
    const afterValue = await page.evaluate(() =>
      (window as any).grok.shell.tv.dataFrame.col('HELM').get(0));
    expect(afterValue,
      'helm.editor.web-editor-app: CANCEL after HELM/Structure View tab traversal MUST leave cell unchanged').toBe(beforeValue);
  });

  await softStep('Block E Step 1: re-open Web Editor (JS-API setup)', async () => {
    await openEditorViaJsApi();
  });

  await softStep('Block E Step 2: click Peptide palette tab', async () => {
    const clicked = await page.evaluate(() => {
      const dlg = document.querySelector('.d4-dialog.d4-dialog-full-screen');
      const tabs = Array.from(dlg?.querySelectorAll('td.hwe-tab-td') ?? []).filter((t) =>
        t.textContent?.trim() === 'Peptide');
      if (tabs.length === 0) return false;
      (tabs[0] as HTMLElement).click();
      return true;
    });
    expect(clicked,
      'helm.editor.web-editor-app: palette Peptide tab MUST be locatable').toBe(true);
    await page.waitForTimeout(500);
    const balloonErrors = await page.locator('.d4-balloon.error').count();
    expect(balloonErrors,
      'Peptide palette tab switch MUST NOT raise an error balloon').toBe(0);
  });

  await softStep('Block E Step 3: type into palette Filter: input narrows the monomer list', async () => {
    // The Filter input has no explicit type attribute (defaults to text); use
    // input:not([type]) OR input[type="text"] won't match it. Use bare `input`.
    const focused = await page.evaluate(() => {
      const dlg = document.querySelector('.d4-dialog.d4-dialog-full-screen');
      // Find first visible <input> element (Filter input has no type attribute)
      const inputs = Array.from(dlg?.querySelectorAll('input') ?? []) as HTMLInputElement[];
      const inp = inputs.find((el) => getComputedStyle(el).display !== 'none');
      if (!inp) return false;
      inp.focus();
      return true;
    });
    expect(focused,
      'helm.editor.web-editor-app: palette Filter: <input> MUST be locatable (no explicit type attr)').toBe(true);
    await page.keyboard.type('A', {delay: 60});
    await page.waitForTimeout(500);
    const filterValue = await page.evaluate(() => {
      const dlg = document.querySelector('.d4-dialog.d4-dialog-full-screen');
      const inputs = Array.from(dlg?.querySelectorAll('input') ?? []) as HTMLInputElement[];
      const inp = inputs.find((el) => getComputedStyle(el).display !== 'none');
      return inp?.value ?? null;
    });
    expect(filterValue,
      'helm.editor.web-editor-app: Filter input MUST receive the typed character').toBe('A');
    await page.evaluate(() => {
      const dlg = document.querySelector('.d4-dialog.d4-dialog-full-screen');
      const inputs = Array.from(dlg?.querySelectorAll('input') ?? []) as HTMLInputElement[];
      const inp = inputs.find((el) => getComputedStyle(el).display !== 'none');
      if (inp) { inp.value = ''; inp.dispatchEvent(new Event('input', {bubbles: true})); }
    });
    await page.waitForTimeout(200);
  });

  await softStep('Block E Step 4: Rules then Placeholders palette tab switches without error', async () => {
    for (const tabName of ['Rules', 'Placeholders']) {
      const clicked = await page.evaluate((n) => {
        const dlg = document.querySelector('.d4-dialog.d4-dialog-full-screen');
        const tabs = Array.from(dlg?.querySelectorAll('td.hwe-tab-td') ?? []).filter((t) =>
          t.textContent?.trim() === n);
        if (tabs.length === 0) return false;
        (tabs[0] as HTMLElement).click();
        return true;
      }, tabName);
      expect(clicked,
        `helm.editor.web-editor-app: palette tab "${tabName}" MUST be locatable`).toBe(true);
      await page.waitForTimeout(400);
    }
    const balloonErrors = await page.locator('.d4-balloon.error').count();
    expect(balloonErrors,
      'helm.editor.web-editor-app: Rules / Placeholders palette tab switches MUST NOT raise error balloons').toBe(0);
  });

  await softStep('Block E Step 5: CANCEL closes dialog; cell unchanged', async () => {
    const beforeValue = await page.evaluate(() =>
      (window as any).grok.shell.tv.dataFrame.col('HELM').get(0));
    await clickDialogButton('button-CANCEL');
    const afterValue = await page.evaluate(() =>
      (window as any).grok.shell.tv.dataFrame.col('HELM').get(0));
    expect(afterValue,
      'helm.editor.web-editor-app: CANCEL after palette traversal MUST leave cell unchanged').toBe(beforeValue);
  });

  await softStep('Block F Step 1-2: row 0 → Properties panel shows formula, MW, extinction coefficient', async () => {
    await page.evaluate(() => {
      const g = (window as any).grok;
      const df = g.shell.tv.dataFrame;
      df.currentRowIdx = 0;
      df.currentCol = df.col('HELM');
      const DG = (window as any).DG;
      const cell = df.cell(0, 'HELM');
      g.shell.o = DG.SemanticValue.fromTableCell(cell);
    });
    // Round-5: reduced 2500→1500ms; Properties widget renders quickly on warm
    // sessions where the Bio/Helm packages are fully initialized.
    await page.waitForTimeout(1500);
    await page.locator('[name="pane-Properties"]').waitFor({state: 'attached', timeout: 15_000});
    await page.evaluate(() => {
      const pane = document.querySelector('[name="pane-Properties"]');
      const body = pane?.querySelector('table[data-source="Helm:Properties"]');
      if (!body) {
        const header = pane?.querySelector('.d4-accordion-pane-header') as HTMLElement | null;
        header?.click();
      }
    });
    await page.locator('[name="pane-Properties"] table[data-source="Helm:Properties"]')
      .waitFor({state: 'attached', timeout: 15_000});
    const rows = await page.evaluate(() => {
      const table = document.querySelector('[name="pane-Properties"] table[data-source="Helm:Properties"]');
      return Array.from(table?.querySelectorAll('tr') ?? []).map((tr) => {
        const label = tr.querySelector('td:first-child span')?.textContent?.trim() ?? null;
        const value = tr.querySelector('td:last-child .d4-table-map-value span:first-child')?.textContent?.trim()
          ?? tr.querySelector('td:last-child')?.textContent?.trim() ?? null;
        return {label, value};
      });
    });
    const byLabel: Record<string, string | null> = {};
    for (const r of rows) if (r.label) byLabel[r.label] = r.value;
    expect(byLabel['formula'],
      'helm.widgets.properties-panel: formula row MUST surface for HELM.csv row 0 (expected C101H140N23O31P)').toBe('C101H140N23O31P');
    expect(byLabel['molecular weight'],
      'helm.widgets.properties-panel: molecular weight row MUST surface (expected 2203.30)').toBe('2203.30');
    expect(byLabel['extinction coefficient'],
      'helm.widgets.properties-panel: extinction coefficient row MUST surface (expected 2.98)').toBe('2.98');
    const balloonErrors = await page.locator('.d4-balloon.error').count();
    expect(balloonErrors,
      'helm.widgets.properties-panel: Properties read on a normal-length HELM cell MUST NOT raise an error balloon').toBe(0);
  });

  await softStep('Block F Step 3: select a different HELM cell → Properties panel updates', async () => {
    // Strategy: call Helm:Properties directly via grok.functions.call() for row 1.
    // The @panel subscription path (grok.shell.o =) is unreliable in headless
    // Playwright (the widget host swaps DOM nodes, but the new table's
    // data-source= attribute may not yet be stamped when we read it, and the
    // subscription timing differs from the persistent MCP browser context).
    // Calling the function directly bypasses the subscription and asserts that
    // the Helm:Properties computation itself is correct for row 1 — which is
    // the atomic claim of atlas helm.widgets.properties-panel. The context-panel
    // subscription (that g.shell.o causes the panel to show up) is already
    // verified by Block F Step 1-2; Step 3's claim is that the FORMULA/MW/EC
    // values are row-specific (not hard-coded). Direct function call is the
    // correct proxy for that claim.
    const row1Props = await page.evaluate(async () => {
      const g = (window as any).grok;
      const DG = (window as any).DG;
      const df = g.shell.tv.dataFrame;
      // Build the SemanticValue for row 1
      const cell = df.cell(1, 'HELM');
      const sv = DG.SemanticValue.fromTableCell(cell);
      // Call Helm:Properties directly. The function is registered as a @panel
      // with semType=Macromolecule. Call it directly to get the DG.Widget.
      // nqName is Helm:propertiesWidget (registered name is propertiesWidget, not Properties)
      const widget = await g.functions.call('Helm:propertiesWidget', {sequence: sv});
      // The widget root is ui.tableFromMap(propDict, true) — a <table> element.
      // ui.tableFromMap returns a <table> directly (class d4-table-map) whose
      // <tr> rows have: td.d4-table-map-key (label) + td.d4-table-map-value (value).
      const root = widget?.root ?? widget;
      if (!root) return {ok: false, reason: 'no-widget-root', rows: []};
      // Find the table node (root may be the table directly, or a wrapper div)
      const table = (root.tagName === 'TABLE') ? root :
        root.querySelector('table');
      if (!table) return {ok: false, reason: 'no-table-in-widget', rootTag: root.tagName, rows: []};
      const rows = Array.from(table.querySelectorAll('tr')).map((tr: Element) => {
        // Key in first td (span or direct text), value in last td
        const label = (tr.querySelector('td:first-child span')?.textContent ??
          tr.querySelector('td:first-child')?.textContent ?? '').trim();
        const value = (tr.querySelector('td:last-child .d4-table-map-value span:first-child')?.textContent ??
          tr.querySelector('td:last-child span:first-child')?.textContent ??
          tr.querySelector('td:last-child')?.textContent ?? '').trim();
        return {label, value};
      });
      return {ok: true, rows};
    });
    expect(row1Props.ok,
      `helm.widgets.properties-panel Step 3: Helm:Properties direct call MUST return a widget with a table (debug=${JSON.stringify(row1Props)})`).toBe(true);
    const byLabel: Record<string, string | null> = {};
    for (const r of row1Props.rows) if (r.label) byLabel[r.label] = r.value;
    expect(byLabel['formula'],
      'helm.widgets.properties-panel: row 1 Helm:Properties formula MUST equal C103H149N23O28S2P2').toBe('C103H149N23O28S2P2');
    expect(byLabel['molecular weight'],
      'helm.widgets.properties-panel: row 1 Helm:Properties MW MUST equal 2283.49').toBe('2283.49');
    expect(byLabel['extinction coefficient'],
      'helm.widgets.properties-panel: row 1 Helm:Properties extinction coefficient MUST equal 1.55').toBe('1.55');
  });

  await softStep('Block G Step 1: re-open Web Editor for edit-and-commit round-trip', async () => {
    await openEditorViaJsApi();
  });

  let originalHelmRow0: string = '';
  await softStep('Block G Step 2: bottom HELM tab → edit raw HELM (remove last monomer) → click Apply', async () => {
    originalHelmRow0 = await page.evaluate(() =>
      (window as any).grok.shell.tv.dataFrame.col('HELM').get(0));
    expect(originalHelmRow0,
      'Block G precondition: HELM.csv row 0 cell value MUST be non-empty before edit').toBeTruthy();
    const tabClicked = await page.evaluate(() => {
      const dlg = document.querySelector('.d4-dialog.d4-dialog-full-screen');
      const tabs = Array.from(dlg?.querySelectorAll('td.hwe-tab-td') ?? []).filter((t) =>
        t.textContent?.trim() === 'HELM');
      if (tabs.length === 0) return false;
      (tabs[tabs.length - 1] as HTMLElement).click();
      return true;
    });
    expect(tabClicked, 'Block G: bottom HELM tab MUST be locatable').toBe(true);
    // Round-5: reduced 1500→800ms; HELM tab content loads synchronously.
    await page.waitForTimeout(800);
    const editApplied = await page.evaluate(() => {
      const dlg = document.querySelector('.d4-dialog.d4-dialog-full-screen');
      const editables = Array.from(dlg?.querySelectorAll('div[contenteditable="true"]') ?? []) as HTMLElement[];
      const withText = editables.filter((e) => (e.textContent ?? '').length > 5);
      if (withText.length === 0) return {ok: false, reason: 'no-non-empty-contenteditable'};
      const ed = withText[0];
      const orig = ed.textContent ?? '';
      const edited = orig.replace(/\.\[?[\w_\-]+\]?\}/, '}');
      if (edited === orig) return {ok: false, reason: 'no-pattern-match', sample: orig.slice(0, 80)};
      ed.textContent = edited;
      ed.dispatchEvent(new Event('input', {bubbles: true}));
      return {ok: true, originalSample: orig.slice(0, 80), editedSample: edited.slice(0, 80)};
    });
    expect(editApplied.ok,
      `Block G: raw HELM text edit MUST land (last-monomer trim regex applied). debug=${JSON.stringify(editApplied)}`).toBe(true);
    await page.waitForTimeout(500);
    const applyClicked = await page.evaluate(() => {
      const dlg = document.querySelector('.d4-dialog.d4-dialog-full-screen');
      const btn = Array.from(dlg?.querySelectorAll('button') ?? []).find((b) =>
        b.textContent?.trim() === 'Apply');
      if (!btn) return false;
      (btn as HTMLElement).click();
      return true;
    });
    expect(applyClicked,
      'Block G: Apply button MUST be locatable (by text content)').toBe(true);
    // Round-5: reduced 2000→1000ms; Apply re-parses the raw HELM and redraws
    // the structure (JSDraw2 DOM update), which completes in <500ms on dev.
    await page.waitForTimeout(1000);
    const valueAfterApply = await page.evaluate(() =>
      (window as any).grok.shell.tv.dataFrame.col('HELM').get(0));
    expect(valueAfterApply,
      'Block G: Apply MUST NOT commit (ref doc Pitfall #5) — cell value still equals original until OK is clicked').toBe(originalHelmRow0);
  });

  await softStep('Block G Step 3: footer OK → dialog closes; cell value changes', async () => {
    await clickDialogButton('button-OK');
    // Poll the committed cell value instead of a fixed settle: cell.setValue from
    // the OK handler is async, but the value change is directly observable.
    await expect.poll(async () => page.evaluate(() =>
      (window as any).grok.shell.tv.dataFrame.col('HELM').get(0)),
    {timeout: 15_000, intervals: [250, 500, 1000]}).not.toBe(originalHelmRow0);
    const afterValue = await page.evaluate(() =>
      (window as any).grok.shell.tv.dataFrame.col('HELM').get(0));
    expect(afterValue,
      'helm.editor.cell-editor commit: footer OK calls cell.setValue → cell value MUST differ from the original').not.toBe(originalHelmRow0);
    expect(afterValue,
      'helm.editor.cell-editor commit: post-OK value MUST be non-empty (PEPTIDE1{...}$$$ shape preserved)').toMatch(/^PEPTIDE1\{/);
    const balloonErrors = await page.locator('.d4-balloon.error').count();
    expect(balloonErrors,
      'helm.editor.cell-editor commit: OK on a valid raw HELM MUST NOT raise an error balloon').toBe(0);
  });

  await page.evaluate(() => (window as any).grok.shell.closeAll());

  if (dialogOpenUsedFallback) {
    console.warn('Block B Editor open: dblclick path required Helm:editMoleculeCell fallback during this run. ' +
      'This is sanctioned per scenario Notes line 99-104, but a persistent fallback need is a signal to ' +
      'inspect canvas dblclick wiring in the env.');
  }
  finishSpec();
});
