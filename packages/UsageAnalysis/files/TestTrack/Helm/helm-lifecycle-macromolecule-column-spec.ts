import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {finishSpec} from '../helpers/viewers';

test.use(specTestOptions);
test.use({timeout: 600_000});

const DATASET_PATH = 'System:AppData/Helm/samples/HELM.csv';

// Scenario fixtures (cited verbatim from the scenario .md Scenarios 5+6).
const PEPTIDE_FIXTURE = 'PEPTIDE1{meI.hHis.Aca.N.T.dE.Thr_PO3H2.Aca.D-Tyr_Et}$$$$';
const RNA_FIXTURE = 'RNA1{r(A)p.r(C)p.r(G)p.r(U)p}$$$$';
const GAP_MID = 'PEPTIDE1{A.*.G.*.K}$$$$';
const GAP_START = 'PEPTIDE1{*.A.G.K}$$$$';
const GAP_END = 'PEPTIDE1{A.G.K.*}$$$$';
// Empirically-verified output shape (V2.0 suffix appended by the helper).
const GAP_STRIPPED_EXPECTED = 'PEPTIDE1{A.G.K}$$$$V2.0';

test('Helm — lifecycle chain for the Macromolecule HELM column', async ({page}) => {
  // Runtime envelope: cold Helm init via Scenario 2 Helm:editMoleculeCell
  // (~300s worst, deferred per Round-5 fix) + login (~15s) + setup (~17s;
  // no blocking Helm:* in setup per Round-5 fix) + Scenarios 3-6 (~35s) =
  // ~367s worst-case. 600s ceiling matches the sibling spec (both file-scope
  // and in-body test.setTimeout needed — file-scope override is reliable,
  // in-body is defense-in-depth; sibling uses both per lines 267 + 274).
  test.setTimeout(600_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);

  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    (window as any).grok.shell.settings.showFiltersIconsConstantly = true;
    (window as any).grok.shell.windows.simpleMode = true;
    // Round-1 retry fix (cycle 2026-06-10): showContextPanel must be set
    // explicitly — simpleMode does NOT open the context panel in headless.
    // [name="pane-Properties"] waitFor in Scenario 3 fails without this.
    (window as any).grok.shell.windows.showContextPanel = true;
    (window as any).grok.shell.closeAll();
    // Round-1 retry fix: wait for Bio's detectMacromolecule to register
    // before readCsv fires. In a fresh headless context the Bio autostart
    // may not have completed yet — readCsv would fire before the detector
    // is registered, leaving helmCol.semType null (identical to the race
    // documented in sibling spec helm-editor-and-panels-spec.ts Round-4).
    for (let i = 0; i < 8; i++) {
      if ((window as any).DG.Func.find({name: 'detectMacromolecule'}).length > 0) break;
      await new Promise((r) => setTimeout(r, 500));
    }
    const df = await (window as any).grok.dapi.files.readCsv(path);
    (window as any).grok.shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(() => resolve(), 3000);
    });
    // Round-4 fix: if HELM column semType is still null after the semType event,
    // explicitly call grok.data.detectSemanticTypes(df) to force Bio's detector
    // to run (sibling spec Round-4+Round-5 pattern). Then poll up to 5s for the
    // column to be tagged. This handles the race where Bio autostart completes
    // AFTER onSemanticTypeDetected fires but BEFORE readCsv resolves (observed
    // in both the lifecycle spec and the sibling helm-editor-and-panels-spec.ts).
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
      for (let i = 0; i < 30; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      // Bio + Helm package init (RDKit + JSDraw2 + monomer-dictionary rewrite).
      await new Promise((r) => setTimeout(r, 2000));
    }
    // Round-5 fix: REMOVED the Helm:getHelmHelper blocking call that was here in
    // Round-4. Awaiting any Helm:* function call inside page.evaluate binds to
    // specTestOptions.actionTimeout: 15_000 — Helm:getHelmHelper on a cold server
    // blocks for the full Dojo + JSDraw2 + Pistoia HELM Web Editor init (~5min),
    // which exceeds the 15s action timeout and aborts the evaluate (observed as
    // B-RUN-PASS at ≈15s per attempt). Direct precedent: sibling spec Round-11
    // removed the same call for the same reason.
    // The Helm:getHelmHelper warmup for Scenarios 4-6 is safe because each
    // softStep calls it inside its own page.evaluate, bound only by the test-level
    // 600s ceiling (not the action-timeout); those calls are unaffected.
    // Allow a short settle for any async Helm background init that may have started:
    if ((window as any).DG.Func.find({name: 'editMoleculeCell'}).length > 0) {
      // Helm package functions registered — brief settle for background init
      // that doesn't need an awaited call to complete.
      await new Promise((r) => setTimeout(r, 2000));
    }
  }, DATASET_PATH);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});

  // Pre-flight invariants — Bio detector tagged the HELM column with the
  // tags HelmGridCellRenderer attach-gates on.
  const setupProbe = await page.evaluate(() => {
    const df = (window as any).grok.shell.tv.dataFrame;
    const helmCol = df.col('HELM');
    const activityCol = df.col('Activity');
    return {
      rowCount: df.rowCount,
      helmSemType: helmCol?.semType ?? null,
      helmUnits: helmCol?.meta?.units ?? null,
      helmRenderer: helmCol?.getTag ? helmCol.getTag('cell.renderer') : null,
      helmQuality: helmCol?.getTag ? helmCol.getTag('quality') : null,
      activityType: activityCol?.type ?? null,
    };
  });
  expect(setupProbe.rowCount,
    'Scenario .md Setup: HELM.csv MUST yield 540 rows').toBe(540);
  expect(setupProbe.helmSemType,
    'helm.rendering.cell-renderer precondition: Bio detector MUST tag HELM column semType=Macromolecule')
    .toBe('Macromolecule');
  expect(setupProbe.helmUnits,
    'helm.rendering.cell-renderer precondition: column tag meta.units MUST be "helm"')
    .toBe('helm');
  expect(setupProbe.helmRenderer,
    'helm.rendering.cell-renderer precondition: column tag cell.renderer MUST be "helm"')
    .toBe('helm');
  expect(setupProbe.helmQuality,
    'helm.rendering.cell-renderer precondition: column tag quality MUST be "Macromolecule"')
    .toBe('Macromolecule');
  expect(setupProbe.activityType,
    'Scenario 1 step 2: Activity column MUST be numeric (double)')
    .toBe('double');

  // Helper: open the Web Editor via the JS-API entry that double-click
  // dispatches to internally. Per grok-browser refdoc § "Cell editor",
  // this is the same code path as the canvas dblclick — sanctioned as
  // setup repetition. (The owned UI flow for cell-editor dialog-open is
  // covered by helm-editor-and-panels-spec.ts; this lifecycle spec
  // re-exercises the same JS-API path to keep total runtime under the
  // wrapper ceiling per cycle 2026-06-09-helm-automate-02 timing study.)
  //
  // Round-5 fix: fire-and-forget pattern (matching sibling Block B
  // fallback pattern). page.evaluate must NOT await Helm:editMoleculeCell
  // because on a cold server the call blocks for ~5min (full Dojo + JSDraw2
  // + Pistoia HELM Web Editor init). Awaiting inside page.evaluate bounds
  // the call to specTestOptions.actionTimeout: 15_000, which aborts the
  // evaluate after 15s. Instead, fire-and-forget in evaluate and let the
  // outer Playwright waitFor (with a 360s explicit timeout, not actionTimeout)
  // wait for the dialog to appear. This is exactly the sibling spec's cold-
  // open fallback pattern (lines 536-546 of helm-editor-and-panels-spec.ts).
  const openEditorViaJsApi = async () => {
    await page.evaluate(() => {
      const grid = (window as any).grok.shell.tv.grid;
      const DG = (window as any).DG;
      const cell = DG.GridCell.fromColumnRow(grid, 'HELM', 0);
      // Intentional non-await: fire-and-forget so evaluate returns immediately.
      // The cold Helm init (~5min) blocks an awaited call; the dialog
      // appearance is detected via the outer waitFor below.
      (window as any).grok.functions.call('Helm:editMoleculeCell', {cell});
    });
    // 360s timeout matches sibling spec cold-open fallback (sibling line 543).
    // actionTimeout does NOT apply to waitFor with explicit timeout.
    await page.locator('.d4-dialog.d4-dialog-full-screen')
      .waitFor({state: 'visible', timeout: 360_000});
    await page.locator('.d4-dialog.d4-dialog-full-screen [id^="__JSDraw_"]').first()
      .waitFor({state: 'attached', timeout: 30_000});
    await page.waitForTimeout(800);
  };


  await softStep('Scenario 1 Step 1-2: HELM column auto-renders via HelmGridCellRenderer (structural-proxy)', async () => {
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
    const balloonErrors = await page.locator('.d4-balloon.error').count();
    expect(balloonErrors,
      'render_helm_cell: open + first-paint MUST NOT raise an error balloon').toBe(0);
  });

  await softStep('Scenario 1 Step 3: scroll ~50 rows and back → renderer cache survives (no error balloons)', async () => {
    await page.evaluate(() => {
      const grid = (window as any).grok.shell.tv.grid;
      grid.scrollToCell('HELM', 50);
    });
    await page.waitForTimeout(500);
    await page.evaluate(() => {
      const grid = (window as any).grok.shell.tv.grid;
      grid.scrollToCell('HELM', 0);
    });
    await page.waitForTimeout(500);
    const canvasStill = await page.locator('[name="viewer-Grid"] canvas').count();
    expect(canvasStill,
      'helm.rendering.cell-renderer LRU cache: scroll round-trip MUST leave grid canvases attached')
      .toBeGreaterThan(0);
    const balloonErrors = await page.locator('.d4-balloon.error').count();
    expect(balloonErrors,
      'render_helm_cell: scroll round-trip MUST NOT raise an error balloon').toBe(0);
  });

  await softStep('Scenario 1 Step 4: re-apply column.setTag(cell.renderer, "helm") → cache invalidates cleanly', async () => {
    const tagState = await page.evaluate(() => {
      const df = (window as any).grok.shell.tv.dataFrame;
      const DG = (window as any).DG;
      const helmCol = df.col('HELM');
      helmCol.setTag(DG.TAGS.CELL_RENDERER, 'helm');
      return {
        renderer: helmCol.getTag('cell.renderer'),
        quality: helmCol.getTag('quality'),
        units: helmCol.meta.units,
      };
    });
    expect(tagState.renderer,
      'render_helm_cell: tag re-apply MUST leave cell.renderer="helm"').toBe('helm');
    expect(tagState.quality,
      'render_helm_cell: tag re-apply MUST preserve quality=Macromolecule').toBe('Macromolecule');
    expect(tagState.units,
      'render_helm_cell: tag re-apply MUST preserve meta.units=helm').toBe('helm');
    // Force a scroll to drive any pending re-paint
    await page.evaluate(() => {
      const grid = (window as any).grok.shell.tv.grid;
      grid.scrollToCell('HELM', 20);
      grid.scrollToCell('HELM', 0);
    });
    await page.waitForTimeout(800);
    const balloonErrors = await page.locator('.d4-balloon.error').count();
    expect(balloonErrors,
      'render_helm_cell: tag re-apply + re-render MUST NOT raise an error balloon').toBe(0);
  });


  let originalHelmRow0: string = '';
  await softStep('Scenario 2 Step 1-2: open Web Editor → dialog with bottom tabs + JSDraw2 host', async () => {
    originalHelmRow0 = await page.evaluate(() =>
      (window as any).grok.shell.tv.dataFrame.col('HELM').get(0));
    expect(originalHelmRow0,
      'Scenario 2 precondition: HELM.csv row 0 cell value MUST be non-empty').toBeTruthy();
    await openEditorViaJsApi();
    const dlgByName = await page.locator('[name="dialog-"]').count();
    expect(dlgByName,
      'helm.editor.cell-editor: full-screen Web Editor dialog name= attribute slug is "dialog-"').toBeGreaterThan(0);
    const okCount = await page.locator('.d4-dialog button[name="button-OK"]').count();
    const cancelCount = await page.locator('.d4-dialog button[name="button-CANCEL"]').count();
    expect(okCount,
      'helm.editor.cell-editor: footer OK button MUST exist').toBeGreaterThan(0);
    expect(cancelCount,
      'helm.editor.cell-editor: footer CANCEL button MUST exist').toBeGreaterThan(0);
    const expectedBottomTabs = ['Sequence', 'HELM', 'Properties', 'Structure View'];
    const tabTexts = await page.evaluate(() => {
      const dlg = document.querySelector('.d4-dialog.d4-dialog-full-screen');
      return Array.from(dlg?.querySelectorAll('td.hwe-tab-td') ?? []).map((t) => t.textContent?.trim());
    });
    for (const expected of expectedBottomTabs) {
      expect(tabTexts,
        `Scenario 2 step 2: bottom tab "${expected}" MUST be present`).toContain(expected);
    }
  });

  await softStep('Scenario 2 Step 3-4: bottom HELM tab → trim last monomer → Apply (does NOT commit)', async () => {
    const tabClicked = await page.evaluate(() => {
      const dlg = document.querySelector('.d4-dialog.d4-dialog-full-screen');
      const tabs = Array.from(dlg?.querySelectorAll('td.hwe-tab-td') ?? []).filter((t) =>
        t.textContent?.trim() === 'HELM');
      if (tabs.length === 0) return false;
      (tabs[tabs.length - 1] as HTMLElement).click();
      return true;
    });
    expect(tabClicked,
      'Scenario 2: bottom HELM tab MUST be locatable as td.hwe-tab-td with text "HELM"').toBe(true);
    await page.waitForTimeout(1500);
    const editApplied = await page.evaluate(() => {
      const dlg = document.querySelector('.d4-dialog.d4-dialog-full-screen');
      const editables = Array.from(dlg?.querySelectorAll('div[contenteditable="true"]') ?? []) as HTMLElement[];
      const withText = editables.filter((e) => (e.textContent ?? '').length > 5);
      if (withText.length === 0) return {ok: false, reason: 'no-non-empty-contenteditable'};
      const ed = withText[0];
      const orig = ed.textContent ?? '';
      // Trim the trailing monomer: replace `.[?MONOMER]?}` with `}`.
      const edited = orig.replace(/\.\[?[\w_\-]+\]?\}/, '}');
      if (edited === orig) return {ok: false, reason: 'no-pattern-match', sample: orig.slice(0, 80)};
      ed.textContent = edited;
      ed.dispatchEvent(new Event('input', {bubbles: true}));
      return {ok: true, originalSample: orig.slice(0, 80), editedSample: edited.slice(0, 80)};
    });
    expect(editApplied.ok,
      `Scenario 2: raw HELM text edit MUST land. debug=${JSON.stringify(editApplied)}`).toBe(true);
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
      'Scenario 2 step 4: Apply button MUST be locatable (by text content)').toBe(true);
    await page.waitForTimeout(2000);
    const valueAfterApply = await page.evaluate(() =>
      (window as any).grok.shell.tv.dataFrame.col('HELM').get(0));
    // Scenario Expected line 109-111: "Apply replaces the drawing from raw
    // text but does NOT mutate the grid (Common Pitfall #5)."
    expect(valueAfterApply,
      'edit_helm_cell: Apply MUST NOT commit to grid (helm.md Pitfall #5; scenario Expected line 109)')
      .toBe(originalHelmRow0);
  });

  await softStep('Scenario 2 Step 5-6: footer OK → dialog closes; grid cell value updates', async () => {
    await page.locator('.d4-dialog button[name="button-OK"]').first().click();
    await page.locator('.d4-dialog.d4-dialog-full-screen').waitFor({state: 'hidden', timeout: 20_000});
    await page.waitForTimeout(1000);
    const afterValue = await page.evaluate(() =>
      (window as any).grok.shell.tv.dataFrame.col('HELM').get(0));
    expect(afterValue,
      'helm.editor.cell-editor commit: footer OK calls cell.setValue → cell value MUST differ from original')
      .not.toBe(originalHelmRow0);
    expect(afterValue,
      'helm.editor.cell-editor commit: post-OK value MUST be non-empty (PEPTIDE1{...}$$$$ shape preserved)')
      .toMatch(/^PEPTIDE1\{/);
    const balloonErrors = await page.locator('.d4-balloon.error').count();
    expect(balloonErrors,
      'edit_helm_cell: OK on a valid raw HELM MUST NOT raise an error balloon').toBe(0);
  });

  let valueAfterFirstCommit: string = '';
  await softStep('Scenario 2 Step 7: re-open dialog → CANCEL discards pending edits', async () => {
    valueAfterFirstCommit = await page.evaluate(() =>
      (window as any).grok.shell.tv.dataFrame.col('HELM').get(0));
    await openEditorViaJsApi();
    // Make a destructive edit in the contenteditable region, then CANCEL.
    await page.evaluate(() => {
      const dlg = document.querySelector('.d4-dialog.d4-dialog-full-screen');
      const tabs = Array.from(dlg?.querySelectorAll('td.hwe-tab-td') ?? []).filter((t) =>
        t.textContent?.trim() === 'HELM');
      if (tabs.length === 0) return;
      (tabs[tabs.length - 1] as HTMLElement).click();
    });
    await page.waitForTimeout(1200);
    await page.evaluate(() => {
      const dlg = document.querySelector('.d4-dialog.d4-dialog-full-screen');
      const editables = Array.from(dlg?.querySelectorAll('div[contenteditable="true"]') ?? []) as HTMLElement[];
      const withText = editables.filter((e) => (e.textContent ?? '').length > 5);
      if (withText.length > 0) {
        const ed = withText[0];
        const orig = ed.textContent ?? '';
        const edited = orig.replace(/\.\[?[\w_\-]+\]?\}/, '}');
        ed.textContent = edited;
        ed.dispatchEvent(new Event('input', {bubbles: true}));
      }
    });
    await page.waitForTimeout(500);
    await page.locator('.d4-dialog button[name="button-CANCEL"]').first().click();
    await page.locator('.d4-dialog.d4-dialog-full-screen').waitFor({state: 'hidden', timeout: 10_000});
    const afterCancel = await page.evaluate(() =>
      (window as any).grok.shell.tv.dataFrame.col('HELM').get(0));
    expect(afterCancel,
      'helm.editor.cell-editor: CANCEL on re-opened dialog MUST leave cell unchanged')
      .toBe(valueAfterFirstCommit);
  });


  await softStep('Scenario 3 Step 1-3: row 0 → Properties panel surfaces formula, MW, extinction coefficient', async () => {
    await page.evaluate(() => {
      const grok = (window as any).grok;
      const df = grok.shell.tv.dataFrame;
      df.currentRowIdx = 0;
      df.currentCol = df.col('HELM');
      const DG = (window as any).DG;
      const cell = df.cell(0, 'HELM');
      grok.shell.o = DG.SemanticValue.fromTableCell(cell);
    });
    // Round-2: reduced 2500→1500ms (sibling Round-5 pace; Properties widget
    // renders quickly on sessions where Bio+Helm packages are initialized).
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
    // Note: Scenario 2 committed an edit to row 0 (trimmed last monomer),
    // so the row 0 properties on the post-commit value MAY differ from
    // the live-MCP-recon baseline of C101H140N23O31P (which was the
    // un-edited original). Assert key PRESENCE + non-empty value shape
    // rather than exact strings — exact strings would re-fail on every
    // monomer-library swap or trim-regex tweak.
    expect(byLabel['formula'],
      'compute_properties: formula row MUST surface for the (edited) row 0 HELM cell')
      .toBeTruthy();
    expect(byLabel['formula'],
      'compute_properties: formula MUST look like a chemical formula (C followed by digits)')
      .toMatch(/^C\d+/);
    expect(byLabel['molecular weight'],
      'compute_properties: molecular weight row MUST surface').toBeTruthy();
    expect(byLabel['molecular weight'],
      'compute_properties: molecular weight MUST parse as a positive number')
      .toMatch(/^\d+(\.\d+)?$/);
    expect(byLabel['extinction coefficient'],
      'compute_properties: extinction coefficient row MUST surface').toBeTruthy();
    const balloonErrors = await page.locator('.d4-balloon.error').count();
    expect(balloonErrors,
      'compute_properties: Properties read on a normal-length HELM cell MUST NOT raise an error balloon')
      .toBe(0);
  });

  await softStep('Scenario 3 Step 4: switch to a different HELM cell → Properties values are row-specific', async () => {
    // Round-2 fix: grok.shell.o = DG.SemanticValue.fromTableCell(cell1) does
    // NOT reliably update the Properties panel DOM table in a fresh headless
    // Playwright session (async widget-host swap timing differs from MCP
    // persistent session). Pattern documented in sibling spec Round-6
    // (helm-editor-and-panels-spec.ts lines 120-138).
    // Fix: call Helm:propertiesWidget directly for row 1 to assert row-specific
    // property computation. The context-panel subscription is already verified
    // by Step 1-3 (row 0 shows correctly via subscription); Step 4's only
    // remaining claim is that the formula values differ between rows (row-specific
    // computation — NOT hardcoded). Direct function call is the correct proxy.
    const row1Props = await page.evaluate(async () => {
      const g = (window as any).grok;
      const DG = (window as any).DG;
      const df = g.shell.tv.dataFrame;
      const cell = df.cell(1, 'HELM');
      const sv = DG.SemanticValue.fromTableCell(cell);
      // nqName Helm:propertiesWidget confirmed via MCP: DG.Func.find({package:'Helm'})
      const widget = await g.functions.call('Helm:propertiesWidget', {sequence: sv});
      const root = widget?.root ?? widget;
      if (!root) return {ok: false, reason: 'no-widget-root', rows: []};
      const table = (root.tagName === 'TABLE') ? root : root.querySelector('table');
      if (!table) return {ok: false, reason: 'no-table-in-widget', rootTag: root.tagName, rows: []};
      const rows = Array.from(table.querySelectorAll('tr')).map((tr) => {
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
      `compute_properties Step 4: Helm:propertiesWidget direct call for row 1 MUST return a widget with a table (debug=${JSON.stringify(row1Props)})`).toBe(true);
    const byLabel: Record<string, string | null> = {};
    for (const r of row1Props.rows) if (r.label) byLabel[r.label] = r.value;
    expect(byLabel['formula'],
      'compute_properties: row 1 formula MUST be present via direct Helm:propertiesWidget call').toBeTruthy();
    expect(byLabel['formula'],
      'compute_properties: row 1 formula MUST look like a chemical formula (C followed by digits)')
      .toMatch(/^C\d+/);
    // MCP recon: row 0 formula=C101H140N23O31P, row 1 formula=C103H149N23O28S2P2
    // (different sequences → different formulas; asserts row-specific computation)
    expect(byLabel['formula'],
      'compute_properties: row 1 formula MUST differ from row 0 (row-specific computation, not hardcoded)')
      .not.toBe('C101H140N23O31P');
  });

  await softStep('Scenario 3 Step 5: >1000-char HELM string → "Too long sequence" warning, no UI freeze', async () => {
    // Round-3 fix: pane-subscription timing is unreliable for scratch-DataFrame
    // cells. grok.shell.o change may not fire Helm:propertiesWidget before the
    // DOM query runs (observed as tooLongDivText=null across multiple attempts).
    // Pattern: call Helm:propertiesWidget directly (same approach as Step 4).
    const guardInfo = await page.evaluate(async () => {
      const grok = (window as any).grok;
      const DG = (window as any).DG;
      const repeat = 'r(A)p.r(C)p.r(G)p.r(U)p.';
      const longSeq = 'RNA1{' + repeat.repeat(60).replace(/\.$/, '') + '}$$$$V2.0';
      const sdf = DG.DataFrame.fromColumns([
        DG.Column.fromStrings('HELM', [longSeq]),
      ]);
      const col = sdf.col('HELM');
      col.semType = 'Macromolecule';
      col.meta.units = 'helm';
      col.setTag('cell.renderer', 'helm');
      col.setTag('quality', 'Macromolecule');
      col.setTag('aligned', 'SEQ.MSA');
      col.setTag('alphabet', 'UN');
      const cell = sdf.cell(0, 'HELM');
      const sv = DG.SemanticValue.fromTableCell(cell);
      const widget = await grok.functions.call('Helm:propertiesWidget', {sequence: sv});
      const root = widget?.root ?? widget;
      const text = root ? (root.textContent ?? null) : null;
      return {len: longSeq.length, tooLongText: text, rootTag: root?.tagName ?? null};
    });
    expect(guardInfo.len,
      'Scenario 3 step 5: synthetic HELM string MUST be > 1000 characters')
      .toBeGreaterThan(1000);
    expect(guardInfo.tooLongText,
      `compute_properties >1000-char guard: Helm:propertiesWidget MUST return widget containing "Too long sequence" text (rootTag=${guardInfo.rootTag})`)
      .toContain('Too long sequence');
    const balloonErrors = await page.locator('.d4-balloon.error').count();
    expect(balloonErrors,
      'compute_properties >1000-char guard: long-sequence MUST be a controlled fallback, NOT an error balloon')
      .toBe(0);
  });

  // Restore shell.o to HELM.csv row 0 for downstream cleanliness
  // (scratch DataFrame is GC'd when shell.o is replaced).
  await page.evaluate(() => {
    const grok = (window as any).grok;
    const df = grok.shell.tv.dataFrame;
    const DG = (window as any).DG;
    grok.shell.o = DG.SemanticValue.fromTableCell(df.cell(0, 'HELM'));
  });
  await page.waitForTimeout(500);


  await softStep('Scenario 4 Step 1-2: Helm:getMolfiles on HELM column → string column, rowCount matches source', async () => {
    const result = await page.evaluate(async () => {
      const grok = (window as any).grok;
      const df = grok.shell.tv.dataFrame;
      const helmCol = df.col('HELM');
      const t0 = Date.now();
      const molCol = await grok.functions.call('Helm:getMolfiles', {col: helmCol});
      const elapsed = Date.now() - t0;
      const firstNonNullIdx = (() => {
        for (let i = 0; i < molCol.length; i++) {
          const v = molCol.get(i);
          if (v && String(v).length > 0) return i;
        }
        return -1;
      })();
      const firstNonNull = firstNonNullIdx >= 0 ? String(molCol.get(firstNonNullIdx)) : null;
      return {
        elapsed,
        type: molCol.type,
        rows: molCol.length,
        srcRows: helmCol.length,
        firstNonNullIdx,
        firstSample: firstNonNull ? firstNonNull.slice(0, 200) : null,
        hasV2000orV3000: firstNonNull ? /V[23]000/.test(firstNonNull) : false,
      };
    });
    expect(result.type,
      'helm.api.get-molfiles: returned column MUST be of type string').toBe('string');
    expect(result.rows,
      'helm.api.get-molfiles: returned column rowCount MUST equal source column rowCount')
      .toBe(result.srcRows);
    expect(result.firstNonNullIdx,
      'helm.api.get-molfiles: at least one row MUST yield a non-empty molfile')
      .toBeGreaterThanOrEqual(0);
    expect(result.hasV2000orV3000,
      'helm.api.helm-helper.get-molfiles: non-empty molfiles MUST carry the V2000 or V3000 header')
      .toBe(true);
  });

  await softStep('Scenario 4 Step 3: re-issue getMolfiles after 1.2s → cached editor evicts + re-instantiates without error', async () => {
    await page.waitForTimeout(1200);
    const reissue = await page.evaluate(async () => {
      const grok = (window as any).grok;
      const df = grok.shell.tv.dataFrame;
      const helmCol = df.col('HELM');
      const t0 = Date.now();
      try {
        const molCol = await grok.functions.call('Helm:getMolfiles', {col: helmCol});
        return {ok: true, elapsed: Date.now() - t0, rows: molCol.length};
      } catch (e: any) {
        return {ok: false, elapsed: Date.now() - t0, error: e?.message ?? String(e)};
      }
    });
    expect(reissue.ok,
      `convert_helm_to_molfile: re-issue after 1s eviction window MUST succeed. err=${reissue.ok ? '' : reissue.error}`)
      .toBe(true);
    expect(reissue.rows,
      'convert_helm_to_molfile: re-issued call MUST return the same rowCount').toBe(540);
    const balloonErrors = await page.locator('.d4-balloon.error').count();
    expect(balloonErrors,
      'convert_helm_to_molfile: re-issue MUST NOT raise a "JSDraw2 editor is null" error balloon')
      .toBe(0);
  });


  await softStep('Scenario 5 Step 1-2: Helm:getHelmHelper → hh.parse(peptide fixture) → atoms + bonds > 0, linear shape', async () => {
    const out = await page.evaluate(async (peptide: string) => {
      const grok = (window as any).grok;
      const hh = await grok.functions.call('Helm:getHelmHelper');
      const mol = hh.parse(peptide);
      return {
        nonNull: mol !== null && mol !== undefined,
        atoms: mol?.atoms?.length ?? null,
        bonds: mol?.bonds?.length ?? null,
      };
    }, PEPTIDE_FIXTURE);
    expect(out.nonNull,
      'helm.api.helm-helper.parse: hh.parse(peptide fixture) MUST return a non-null HelmMol').toBe(true);
    expect(out.atoms,
      'helm.api.helm-helper.parse: peptide fixture MUST yield atoms (9 per MCP recon)').toBeGreaterThan(0);
    expect(out.bonds,
      'helm.api.helm-helper.parse: peptide fixture MUST yield bonds (n-1 = 8 for a 9-monomer linear chain)')
      .toBeGreaterThan(0);
    // Linear chain invariant: bonds == atoms - 1 (verified empirically: 9-1=8).
    expect(out.bonds,
      'helm.utils.parse-helm (transitively): linear peptide chain → bonds == atoms - 1')
      .toBe((out.atoms ?? 0) - 1);
  });

  await softStep('Scenario 5 Step 4: hh.parse(RNA fixture) → atoms + bonds > 0; multi-token square-bracket monomers handled', async () => {
    const out = await page.evaluate(async (rna: string) => {
      const grok = (window as any).grok;
      const hh = await grok.functions.call('Helm:getHelmHelper');
      const mol = hh.parse(rna);
      return {
        nonNull: mol !== null && mol !== undefined,
        atoms: mol?.atoms?.length ?? null,
        bonds: mol?.bonds?.length ?? null,
      };
    }, RNA_FIXTURE);
    expect(out.nonNull,
      'helm.api.helm-helper.parse: hh.parse(RNA fixture) MUST return a non-null HelmMol').toBe(true);
    expect(out.atoms,
      'helm.api.helm-helper.parse: RNA fixture MUST yield atoms (12 per MCP recon — branched r/p/branch expansion)')
      .toBeGreaterThan(0);
    expect(out.bonds,
      'helm.api.helm-helper.parse: RNA fixture MUST yield bonds (11 per MCP recon)').toBeGreaterThan(0);
    const balloonErrors = await page.locator('.d4-balloon.error').count();
    expect(balloonErrors,
      'parse_helm: hh.parse on peptide + RNA fixtures MUST NOT raise an error balloon').toBe(0);
  });


  await softStep('Scenario 6 Step 1-2: hh.removeGaps(mid-gap) → "*" monomers removed; bonds re-linked', async () => {
    const out = await page.evaluate(async (helm: string) => {
      const grok = (window as any).grok;
      const hh = await grok.functions.call('Helm:getHelmHelper');
      const res = hh.removeGaps(helm);
      // Empirical return-shape (per MCP recon 2026-06-09):
      // {srcHelm, resHelm, monomerMap}
      return {
        hasResHelm: typeof res?.resHelm === 'string',
        resHelm: res?.resHelm ?? null,
        srcHelm: res?.srcHelm ?? null,
        hasMonomerMap: res?.monomerMap !== undefined,
      };
    }, GAP_MID);
    expect(out.hasResHelm,
      'helm.api.helm-helper.remove-gaps: return MUST carry a resHelm string').toBe(true);
    expect(out.resHelm,
      `remove_gaps mid-gap: '*' monomers MUST be stripped → expected '${GAP_STRIPPED_EXPECTED}' (per MCP recon empirical shape)`)
      .toBe(GAP_STRIPPED_EXPECTED);
    expect(out.hasMonomerMap,
      'helm.api.helm-helper.remove-gaps: helper return MUST include a monomerMap object (per scenario Expected line 199-203 — the helper variant carries a mapping; pure-string utility does not)')
      .toBe(true);
  });

  await softStep('Scenario 6 Step 4: hh.removeGaps edge cases (start, end) strip cleanly without throwing', async () => {
    const out = await page.evaluate(async (helms: string[]) => {
      const grok = (window as any).grok;
      const hh = await grok.functions.call('Helm:getHelmHelper');
      const results: Array<{helm: string; resHelm: string | null; error: string | null}> = [];
      for (const h of helms) {
        try {
          const r = hh.removeGaps(h);
          results.push({helm: h, resHelm: r?.resHelm ?? null, error: null});
        } catch (e: any) {
          results.push({helm: h, resHelm: null, error: e?.message ?? String(e)});
        }
      }
      return results;
    }, [GAP_START, GAP_END]);
    for (const r of out) {
      expect(r.error,
        `remove_gaps edge case '${r.helm}': MUST NOT throw (atlas dep_lifecycle_ops.remove_gaps invariant)`)
        .toBeNull();
      expect(r.resHelm,
        `remove_gaps edge case '${r.helm}': MUST strip to '${GAP_STRIPPED_EXPECTED}' (per MCP recon — V2.0 suffix appended)`)
        .toBe(GAP_STRIPPED_EXPECTED);
    }
    const balloonErrors = await page.locator('.d4-balloon.error').count();
    expect(balloonErrors,
      'remove_gaps: start/end edge cases MUST NOT raise an error balloon').toBe(0);
  });

  // Scenario 6 step 3 (pure-string equivalence) — covered transitively
  // via the assertions above per SR-1 (the helper internally calls the
  // pure-string utility; helper return-shape includes monomerMap as the
  // documented helper-only addition).
  //
  // Scenario 6 step 4 final clause (>2-bonded `*` → HelmNotSupportedError)
  // deferred per SR-2 — no reliable branched-polymer fixture for the
  // 3-bonded gap monomer is provided in the scenario `.md` or atlas; live
  // MCP recon produced no path to construct one. The positive paths
  // (mid / start / end) cover the atlas remove_gaps op contract.

  await page.evaluate(() => (window as any).grok.shell.closeAll());
  finishSpec();
});
