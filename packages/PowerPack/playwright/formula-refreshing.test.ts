/* ---
sub_features_covered: [powerpack.dialogs, powerpack.dialogs.add-new-column, powerpack.dialogs.add-new-column-func, powerpack.dialogs.prepare-add-column-call, powerpack.formula.is-formula-column, powerpack.formula.widget]
--- */
// GROK-17109: 3-step calc-column chain (Weight2→Weight3→Weight4) + Formula info panel edits +
// save+reopen persistence of formula tags.
//
// Load-bearing facts:
//   - The CM6 EditorView is exposed at `cmContent.cmTile.view` (NOT cmDiv.cmView.view, which is always
//     undefined and falls through to the brittle keyboard fallback). Same surface for dialog + panel widget.
//   - The Formula info panel widget renders `.add-new-column-dialog-cm-div` (NOT -widget-cm-div): the
//     constructor reads this.widget BEFORE assigning it, so the -dialog- class wins. Scope to .grok-prop-panel.
//   - Widget mode does NOT call prepareForSeleniumTests(), so button-Add-New-Column---OK is unset; the
//     widget's Apply button carries [name="button-Apply"] (text-content lookup as fallback).
//   - The Formula accordion-pane header click TOGGLES (not idempotent) and persists across column switches;
//     read .expanded/aria-expanded before clicking, else a blind click collapses it (cm-content hidden in DOM).
//   - The JS API path PowerPack:formulaWidget(col) is a deterministic fallback exercising the same widget code.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '@datagrok-libraries/test/src/playwright/spec-login';
import {finishSpec} from '@datagrok-libraries/test/src/playwright/viewers';
import {openTableFromFile, assertProvenanceScript} from '@datagrok-libraries/test/src/playwright/openers';
import {saveProjectWithProvenance, deleteProjectWithCleanup} from '@datagrok-libraries/test/src/playwright/projects';

test.use(specTestOptions);

test('PowerPack: Formula refreshing — 3-step calc-column chain + Formula info panel edits + GROK-17109 save+reopen persistence', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  const stamp = Date.now();
  const projectName = `AutoTest-FormulaRefreshing-${stamp}`;
  let projectId: string | null = null;
  let tableInfoId: string | null = null;

  await loginToDatagrok(page);

  await page.evaluate(() => {
    const grok = (window as any).grok;
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    try { grok.shell.closeAll(); } catch (_) {}
  });
  await page.waitForTimeout(500);

  try {
    // Setup: open demog.csv with datasync provenance so the GROK-17109 invariant can be tested.
    await softStep('Setup: open System:DemoFiles/demog.csv with datasync provenance', async () => {
      const opened = await openTableFromFile(page, 'System:DemoFiles/demog.csv');
      await page.locator('[name="viewer-Grid"]').waitFor({timeout: 60_000});
      await page.waitForTimeout(1000);
      // Without wired provenance, save-with-datasync degrades to snapshot-only and GROK-17109 can't be tested.
      await assertProvenanceScript(page, 'files', opened.script);
      const cols = await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        return df ? df.columns.names() : [];
      });
      expect(cols).toContain('WEIGHT');
    });

    // Block 1: Build a 3-step calc-column chain Weight2 → Weight3 → Weight4 via the Add-New-Column dialog.

    await softStep('B1.1: add Weight2 = ${WEIGHT}+100 via Add New Column dialog', async () => {
      await openAddNewColumnDialog(page);
      await composeAddNewColumn(page, 'Weight2', '${WEIGHT} + 100');
      await clickAddNewColumnOK(page);
      await waitForColumnPresent(page, 'Weight2');
      const check = await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        if (!df) return null;
        const w = df.col('WEIGHT'); const w2 = df.col('Weight2');
        for (let i = 0; i < Math.min(df.rowCount, 50); i++) {
          const wv = w.get(i); const w2v = w2.get(i);
          if (wv !== null && w2v !== null && Number.isFinite(wv) && Number.isFinite(w2v))
            return {wv, w2v, diff: w2v - wv};
        }
        return null;
      });
      expect(check).not.toBeNull();
      expect(check!.diff).toBeCloseTo(100, 1);
      const tag = await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        const w2 = df.col('Weight2');
        return w2?.tags?.get?.('formula') ?? w2?.tags?.get?.('.formula') ?? '';
      });
      expect(tag).toContain('${WEIGHT}');
    });

    await softStep('B1.2: add Weight3 = ${Weight2}+100 via Add New Column dialog', async () => {
      await openAddNewColumnDialog(page);
      await composeAddNewColumn(page, 'Weight3', '${Weight2} + 100');
      await clickAddNewColumnOK(page);
      await waitForColumnPresent(page, 'Weight3');
      const check = await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        if (!df) return null;
        const w = df.col('WEIGHT'); const w3 = df.col('Weight3');
        for (let i = 0; i < Math.min(df.rowCount, 50); i++) {
          const wv = w.get(i); const w3v = w3.get(i);
          if (wv !== null && w3v !== null && Number.isFinite(wv) && Number.isFinite(w3v))
            return {wv, w3v, diff: w3v - wv};
        }
        return null;
      });
      expect(check).not.toBeNull();
      expect(check!.diff).toBeCloseTo(200, 1);
    });

    await softStep('B1.3: add Weight4 = Log10(${Weight3})-0.2 via Add New Column dialog', async () => {
      await openAddNewColumnDialog(page);
      await composeAddNewColumn(page, 'Weight4', 'Log10(${Weight3}) - 0.2');
      await clickAddNewColumnOK(page);
      await waitForColumnPresent(page, 'Weight4');
      const check = await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        if (!df) return null;
        const w3 = df.col('Weight3'); const w4 = df.col('Weight4');
        for (let i = 0; i < Math.min(df.rowCount, 50); i++) {
          const w3v = w3.get(i); const w4v = w4.get(i);
          if (w3v !== null && w4v !== null && Number.isFinite(w3v) && Number.isFinite(w4v) && w3v > 0)
            return {w3v, w4v, expected: Math.log10(w3v) - 0.2};
        }
        return null;
      });
      expect(check).not.toBeNull();
      expect(check!.w4v).toBeCloseTo(check!.expected, 3);
    });

    // Block 2: Verify formula dependency recalc via the Formula info panel; downstream-only propagation.

    await softStep('B2.1+2: edit Weight2 formula via Formula info panel; verify Weight2/Weight3/Weight4 recompute', async () => {
      const pre = await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        return {
          w: df.col('WEIGHT').get(0),
          w2: df.col('Weight2').get(0),
          w3: df.col('Weight3').get(0),
          w4: df.col('Weight4').get(0),
        };
      });
      await editFormulaViaInfoPanel(page, 'Weight2', '${WEIGHT} + 200');
      await page.waitForTimeout(1000); // recalc settle.
      const post = await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        return {
          w: df.col('WEIGHT').get(0),
          w2: df.col('Weight2').get(0),
          w3: df.col('Weight3').get(0),
          w4: df.col('Weight4').get(0),
          w2Tag: df.col('Weight2')?.tags?.get?.('formula') ?? df.col('Weight2')?.tags?.get?.('.formula') ?? '',
        };
      });
      expect(post.w).toBeCloseTo(pre.w, 6); // source WEIGHT unaffected
      expect(post.w2).toBeCloseTo(post.w + 200, 1);
      expect(post.w2).not.toBeCloseTo(pre.w2, 1);
      // Weight3/Weight4 recompute automatically downstream.
      expect(post.w3).toBeCloseTo(post.w2 + 100, 1);
      expect(post.w3).not.toBeCloseTo(pre.w3, 1);
      expect(Number.isFinite(post.w4)).toBe(true);
      expect(post.w4).toBeCloseTo(Math.log10(post.w3) - 0.2, 3);
      expect(post.w2Tag).toContain('+ 200');
      const warnings = await page.evaluate(() => {
        try { return ((window as any).grok.shell.warnings || []).length; } catch { return 0; }
      });
      expect(warnings).toBe(0);
    });

    await softStep('B2.3: edit Weight3 formula via Formula info panel; Weight3/Weight4 recompute, Weight2 unaffected', async () => {
      const pre = await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        return {
          w2: df.col('Weight2').get(0),
          w3: df.col('Weight3').get(0),
          w4: df.col('Weight4').get(0),
        };
      });
      await editFormulaViaInfoPanel(page, 'Weight3', '${Weight2} + 50');
      await page.waitForTimeout(1000);
      const post = await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        return {
          w2: df.col('Weight2').get(0),
          w3: df.col('Weight3').get(0),
          w4: df.col('Weight4').get(0),
          w3Tag: df.col('Weight3')?.tags?.get?.('formula') ?? df.col('Weight3')?.tags?.get?.('.formula') ?? '',
        };
      });
      expect(post.w2).toBeCloseTo(pre.w2, 6); // Weight2 upstream — unaffected
      expect(post.w3).toBeCloseTo(post.w2 + 50, 1);
      expect(post.w3).not.toBeCloseTo(pre.w3, 1);
      expect(Number.isFinite(post.w4)).toBe(true);
      expect(post.w4).toBeCloseTo(Math.log10(post.w3) - 0.2, 3);
      expect(post.w3Tag).toContain('+ 50');
      const warnings = await page.evaluate(() => {
        try { return ((window as any).grok.shell.warnings || []).length; } catch { return 0; }
      });
      expect(warnings).toBe(0);
    });

    await softStep('B2.4: edit Weight4 formula via Formula info panel; Weight4 recomputes, Weight2/Weight3 unaffected', async () => {
      const pre = await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        return {
          w2: df.col('Weight2').get(0),
          w3: df.col('Weight3').get(0),
          w4: df.col('Weight4').get(0),
        };
      });
      await editFormulaViaInfoPanel(page, 'Weight4', 'Log10(${Weight3}) - 0.1');
      await page.waitForTimeout(1000);
      const post = await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        return {
          w2: df.col('Weight2').get(0),
          w3: df.col('Weight3').get(0),
          w4: df.col('Weight4').get(0),
          w4Tag: df.col('Weight4')?.tags?.get?.('formula') ?? df.col('Weight4')?.tags?.get?.('.formula') ?? '',
        };
      });
      // Terminal-node edit must NOT propagate upstream to Weight2/Weight3.
      expect(post.w2).toBeCloseTo(pre.w2, 6);
      expect(post.w3).toBeCloseTo(pre.w3, 6);
      expect(Number.isFinite(post.w4)).toBe(true);
      expect(post.w4).toBeCloseTo(Math.log10(post.w3) - 0.1, 3);
      expect(post.w4).not.toBeCloseTo(pre.w4, 3);
      expect(post.w4Tag).toContain('- 0.1');
      const warnings = await page.evaluate(() => {
        try { return ((window as any).grok.shell.warnings || []).length; } catch { return 0; }
      });
      expect(warnings).toBe(0);
    });

    // Block 3: Save, close, reopen — verify chained calc columns persist with formula tags (GROK-17109).

    let preSaveSnapshot: {w: number; w2: number; w3: number; w4: number} | null = null;

    await softStep('B3.1: save project with datasync provenance (chained calc columns + formula tags persisted)', async () => {
      preSaveSnapshot = await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        return {
          w: df.col('WEIGHT').get(0),
          w2: df.col('Weight2').get(0),
          w3: df.col('Weight3').get(0),
          w4: df.col('Weight4').get(0),
        };
      });
      const saved = await saveProjectWithProvenance(page, projectName);
      projectId = saved.projectId;
      tableInfoId = saved.tableInfoId;
      expect(projectId).toBeTruthy();
      const exists = await page.evaluate(async (pid) => {
        const grok = (window as any).grok;
        const p = await grok.dapi.projects.find(pid);
        return p != null;
      }, projectId);
      expect(exists).toBe(true);
    });

    await softStep('B3.2: close the project / working state cleared', async () => {
      await page.evaluate(() => {
        try { (window as any).grok.shell.closeAll(); } catch (_) {}
      });
      await page.waitForTimeout(1000);
      const tableCount = await page.evaluate(() => {
        try { return Number((window as any).grok.shell.tables?.length) || 0; }
        catch { return 0; }
      });
      expect(tableCount).toBe(0);
    });

    await softStep('B3.3+4: reopen project; verify Weight2/Weight3/Weight4 persist with formula tags + values intact (GROK-17109)', async () => {
      if (!projectId) throw new Error('B3.1 did not produce a projectId');
      if (!preSaveSnapshot) throw new Error('B3.1 did not capture pre-save snapshot');
      const reopen = await page.evaluate(async (pid) => {
        const grok = (window as any).grok;
        const p = await grok.dapi.projects.find(pid);
        await p.open();
        // File source re-executes OpenFile, slower than snapshot — wait for tables to re-materialize.
        for (let i = 0; i < 40; i++) {
          const tv = grok.shell.tv;
          if (tv?.dataFrame) break;
          await new Promise((r) => setTimeout(r, 500));
        }
        await new Promise((r) => setTimeout(r, 2000));
        const df = grok.shell.tv?.dataFrame;
        if (!df) return {ok: false, why: 'no df after reopen'};
        const names = df.columns.names();
        const w2 = df.col('Weight2');
        const w3 = df.col('Weight3');
        const w4 = df.col('Weight4');
        const tagOf = (c: any) => c?.tags?.get?.('formula') ?? c?.tags?.get?.('.formula') ?? '';
        return {
          ok: true,
          names,
          hasWeight2: names.includes('Weight2'),
          hasWeight3: names.includes('Weight3'),
          hasWeight4: names.includes('Weight4'),
          w2Formula: tagOf(w2),
          w3Formula: tagOf(w3),
          w4Formula: tagOf(w4),
          w: df.col('WEIGHT')?.get?.(0) ?? null,
          w2v: w2?.get?.(0) ?? null,
          w3v: w3?.get?.(0) ?? null,
          w4v: w4?.get?.(0) ?? null,
          warnings: (() => { try { return (grok.shell.warnings || []).length; } catch { return 0; } })(),
        };
      }, projectId);
      expect(reopen.ok).toBe(true);
      // GROK-17109: all three calc columns persist on reopen (pre-fix they disappeared).
      expect(reopen.hasWeight2).toBe(true);
      expect(reopen.hasWeight3).toBe(true);
      expect(reopen.hasWeight4).toBe(true);
      // Formula tags preserved (else the column reverts to a plain snapshot column).
      expect(reopen.w2Formula.length).toBeGreaterThan(0);
      expect(reopen.w3Formula.length).toBeGreaterThan(0);
      expect(reopen.w4Formula.length).toBeGreaterThan(0);
      // Tags carry the last-applied (block 2) formula state, not the original block-1 formulas.
      expect(reopen.w2Formula).toContain('+ 200');
      expect(reopen.w3Formula).toContain('+ 50');
      expect(reopen.w4Formula).toContain('- 0.1');
      // Values reflect the last-applied state; Weight2/3/4 recompute from post-reopen WEIGHT.
      expect(reopen.w).toBeCloseTo(preSaveSnapshot.w, 6);
      expect(reopen.w2v).toBeCloseTo(reopen.w + 200, 1);
      expect(reopen.w3v).toBeCloseTo(reopen.w2v + 50, 1);
      expect(reopen.w4v).toBeCloseTo(Math.log10(reopen.w3v) - 0.1, 3);
      expect(reopen.warnings).toBe(0);
    });
  } finally {
    await deleteProjectWithCleanup(page, {
      projectId: projectId ?? undefined,
      tableInfoId: tableInfoId ?? undefined,
    });
    await page.evaluate(() => {
      try {
        const cancel = document.querySelector('.d4-dialog [name="button-Add-New-Column---CANCEL"]') as HTMLElement | null;
        if (cancel) cancel.click();
        const anyCancel = document.querySelector('.d4-dialog [name="button-CANCEL"]') as HTMLElement | null;
        if (anyCancel) anyCancel.click();
      } catch (_) { /* best effort */ }
      try { (window as any).grok.shell.closeAll(); } catch (_) {}
    }).catch(() => {});
  }

  finishSpec();
});

// File-local helpers (kept inline; no powerpack helper module exists yet).

async function openAddNewColumnDialog(page: any): Promise<void> {
  await page.locator('.d4-dialog').first()
    .waitFor({state: 'detached', timeout: 5_000}).catch(() => {});
  const icon = page.locator('[name="icon-add-new-column"]').first();
  await icon.waitFor({timeout: 30_000, state: 'visible'});
  await icon.click({timeout: 10_000});
  const dlg = page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).first();
  await dlg.waitFor({timeout: 30_000});
  await expect(dlg).toBeVisible();
}

async function composeAddNewColumn(page: any, name: string, formula: string): Promise<void> {
  const dlg = page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).first();
  // Fill Name via native setter + input/change events (Dart InputBase listens on these).
  await page.evaluate((n: string) => {
    const input = document.querySelector('[name="input-Add-New-Column---Name"]') as HTMLInputElement | null;
    if (!input) throw new Error('Name input not found');
    const setter = Object.getOwnPropertyDescriptor(window.HTMLInputElement.prototype, 'value')!.set!;
    setter.call(input, n);
    input.dispatchEvent(new Event('input', {bubbles: true}));
    input.dispatchEvent(new Event('change', {bubbles: true}));
  }, name);
  await page.waitForTimeout(150);
  // Compose via CM6 view.dispatch; the EditorView is at cmContent.cmTile.view. force:true click bypasses overlay intercept.
  const cm = dlg.locator('.add-new-column-dialog-cm-div .cm-content').first();
  await cm.waitFor({timeout: 15_000, state: 'visible'});
  await cm.click({force: true});
  await page.waitForTimeout(200);
  let composed: {ok: boolean; doc?: string} = {ok: false};
  for (let i = 0; i < 10; i++) {
    composed = await page.evaluate((f: string) => {
      const cmContent = document.querySelector(
        '.d4-dialog .add-new-column-dialog-cm-div .cm-content') as HTMLElement | null;
      if (!cmContent) return {ok: false};
      const view = (cmContent as any).cmTile?.view ?? null;
      if (!view) return {ok: false};
      view.dispatch({changes: {from: 0, to: view.state.doc.length, insert: f}});
      return {ok: true, doc: view.state.doc.toString()};
    }, formula);
    if (composed.ok) break;
    await page.waitForTimeout(200);
  }
  // Keyboard fallback if cmTile.view never surfaced.
  if (!composed.ok) {
    await cm.click({force: true});
    await page.keyboard.press('Control+A');
    await page.keyboard.press('Delete');
    await page.waitForTimeout(100);
    await page.keyboard.type(formula, {delay: 30});
    await page.waitForTimeout(200);
    composed = await page.evaluate(() => {
      const cmContent = document.querySelector(
        '.d4-dialog .add-new-column-dialog-cm-div .cm-content') as HTMLElement | null;
      if (!cmContent) return {ok: false};
      const view = (cmContent as any).cmTile?.view ?? null;
      const doc = view ? view.state.doc.toString() : (cmContent.innerText || '');
      return {ok: true, doc};
    });
  }
  if (!composed.ok)
    throw new Error('composeAddNewColumn: CodeMirror cmTile.view not exposed even after keyboard fallback');
  // Contains-check (keyboard fallback may collapse whitespace); downstream value checks are load-bearing.
  const doc = composed.doc || '';
  const firstToken = formula.split(/\s+/)[0];
  if (firstToken)
    expect(doc).toContain(firstToken);
}

async function clickAddNewColumnOK(page: any): Promise<void> {
  const dlg = page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).first();
  await dlg.locator('[name="button-Add-New-Column---OK"]').first().click();
  // Calc columns finish evaluating asynchronously — see waitForColumnPresent for the arrival poll.
  await page.locator('.d4-dialog').first()
    .waitFor({state: 'detached', timeout: 10_000}).catch(() => {});
}

async function waitForColumnPresent(page: any, columnName: string): Promise<void> {
  let added = false;
  for (let i = 0; i < 40; i++) {
    added = await page.evaluate((n: string) => {
      const df = (window as any).grok.shell.tv?.dataFrame;
      return df ? df.columns.names().includes(n) : false;
    }, columnName);
    if (added) break;
    await page.waitForTimeout(250);
  }
  expect(added).toBe(true);
}

// Edit a calc-column formula via the Formula info panel widget (PowerPack:formulaWidget spawns the same
// AddNewColumnDialog pre-bound to the column). Accordion-pane-click path with a direct-JS-API fallback.
async function editFormulaViaInfoPanel(page: any, columnName: string, newFormula: string): Promise<void> {
  // Set grok.shell.o to the column (same end state as a header click, deterministic under headless).
  await page.evaluate((n: string) => {
    const grok = (window as any).grok;
    const df = grok.shell.tv?.dataFrame;
    if (!df) throw new Error('editFormulaViaInfoPanel: no active dataframe');
    const col = df.col(n);
    if (!col) throw new Error(`editFormulaViaInfoPanel: column "${n}" not found`);
    grok.shell.o = col;
  }, columnName);
  await page.waitForTimeout(500); // let Context Panel populate

  // Find the "Formula" pane header by text in .grok-prop-panel. Header click TOGGLES (not idempotent) and
  // expansion persists across column switches; read .expanded/aria-expanded first, click only when not expanded.
  const accordionPathWorked = await page.evaluate(() => {
    const propPanel = document.querySelector('.grok-prop-panel');
    if (!propPanel) return false;
    const headers = Array.from(propPanel.querySelectorAll('.d4-accordion-pane-header, .d4-accordion-title'));
    const target = headers.find((h) => (h.textContent || '').trim() === 'Formula') as HTMLElement | undefined;
    if (!target) return false;
    const isExpanded = target.classList.contains('expanded') ||
      target.getAttribute('aria-expanded') === 'true';
    if (!isExpanded) target.click();
    return true;
  });
  await page.waitForTimeout(500); // accordion expansion + widget re-render settle

  // Wait for the widget's CM host. It renders .add-new-column-dialog-cm-div (the constructor reads
  // this.widget before assigning it, so -dialog- wins); scope to .grok-prop-panel.
  let widgetCmFound = false;
  if (accordionPathWorked) {
    for (let i = 0; i < 25; i++) {
      widgetCmFound = await page.evaluate(() => {
        const propPanel = document.querySelector('.grok-prop-panel');
        return !!propPanel?.querySelector('.add-new-column-dialog-cm-div .cm-content');
      });
      if (widgetCmFound) break;
      await page.waitForTimeout(200);
    }
  }

  // Fallback: spawn the widget directly via JS API (exercises the same panel render function).
  if (!widgetCmFound) {
    await page.evaluate(async (n: string) => {
      const grok = (window as any).grok;
      const ui = (window as any).ui;
      const DG = (window as any).DG;
      const df = grok.shell.tv?.dataFrame;
      const col = df.col(n);
      const widget: any = await grok.functions.call('PowerPack:formulaWidget', {col});
      const propPanel = document.querySelector('.grok-prop-panel');
      const host = ui.div([widget.root]);
      host.style.padding = '8px';
      host.setAttribute('data-fr-fallback', '1'); // tag for cleanup
      propPanel?.appendChild(host);
      void DG; // silence unused
    }, columnName);
    for (let i = 0; i < 25; i++) {
      widgetCmFound = await page.evaluate(() => {
        const propPanel = document.querySelector('.grok-prop-panel');
        return !!propPanel?.querySelector('.add-new-column-dialog-cm-div .cm-content');
      });
      if (widgetCmFound) break;
      await page.waitForTimeout(200);
    }
  }

  if (!widgetCmFound)
    throw new Error(`editFormulaViaInfoPanel: Formula widget CM host (.grok-prop-panel .add-new-column-dialog-cm-div .cm-content) did not render for column "${columnName}"`);

  // Dispatch the new formula via cmTile.view (same as composeAddNewColumn). Self-heal if the pane collapsed.
  await page.evaluate(() => {
    const pp = document.querySelector('.grok-prop-panel');
    if (!pp) return;
    const cm = pp.querySelector('.add-new-column-dialog-cm-div .cm-content') as HTMLElement | null;
    if (cm && cm.offsetWidth === 0 && cm.offsetHeight === 0) {
      const headers = Array.from(pp.querySelectorAll('.d4-accordion-pane-header, .d4-accordion-title'));
      const formulaHeader = headers.find((h) => (h.textContent || '').trim() === 'Formula') as HTMLElement | undefined;
      if (formulaHeader) formulaHeader.click();
    }
  });
  await page.waitForTimeout(400);
  const panelCm = page.locator(
    '.grok-prop-panel .add-new-column-dialog-cm-div .cm-content').first();
  await panelCm.waitFor({timeout: 15_000, state: 'visible'});
  await panelCm.click({force: true});
  await page.waitForTimeout(200);
  let composed: {ok: boolean; doc?: string} = {ok: false};
  for (let i = 0; i < 10; i++) {
    composed = await page.evaluate((f: string) => {
      const propPanel = document.querySelector('.grok-prop-panel');
      const cmContent = propPanel?.querySelector(
        '.add-new-column-dialog-cm-div .cm-content') as HTMLElement | null;
      if (!cmContent) return {ok: false};
      const view = (cmContent as any).cmTile?.view ?? null;
      if (!view) return {ok: false};
      view.dispatch({changes: {from: 0, to: view.state.doc.length, insert: f}});
      return {ok: true, doc: view.state.doc.toString()};
    }, newFormula);
    if (composed.ok) break;
    await page.waitForTimeout(200);
  }
  // Keyboard fallback for the in-panel widget.
  if (!composed.ok) {
    await panelCm.click({force: true});
    await page.keyboard.press('Control+A');
    await page.keyboard.press('Delete');
    await page.waitForTimeout(100);
    await page.keyboard.type(newFormula, {delay: 30});
    await page.waitForTimeout(200);
    composed = await page.evaluate(() => {
      const propPanel = document.querySelector('.grok-prop-panel');
      const cmContent = propPanel?.querySelector(
        '.add-new-column-dialog-cm-div .cm-content') as HTMLElement | null;
      if (!cmContent) return {ok: false};
      const view = (cmContent as any).cmTile?.view ?? null;
      const doc = view ? view.state.doc.toString() : (cmContent.innerText || '');
      return {ok: true, doc};
    });
  }
  if (!composed.ok)
    throw new Error('editFormulaViaInfoPanel: CodeMirror cmTile.view not exposed on in-panel widget even after keyboard fallback');
  // Contains-check (keyboard fallback may collapse whitespace); recalc value assertions are load-bearing.
  const doc = composed.doc || '';
  const firstToken = newFormula.split(/\s+/)[0];
  if (firstToken)
    expect(doc).toContain(firstToken);

  // Widget mode has no OK button (prepareForSeleniumTests is skipped); click [name="button-Apply"], text fallback.
  let applyClicked = false;
  for (let i = 0; i < 25; i++) {
    applyClicked = await page.evaluate(() => {
      const propPanel = document.querySelector('.grok-prop-panel');
      if (!propPanel) return false;
      let apply = propPanel.querySelector('[name="button-Apply"]') as HTMLElement | null;
      if (!apply) {
        const buttons = Array.from(propPanel.querySelectorAll('button, .ui-btn, .d4-button'));
        apply = (buttons.find((b) => (b.textContent || '').trim() === 'Apply') as HTMLElement | undefined) || null;
      }
      if (!apply) return false;
      const disabled = (apply as HTMLButtonElement).disabled;
      if (disabled) return false;
      apply.click();
      return true;
    });
    if (applyClicked) break;
    await page.waitForTimeout(200);
  }
  if (!applyClicked)
    throw new Error('editFormulaViaInfoPanel: Apply button not found (or stayed disabled) in Formula widget panel');
  await page.waitForTimeout(400); // addNewColumnAction is async; caller waits a further 1s for recalc

  // Detach the fallback widget host if it was added.
  await page.evaluate(() => {
    document.querySelectorAll('[data-fr-fallback="1"]').forEach((el) => el.remove());
  });
}
