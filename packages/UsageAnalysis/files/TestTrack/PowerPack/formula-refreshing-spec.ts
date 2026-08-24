import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
import {openTableFromFile, assertProvenanceScript} from '../helpers/openers';
import {saveProjectWithProvenance, deleteProjectWithCleanup} from '../helpers/projects';
import {WARNING_BALLOON, proveBalloonChannel} from '../helpers/balloons';

test.use(specTestOptions);

// Reads the absence claim first — a balloon raised by the action under test is still on
// screen at this point — then probes, so that "no warning balloon" is backed by evidence the
// channel could have shown one.
async function assertNoWarningBalloon(page: Page, what: string): Promise<void> {
  // Read instantaneously, not via toHaveCount(0): that form retries for the default timeout,
  // and a real balloon autohides after 5 s (balloon.dart), so retrying could wait one out.
  const warn = page.locator(WARNING_BALLOON);
  const n = await warn.count();
  const seen = n > 0 ? await warn.first().innerText() : '';
  expect(n, `${what} must not raise a warning balloon; first warning: "${seen}"`).toBe(0);

  await proveBalloonChannel(page, 'formula-refreshing');
}

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

    await softStep('Setup: open System:DemoFiles/demog.csv with datasync provenance', async () => {
      const opened = await openTableFromFile(page, 'System:DemoFiles/demog.csv');
      await page.locator('[name="viewer-Grid"]').waitFor({timeout: 60_000});
      await page.waitForTimeout(1000);

      await assertProvenanceScript(page, 'files', opened.script);
      const cols = await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        return df ? df.columns.names() : [];
      });
      expect(cols).toContain('WEIGHT');
    });

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
      await page.waitForTimeout(1000); 
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
      expect(post.w).toBeCloseTo(pre.w, 6); 
      expect(post.w2).toBeCloseTo(post.w + 200, 1);
      expect(post.w2).not.toBeCloseTo(pre.w2, 1);

      expect(post.w3).toBeCloseTo(post.w2 + 100, 1);
      expect(post.w3).not.toBeCloseTo(pre.w3, 1);
      expect(Number.isFinite(post.w4)).toBe(true);
      expect(post.w4).toBeCloseTo(Math.log10(post.w3) - 0.2, 3);
      expect(post.w2Tag).toContain('+ 200');
      // Absence is read first, while a balloon raised by the edit would still be standing; the
      // probe that follows it proves the channel live.
      await assertNoWarningBalloon(page, 'editing the Weight2 formula');
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
      expect(post.w2).toBeCloseTo(pre.w2, 6); 
      expect(post.w3).toBeCloseTo(post.w2 + 50, 1);
      expect(post.w3).not.toBeCloseTo(pre.w3, 1);
      expect(Number.isFinite(post.w4)).toBe(true);
      expect(post.w4).toBeCloseTo(Math.log10(post.w3) - 0.2, 3);
      expect(post.w3Tag).toContain('+ 50');
      await assertNoWarningBalloon(page, 'editing the Weight3 formula');
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

      expect(post.w2).toBeCloseTo(pre.w2, 6);
      expect(post.w3).toBeCloseTo(pre.w3, 6);
      expect(Number.isFinite(post.w4)).toBe(true);
      expect(post.w4).toBeCloseTo(Math.log10(post.w3) - 0.1, 3);
      expect(post.w4).not.toBeCloseTo(pre.w4, 3);
      expect(post.w4Tag).toContain('- 0.1');
      await assertNoWarningBalloon(page, 'editing the Weight4 formula');
    });

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
        };
      }, projectId);
      expect(reopen.ok).toBe(true);

      expect(reopen.hasWeight2).toBe(true);
      expect(reopen.hasWeight3).toBe(true);
      expect(reopen.hasWeight4).toBe(true);

      expect(reopen.w2Formula.length).toBeGreaterThan(0);
      expect(reopen.w3Formula.length).toBeGreaterThan(0);
      expect(reopen.w4Formula.length).toBeGreaterThan(0);

      expect(reopen.w2Formula).toContain('+ 200');
      expect(reopen.w3Formula).toContain('+ 50');
      expect(reopen.w4Formula).toContain('- 0.1');

      expect(reopen.w).toBeCloseTo(preSaveSnapshot.w, 6);
      expect(reopen.w2v).toBeCloseTo(reopen.w + 200, 1);
      expect(reopen.w3v).toBeCloseTo(reopen.w2v + 50, 1);
      expect(reopen.w4v).toBeCloseTo(Math.log10(reopen.w3v) - 0.1, 3);
      await assertNoWarningBalloon(page, 'reopening the saved project');
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
      } catch (_) {  }
      try { (window as any).grok.shell.closeAll(); } catch (_) {}
    }).catch(() => {});
  }

  finishSpec();
});

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

  await page.evaluate((n: string) => {
    const input = document.querySelector('[name="input-Add-New-Column---Name"]') as HTMLInputElement | null;
    if (!input) throw new Error('Name input not found');
    const setter = Object.getOwnPropertyDescriptor(window.HTMLInputElement.prototype, 'value')!.set!;
    setter.call(input, n);
    input.dispatchEvent(new Event('input', {bubbles: true}));
    input.dispatchEvent(new Event('change', {bubbles: true}));
  }, name);
  await page.waitForTimeout(150);

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

  const doc = composed.doc || '';
  const firstToken = formula.split(/\s+/)[0];
  if (firstToken)
    expect(doc).toContain(firstToken);
}

async function clickAddNewColumnOK(page: any): Promise<void> {
  const dlg = page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).first();
  await dlg.locator('[name="button-Add-New-Column---OK"]').first().click();

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

// AddNewColumnDialog picks its CodeMirror host class by MODE — `add-new-column-widget-cm-div`
// when constructed with a widget, `add-new-column-dialog-cm-div` otherwise (add-new-column.ts:211).
// PowerPack:formulaWidget builds it in widget mode (package.ts:195-210), so the panel-hosted editor
// carries the widget class; the dialog helper above keeps the dialog one. Polling for the dialog class
// here matched nothing on either path, which is why the panel fallback appeared to fail too.
async function editFormulaViaInfoPanel(page: any, columnName: string, newFormula: string): Promise<void> {

  await page.evaluate((n: string) => {
    const grok = (window as any).grok;
    const df = grok.shell.tv?.dataFrame;
    if (!df) throw new Error('editFormulaViaInfoPanel: no active dataframe');
    const col = df.col(n);
    if (!col) throw new Error(`editFormulaViaInfoPanel: column "${n}" not found`);
    grok.shell.o = col;
  }, columnName);
  await page.waitForTimeout(500); 

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
  await page.waitForTimeout(500); 

  let widgetCmFound = false;
  if (accordionPathWorked) {
    for (let i = 0; i < 25; i++) {
      widgetCmFound = await page.evaluate(() => {
        const propPanel = document.querySelector('.grok-prop-panel');
        return !!propPanel?.querySelector('.add-new-column-widget-cm-div .cm-content');
      });
      if (widgetCmFound) break;
      await page.waitForTimeout(200);
    }
  }

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
      host.setAttribute('data-fr-fallback', '1'); 
      propPanel?.appendChild(host);
      void DG; 
    }, columnName);
    for (let i = 0; i < 25; i++) {
      widgetCmFound = await page.evaluate(() => {
        const propPanel = document.querySelector('.grok-prop-panel');
        return !!propPanel?.querySelector('.add-new-column-widget-cm-div .cm-content');
      });
      if (widgetCmFound) break;
      await page.waitForTimeout(200);
    }
  }

  if (!widgetCmFound)
    throw new Error(`editFormulaViaInfoPanel: Formula widget CM host (.grok-prop-panel .add-new-column-widget-cm-div .cm-content) did not render for column "${columnName}"`);

  await page.evaluate(() => {
    const pp = document.querySelector('.grok-prop-panel');
    if (!pp) return;
    const cm = pp.querySelector('.add-new-column-widget-cm-div .cm-content') as HTMLElement | null;
    if (cm && cm.offsetWidth === 0 && cm.offsetHeight === 0) {
      const headers = Array.from(pp.querySelectorAll('.d4-accordion-pane-header, .d4-accordion-title'));
      const formulaHeader = headers.find((h) => (h.textContent || '').trim() === 'Formula') as HTMLElement | undefined;
      if (formulaHeader) formulaHeader.click();
    }
  });
  await page.waitForTimeout(400);
  const panelCm = page.locator(
    '.grok-prop-panel .add-new-column-widget-cm-div .cm-content').first();
  await panelCm.waitFor({timeout: 15_000, state: 'visible'});
  await panelCm.click({force: true});
  await page.waitForTimeout(200);
  let composed: {ok: boolean; doc?: string} = {ok: false};
  for (let i = 0; i < 10; i++) {
    composed = await page.evaluate((f: string) => {
      const propPanel = document.querySelector('.grok-prop-panel');
      const cmContent = propPanel?.querySelector(
        '.add-new-column-widget-cm-div .cm-content') as HTMLElement | null;
      if (!cmContent) return {ok: false};
      const view = (cmContent as any).cmTile?.view ?? null;
      if (!view) return {ok: false};
      view.dispatch({changes: {from: 0, to: view.state.doc.length, insert: f}});
      return {ok: true, doc: view.state.doc.toString()};
    }, newFormula);
    if (composed.ok) break;
    await page.waitForTimeout(200);
  }

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
        '.add-new-column-widget-cm-div .cm-content') as HTMLElement | null;
      if (!cmContent) return {ok: false};
      const view = (cmContent as any).cmTile?.view ?? null;
      const doc = view ? view.state.doc.toString() : (cmContent.innerText || '');
      return {ok: true, doc};
    });
  }
  if (!composed.ok)
    throw new Error('editFormulaViaInfoPanel: CodeMirror cmTile.view not exposed on in-panel widget even after keyboard fallback');

  const doc = composed.doc || '';
  const firstToken = newFormula.split(/\s+/)[0];
  if (firstToken)
    expect(doc).toContain(firstToken);

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
  await page.waitForTimeout(400); 

  await page.evaluate(() => {
    document.querySelectorAll('[data-fr-fallback="1"]').forEach((el) => el.remove());
  });
}
