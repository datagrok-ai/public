import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
test.use(specTestOptions);
const CM_SELECTOR = '.d4-dialog .add-new-column-dialog-cm-div .cm-content';
async function dispatchEditorReplace(
  page: Page, text: string,
): Promise<{ ok: boolean; doc: string }> {
  const cm = page.locator(CM_SELECTOR).first();
  await cm.waitFor({timeout: 15_000, state: 'visible'});
  await cm.click();
  await page.waitForTimeout(200);
  let result: {ok: boolean; doc: string} = {ok: false, doc: ''};
  for (let i = 0; i < 10; i++) {
    result = await page.evaluate((sel) => {
      const cmDiv = document.querySelector(sel) as HTMLElement | null;
      if (!cmDiv) return {ok: false, doc: ''};
      const tileView = (cmDiv as any)?.cmTile?.view;
      const legacyOnDiv = (cmDiv as any)?.cmView?.view;
      const legacyOnParent = (cmDiv.parentElement as any)?.cmView?.view;
      const view = tileView ?? legacyOnDiv ?? legacyOnParent ?? null;
      if (!view) return {ok: false, doc: ''};
      const doc = view.state.doc.toString();
      return {ok: true, doc};
    }, CM_SELECTOR);
    if (result.ok) break;
    await page.waitForTimeout(200);
  }
  if (!result.ok) return result;
  result = await page.evaluate((args) => {
    const cmDiv = document.querySelector(args.sel) as HTMLElement | null;
    if (!cmDiv) return {ok: false, doc: ''};
    const tileView = (cmDiv as any)?.cmTile?.view;
    const legacyOnDiv = (cmDiv as any)?.cmView?.view;
    const legacyOnParent = (cmDiv.parentElement as any)?.cmView?.view;
    const view = tileView ?? legacyOnDiv ?? legacyOnParent ?? null;
    if (!view) return {ok: false, doc: ''};
    view.dispatch({changes: {
      from: 0, to: view.state.doc.length, insert: args.text,
    }});
    return {ok: true, doc: view.state.doc.toString()};
  }, {sel: CM_SELECTOR, text});
  return result;
}
async function keyboardFallbackReplace(
  page: Page, text: string,
): Promise<{ doc: string }> {
  const cm = page.locator(CM_SELECTOR).first();
  await cm.click();
  await page.waitForTimeout(150);
  await page.keyboard.press('Control+A');
  await page.keyboard.press('Delete');
  await page.waitForTimeout(100);
  await page.keyboard.type(text, {delay: 30});
  await page.waitForTimeout(200);
  const doc = await page.evaluate((sel) => {
    const cmDiv = document.querySelector(sel) as HTMLElement | null;
    if (!cmDiv) return '';
    const tileView = (cmDiv as any)?.cmTile?.view;
    const legacyOnDiv = (cmDiv as any)?.cmView?.view;
    const legacyOnParent = (cmDiv.parentElement as any)?.cmView?.view;
    const view = tileView ?? legacyOnDiv ?? legacyOnParent ?? null;
    return view ? view.state.doc.toString() : (cmDiv.innerText || '');
  }, CM_SELECTOR);
  return {doc};
}
async function composeFormula(
  page: Page, text: string,
): Promise<string> {
  const cm = page.locator(CM_SELECTOR).first();
  await cm.waitFor({timeout: 15_000, state: 'visible'});
  await cm.click();
  await page.waitForTimeout(200);
  await page.keyboard.press('Control+A');
  await page.keyboard.press('Delete');
  await page.waitForTimeout(100);
  const dispatched = await dispatchEditorReplace(page, text);
  let doc = dispatched.doc;
  if (!dispatched.ok || !doc.includes(text)) {
    const kb = await keyboardFallbackReplace(page, text);
    doc = kb.doc;
  }
  await page.waitForTimeout(500);
  return doc;
}
async function readHighlightedColumnTokens(page: Page): Promise<string[]> {
  return page.evaluate((sel) => {
    const cmDiv = document.querySelector(sel) as HTMLElement | null;
    if (!cmDiv) return [] as string[];
    return Array.from(cmDiv.querySelectorAll('.cm-column-name'))
      .map((el) => (el.textContent || '').trim())
      .filter((t) => t.length > 0);
  }, CM_SELECTOR);
}
async function readFirstHighlightBlueness(page: Page): Promise<
  { hasSpan: boolean; rgb: string | null; isBlue: boolean }
> {
  return page.evaluate((sel) => {
    const cmDiv = document.querySelector(sel) as HTMLElement | null;
    if (!cmDiv) return {hasSpan: false, rgb: null, isBlue: false};
    const span = cmDiv.querySelector('.cm-column-name') as HTMLElement | null;
    if (!span) return {hasSpan: false, rgb: null, isBlue: false};
    const cs = window.getComputedStyle(span);
    const rgb = cs.color;
    const m = rgb.match(/rgba?\((\d+),\s*(\d+),\s*(\d+)/);
    if (!m) return {hasSpan: true, rgb, isBlue: false};
    const r = Number(m[1]); const g = Number(m[2]); const b = Number(m[3]);
    const isBlue = b > r && b >= g && b > 0;
    return {hasSpan: true, rgb, isBlue};
  }, CM_SELECTOR);
}
async function readDistinctHighlightedColumns(page: Page): Promise<string[]> {
  return page.evaluate((sel) => {
    const cmDiv = document.querySelector(sel) as HTMLElement | null;
    if (!cmDiv) return [] as string[];
    const distinct = new Set<string>();
    for (const el of Array.from(cmDiv.querySelectorAll('.cm-column-name'))) {
      const t = (el.textContent || '').trim();
      const m = t.match(/^\$[\{\[](.+)[\}\]]$/);
      if (m) distinct.add(m[1]);
    }
    return Array.from(distinct);
  }, CM_SELECTOR);
}
async function dispatchAndCaptureErrors(
  page: Page, text: string,
): Promise<{ doc: string; errors: string[]; threw: boolean }> {
  const cm = page.locator(CM_SELECTOR).first();
  await cm.waitFor({timeout: 15_000, state: 'visible'});
  await cm.click();
  await page.waitForTimeout(200);
  await page.keyboard.press('Control+A');
  await page.keyboard.press('Delete');
  await page.waitForTimeout(150);
  const result = await page.evaluate(async (args) => {
    const errs: string[] = [];
    const onErr = (ev: ErrorEvent) =>
      errs.push(String(ev.message || (ev as any).error || ev));
    window.addEventListener('error', onErr);
    const origConsoleError = console.error;
    console.error = function(...a: any[]) {
      errs.push(a.map((x) => String(x)).join(' '));
      return origConsoleError.apply(console, a as any);
    };
    let threw = false;
    let doc = '';
    try {
      const cmDiv = document.querySelector(args.sel) as HTMLElement | null;
      if (!cmDiv) { errs.push('no cm-content host'); throw new Error('no host'); }
      const tileView = (cmDiv as any)?.cmTile?.view;
      const legacyOnDiv = (cmDiv as any)?.cmView?.view;
      const legacyOnParent = (cmDiv.parentElement as any)?.cmView?.view;
      const view = tileView ?? legacyOnDiv ?? legacyOnParent ?? null;
      if (!view) { errs.push('no EditorView accessor'); throw new Error('no view'); }
      view.dispatch({changes: {from: 0, to: view.state.doc.length, insert: args.text}});
      doc = view.state.doc.toString();
    } catch (e) {
      threw = true;
      errs.push('DISPATCH THREW: ' + String(e));
    }
    await new Promise((r) => setTimeout(r, 1200));
    window.removeEventListener('error', onErr);
    console.error = origConsoleError;
    return {doc, errors: errs, threw};
  }, {sel: CM_SELECTOR, text});
  await page.waitForTimeout(300);
  return result;
}
test('PowerPack: Add new column - column-name highlight (GROK-17004 invariant)', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);
  await page.evaluate(async () => {
    const grok = (window as any).grok;
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    try { grok.shell.closeAll(); } catch (_) {}
    const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
    grok.shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
  });
  await page.locator('[name="viewer-Grid"]').waitFor({timeout: 60_000});
  await page.waitForTimeout(1000);
  const cols = await page.evaluate(() => {
    const df = (window as any).grok.shell.tv?.dataFrame;
    return df ? df.columns.names() : [];
  });
  expect(cols).toContain('AGE');
  await softStep('Setup: open Add New Column dialog via toolbar icon', async () => {
    const icon = page.locator('[name="icon-add-new-column"]').first();
    await icon.waitFor({timeout: 30_000, state: 'visible'});
    await icon.click({timeout: 10_000});
    const dlg = page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).first();
    await dlg.waitFor({timeout: 30_000});
    await expect(dlg).toBeVisible();
  });
  const cm = page.locator(CM_SELECTOR).first();
  await cm.waitFor({timeout: 15_000, state: 'visible'});
  await softStep('Scenario 1 / Step 1: insert "Abs(${AGE})" into the formula editor', async () => {
    const doc = await composeFormula(page, 'Abs(${AGE})');
    expect(doc).toContain('Abs(${AGE})');
  });
  await softStep('Scenario 1 / Step 2 + GROK-17004 INVARIANT: "AGE" highlighted in blue', async () => {
    const tokens = await readHighlightedColumnTokens(page);
    expect(tokens.length).toBeGreaterThan(0);
    const hasAge = tokens.some((t) => /\$\{AGE\}/i.test(t) || /\bAGE\b/i.test(t));
    expect(hasAge).toBe(true);
    const blueness = await readFirstHighlightBlueness(page);
    expect(blueness.hasSpan).toBe(true);
    expect(blueness.isBlue).toBe(true);
  });
  await softStep('Scenario 2 / Step 1: insert "Avg($[AGE])" into the formula editor', async () => {
    const doc = await composeFormula(page, 'Avg($[AGE])');
    expect(doc).toContain('Avg($[AGE])');
  });
  await softStep('Scenario 2 / Step 2: "AGE" inside $[...] highlighted in blue', async () => {
    const tokens = await readHighlightedColumnTokens(page);
    expect(tokens.length).toBeGreaterThan(0);
    const hasAgeBracket = tokens.some((t) => /\$\[AGE\]/i.test(t) || /\bAGE\b/i.test(t));
    expect(hasAgeBracket).toBe(true);
    const blueness = await readFirstHighlightBlueness(page);
    expect(blueness.hasSpan).toBe(true);
    expect(blueness.isBlue).toBe(true);
  });
  await softStep('Scenario 3 / Step 1-2: insert "Round(${AGE})" via dispatch (autocomplete-completion end-state)', async () => {
    const doc = await composeFormula(page, 'Round(${AGE})');
    expect(doc).toContain('Round(${AGE})');
  });
  await softStep('Scenario 3 / Step 3: inserted ${<column>} reference highlighted in blue', async () => {
    const tokens = await readHighlightedColumnTokens(page);
    expect(tokens.length).toBeGreaterThan(0);
    const hasRef = tokens.some((t) => /\$\{[^}]+\}/.test(t));
    expect(hasRef).toBe(true);
    const blueness = await readFirstHighlightBlueness(page);
    expect(blueness.hasSpan).toBe(true);
    expect(blueness.isBlue).toBe(true);
  });
  await softStep('Scenario 4 / Step 1: clear editor; type function prefix "Sin("', async () => {
    await composeFormula(page, '');
    await cm.click();
    await page.waitForTimeout(150);
    await page.keyboard.type('Sin(', {delay: 50});
    await page.waitForTimeout(200);
    await page.keyboard.press('Escape');
    await page.waitForTimeout(120);
  });
  await softStep('Scenario 4 / Step 2: insert ${AGE} reference at end (drag-drop equivalent end state)', async () => {
    const result = await page.evaluate((sel) => {
      const cmDiv = document.querySelector(sel) as HTMLElement | null;
      if (!cmDiv) return {ok: false, doc: ''};
      const tileView = (cmDiv as any)?.cmTile?.view;
      const legacyOnDiv = (cmDiv as any)?.cmView?.view;
      const legacyOnParent = (cmDiv.parentElement as any)?.cmView?.view;
      const view = tileView ?? legacyOnDiv ?? legacyOnParent ?? null;
      if (!view) return {ok: false, doc: ''};
      view.dispatch({changes: {
        from: view.state.doc.length,
        to: view.state.doc.length,
        insert: '${AGE})',
      }});
      return {ok: true, doc: view.state.doc.toString()};
    }, CM_SELECTOR);
    expect(result.ok).toBe(true);
    expect(result.doc).toContain('${AGE}');
    await page.waitForTimeout(500);
  });
  await softStep('Scenario 4 / Step 3: inserted ${AGE} reference highlighted in blue', async () => {
    const tokens = await readHighlightedColumnTokens(page);
    expect(tokens.length).toBeGreaterThan(0);
    const hasAge = tokens.some((t) => /\$\{AGE\}/i.test(t) || /\bAGE\b/i.test(t));
    expect(hasAge).toBe(true);
    const blueness = await readFirstHighlightBlueness(page);
    expect(blueness.hasSpan).toBe(true);
    expect(blueness.isBlue).toBe(true);
  });
  await softStep('Scenario 5 / Step 1: switch active dataset to SPGI', async () => {
    await page.evaluate(async () => {
      const grok = (window as any).grok;
      try { grok.shell.closeAll(); } catch (_) {  }
      const df = await grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv');
      grok.shell.addTableView(df);
      await new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 4000);
      });
    });
    await page.locator('[name="viewer-Grid"]').waitFor({timeout: 60_000});
    await page.waitForTimeout(2000);
    const spgiCols = await page.evaluate(() => {
      const df = (window as any).grok.shell.tv?.dataFrame;
      return df ? df.columns.names() : [];
    });
    for (const c of ['Whole blood assay 1', 'Route Admin', 'Chemical Space X', 'Average Mass', 'Species'])
      expect(spgiCols).toContain(c);
  });
  await softStep('Scenario 5 / Step 2: open Add New Column dialog against SPGI', async () => {
    const icon = page.locator('[name="icon-add-new-column"]').first();
    await icon.waitFor({timeout: 30_000, state: 'visible'});
    await icon.click({timeout: 10_000});
    const dlg = page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).first();
    await dlg.waitFor({timeout: 30_000});
    await expect(dlg).toBeVisible();
    const cm5 = page.locator(CM_SELECTOR).first();
    await cm5.waitFor({timeout: 15_000, state: 'visible'});
  });
  const GROK_17004_FORMULA =
    'if(${Whole blood assay 1} != null, ${Whole blood assay 1}, ' +
    'if(${Route Admin}=="PO", ${Whole blood assay 1} / ${Chemical Space X} ' +
    '* 100 / 6 / ${Average Mass} * 1000000.0,null))/' +
    "if(Contains(${Species}, 'Rat') || Contains(${Species}, 'Rat Legacy'), 80, " +
    "if(Contains(${Species}, 'Mouse'), 125, " +
    'if(${Species}=="Dog", 30.9, if(${Species}=="Monkey", 43.6, ' +
    'if(${Species}=="Minipig", 39, null)))))*100';
  let pasteResult: {doc: string; errors: string[]; threw: boolean} =
    {doc: '', errors: [], threw: false};
  await softStep('Scenario 5 / Step 3-4 + GROK-17004 INVARIANT: paste complex formula; handler does NOT throw', async () => {
    pasteResult = await dispatchAndCaptureErrors(page, GROK_17004_FORMULA);
    expect(pasteResult.doc).toBe(GROK_17004_FORMULA);
    expect(pasteResult.threw).toBe(false);
    const grokError = pasteResult.errors.find((e) =>
      /TypeError|Cannot read propert|reading 'to'|getColumnNamesAndSelections|DISPATCH THREW/i.test(e));
    expect(grokError, `unexpected error during paste: ${pasteResult.errors.join(' | ')}`).toBeUndefined();
  });
  await softStep('Scenario 5 / Step 5 + GROK-17004 INVARIANT: five distinct columns highlighted in blue', async () => {
    const tokens = await readHighlightedColumnTokens(page);
    // At least one highlight span surfaced (the addColHighlight pipeline ran).
    expect(tokens.length).toBeGreaterThan(0);
    // Each of the five referenced columns is highlighted (distinct set).
    const distinct = await readDistinctHighlightedColumns(page);
    for (const c of ['Whole blood assay 1', 'Route Admin', 'Chemical Space X', 'Average Mass', 'Species'])
      expect(distinct, `column "${c}" not highlighted; distinct=[${distinct.join(', ')}]`).toContain(c);
    // At least one rendered span has a blue computed color (the
    // getColumnNamesAndSelections -> addColHighlight pipeline completed).
    const blueness = await readFirstHighlightBlueness(page);
    expect(blueness.hasSpan).toBe(true);
    expect(blueness.isBlue).toBe(true);
  });
  // ---- Cleanup: cancel the dialog and clear shell state ----
  await page.evaluate(() => {
    const cancel = document.querySelector(
      '.d4-dialog [name="button-Add-New-Column---CANCEL"]',
    ) as HTMLElement | null;
    if (cancel) cancel.click();
    const anyCancel = document.querySelector(
      '.d4-dialog [name="button-CANCEL"]',
    ) as HTMLElement | null;
    if (anyCancel) anyCancel.click();
  }).catch(() => {});
  await page.evaluate(() => {
    try { (window as any).grok.shell.closeAll(); } catch (_) { /* best effort */ }
  }).catch(() => {});
  finishSpec();
});
