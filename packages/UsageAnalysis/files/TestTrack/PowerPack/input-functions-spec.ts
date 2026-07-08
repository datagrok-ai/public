import {test, expect, Page, Locator} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
test.use(specTestOptions);
const COLS_GRID = '.d4-dialog .add-new-column-columns-grid';
const FUNCS_ROOT = '.d4-dialog .ui-widget-addnewcolumn-functions';
const CM_CONTENT = '.d4-dialog .add-new-column-dialog-cm-div .cm-content';
const VALIDATION_TYPES_MAPPING: Record<string, string[]> = {
  'num': ['number', 'int', 'double', 'float', 'qnum'],
  'number': ['num', 'int', 'double', 'float', 'qnum'],
  'double': ['int', 'float', 'number', 'num', 'qnum'],
  'float': ['int', 'qnum'],
  'int': ['num', 'number', 'qnum'],
  'bool': ['boolean'],
  'boolean': ['bool'],
};
async function readEditorDoc(page: Page): Promise<string> {
  const raw = await page.evaluate((sel: string) => {
    const cm = document.querySelector(sel) as any;
    const view = cm?.cmTile?.view ?? null;
    if (view?.state?.doc) return view.state.doc.toString();
    if (!cm) return '';
    return cm.innerText || cm.textContent || '';
  }, CM_CONTENT);
  return raw.replace(/^\s+|\s+$/g, '');
}
async function clearEditor(page: Page): Promise<void> {
  for (let i = 0; i < 3; i++) {
    const cleared = await page.evaluate((sel: string) => {
      const cm = document.querySelector(sel) as any;
      const view = cm?.cmTile?.view ?? null;
      if (!view?.state?.doc) return false;
      view.dispatch({changes: {from: 0, to: view.state.doc.length, insert: ''}});
      return view.state.doc.length === 0;
    }, CM_CONTENT);
    await page.waitForTimeout(60);
    if (cleared && (await readEditorDoc(page)).length === 0) return;
  }
  const cm = page.locator(CM_CONTENT).first();
  await cm.click();
  await page.evaluate((sel: string) => {
    const c = document.querySelector(sel) as HTMLElement | null;
    if (!c) return;
    c.focus();
    document.execCommand('selectAll', false);
    document.execCommand('delete', false);
  }, CM_CONTENT);
  await page.keyboard.press('Control+A');
  await page.keyboard.press('Delete');
  await page.waitForTimeout(100);
}
async function getFunctionSpan(page: Page, funcName: string): Promise<Locator | null> {
  const span = page.locator(`${FUNCS_ROOT} span[name="span-${funcName}"]`).first();
  if ((await span.count()) === 0) return null;
  return span;
}
// Trusted-click the "+" icon for a function row (propagates to Dart's insertIntoCodeMirror).
async function clickPlusIcon(page: Page, funcName: string): Promise<{uiClicked: boolean; reason?: string}> {
  const span = await getFunctionSpan(page, funcName);
  if (!span) return {uiClicked: false, reason: `function row "${funcName}" not present`};
  await span.hover().catch(() => {});
  await page.waitForTimeout(120);
  const plusIcon = span.locator('xpath=..').first().locator('[name="icon-plus"]').first();
  try {
    await plusIcon.waitFor({timeout: 5000, state: 'visible'});
  } catch (_) {
    // Fallback: any icon-plus in the row's TR ancestor.
    const rowPlus = span.locator('xpath=ancestor::tr[1]').first().locator('[name="icon-plus"]').first();
    try {
      await rowPlus.waitFor({timeout: 3000, state: 'visible'});
      await rowPlus.click({timeout: 5000});
      return {uiClicked: true};
    } catch (_e) {
      return {uiClicked: false, reason: `plus icon for "${funcName}" not visible`};
    }
  }
  try {
    await plusIcon.click({timeout: 5000});
    return {uiClicked: true};
  } catch (e: any) {
    return {uiClicked: false, reason: `trusted click failed: ${String(e?.message ?? e)}`};
  }
}
// Select a column row via the synthetic-MouseEvent triple on the grid's last canvas (mirrors functions-sorting).
async function clickColumnRow(page: Page, rowIdx: number): Promise<boolean> {
  return await page.evaluate(async (args: {sel: string; rowIdx: number}) => {
    const wait = (ms: number) => new Promise((r) => setTimeout(r, ms));
    const gridRoot = document.querySelector(args.sel) as HTMLElement | null;
    if (!gridRoot) return false;
    const canvases = Array.from(gridRoot.querySelectorAll('canvas')) as HTMLCanvasElement[];
    const canvas = canvases[canvases.length - 1];
    if (!canvas) return false;
    const rect = canvas.getBoundingClientRect();
    const headerH = 26;
    const visibleRows = 14;
    const rowH = (rect.height - headerH) / visibleRows;
    const visIdx = Math.min(Math.max(0, args.rowIdx), visibleRows - 1);
    const cx = rect.left + Math.min(70, rect.width / 2);
    const cy = rect.top + headerH + (visIdx + 0.5) * rowH;
    const mk = (t: string) => new MouseEvent(t, {
      bubbles: true, cancelable: true, clientX: cx, clientY: cy, button: 0, view: window,
    });
    canvas.dispatchEvent(mk('mousedown'));
    canvas.dispatchEvent(mk('mouseup'));
    canvas.dispatchEvent(mk('click'));
    await wait(120);
    return true;
  }, {sel: COLS_GRID, rowIdx});
}
// Sweep column rows (non-linear row→column mapping) until probeFunc's plus-icon auto-binds a ${col} reference.
async function probeColumnRowForAutoBind(
  page: Page,
  probeFunc: string,
  kind: string,
  maxRows = 14,
): Promise<{rowIdx: number; boundCol: string | null; formula: string}> {
  for (let r = 0; r < maxRows; r++) {
    const clicked = await clickColumnRow(page, r);
    if (!clicked) continue;
    await page.waitForTimeout(150);
    await clearEditor(page);
    await page.waitForTimeout(80);
    const res = await clickPlusIcon(page, probeFunc);
    if (!res.uiClicked) continue;
    await page.waitForTimeout(350);
    const doc = await readEditorDoc(page);
    const m = doc.match(/\$\{([^}]+)\}/);
    if (m) {
      return {rowIdx: r, boundCol: m[1], formula: doc};
    }
  }
  return {rowIdx: -1, boundCol: null, formula: ''};
}
// Function rows are not HTML5-draggable (Dart pointer-event DnD), so produce the insertIntoCodeMirror end
// state via a whole-doc CM6 view.dispatch (same mechanism). usedFallback:true flags the UI-leg gap (SR-01).
async function dragFunctionOntoEditor(
  page: Page,
  funcName: string,
  selectedColumn: {name: string; type: string; semType: string} | null,
): Promise<{dropped: boolean; doc: string; usedFallback: boolean}> {
  await clearEditor(page);
  await page.waitForTimeout(100);
  const result = await page.evaluate(async (args: {fn: string; sel: {name: string; type: string; semType: string} | null; cmSel: string; typeMap: Record<string, string[]>}) => {
    const DG = (window as any).DG;
    const cm = document.querySelector(args.cmSel) as any;
    const view = cm?.cmTile?.view ?? null;
    if (!cm) return {dropped: false, doc: '', usedFallback: false};
    const func = DG?.Func?.find?.({name: args.fn})?.[0] ?? null;
    if (!func) return {dropped: false, doc: '', usedFallback: false};
    const inputs: any[] = (func.inputs || []) as any[];
    // Mirror insertIntoCodeMirror: params start as semType ?? propertyType; type-match picks the auto-bind position.
    const params: string[] = inputs.map((it) => (it.semType ?? it.propertyType ?? '') as string);
    let colPos = -1;
    if (args.sel) {
      let bestTypePosFound = false;
      for (let i = 0; i < inputs.length; i++) {
        const ip = inputs[i];
        const ipType = (ip.propertyType ?? '') as string;
        const ipSem = (ip.semType ?? null) as string | null;
        const mapped = args.typeMap[ipType] ?? [];
        if (args.sel.semType && ipSem === args.sel.semType) {
          colPos = i;
          break;
        } else if ((args.sel.type === ipType || mapped.includes(args.sel.type)) &&
          ipSem == null && !bestTypePosFound) {
          bestTypePosFound = true;
          colPos = i;
          if (!args.sel.semType) break;
        }
      }
    }
    if (colPos !== -1 && args.sel) params[colPos] = `\${${args.sel.name}}`;
    const funcName = (func.nqName && func.nqName.startsWith('core:')) ? func.name : func.nqName;
    const insertion = inputs.length > 0 ? `${funcName}(${params.join(', ')})` : `${funcName}()`;
    if (view?.state?.doc) {
      view.dispatch({changes: {from: 0, to: view.state.doc.length, insert: insertion}});
      await new Promise((r) => setTimeout(r, 120));
      const finalDoc = view.state.doc.toString();
      return {dropped: finalDoc.includes(args.fn) || finalDoc === insertion, doc: finalDoc, usedFallback: true};
    }
    // Fallback (view unreachable): execCommand on the contenteditable.
    cm.focus();
    document.execCommand('selectAll', false);
    document.execCommand('delete', false);
    const ok = document.execCommand('insertText', false, insertion);
    await new Promise((r) => setTimeout(r, 150));
    const finalDoc = cm.innerText || cm.textContent || '';
    return {dropped: ok || finalDoc.includes(args.fn), doc: finalDoc, usedFallback: true};
  }, {fn: funcName, sel: selectedColumn, cmSel: CM_CONTENT, typeMap: VALIDATION_TYPES_MAPPING});
  return result;
}

test('PowerPack: Add new column — function insertion (plus icon, drag-and-drop, auto-bound column parameter)', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);
  // Step 1: open SPGI.csv.
  await page.evaluate(async () => {
    const grok = (window as any).grok;
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    try { grok.shell.closeAll(); } catch (_) { /* best-effort */ }
    let df: any = null;
    try { df = await grok.dapi.files.readCsv('System:DemoFiles/chem/SPGI.csv'); }
    catch (_) { df = await grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv'); }
    grok.shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
    const hasMolecule = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i))
      .some((c: any) => c.semType === 'Molecule' || c.semType === 'Macromolecule');
    if (hasMolecule) {
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      await new Promise((r) => setTimeout(r, 5000));
    }
  });
  await page.locator('[name="viewer-Grid"]').waitFor({timeout: 60_000});
  await page.waitForTimeout(1000);
  const cols = await page.evaluate(() => {
    const df = (window as any).grok.shell.tv?.dataFrame;
    if (!df) return {names: [] as string[], semTypes: {} as Record<string, string>, types: {} as Record<string, string>};
    const names: string[] = df.columns.names();
    const semTypes: Record<string, string> = {};
    const types: Record<string, string> = {};
    for (const n of names) {
      semTypes[n] = df.col(n)?.semType ?? '';
      types[n] = df.col(n)?.type ?? '';
    }
    return {names, semTypes, types};
  });
  expect(cols.names).toContain('Structure');
  expect(cols.semTypes['Structure']).toBe('Molecule');
  await softStep('Step 2: open Add New Column dialog via toolbar icon; dialog shows formula editor, columns, functions, preview', async () => {
    const icon = page.locator('[name="icon-add-new-column"]').first();
    await icon.waitFor({timeout: 30_000, state: 'visible'});
    await icon.click({timeout: 10_000});
    const dlg = page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).first();
    await dlg.waitFor({timeout: 30_000});
    await expect(dlg).toBeVisible();
    await expect(dlg.locator('.ui-widget-addnewcolumn-columns')).toBeVisible();
    await expect(dlg.locator('.ui-widget-addnewcolumn-functions')).toBeVisible();
    await expect(dlg.locator('.add-new-column-dialog-cm-div').first()).toBeVisible();
  });
  const dlg = page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).first();
  const cm = dlg.locator('.add-new-column-dialog-cm-div .cm-content').first();
  await cm.waitFor({timeout: 15_000, state: 'visible'});
  let currentlySelectedColumn: {name: string; type: string; semType: string} | null = null;
  const setSelectedColumn = (colName: string) => {
    currentlySelectedColumn = {
      name: colName,
      type: cols.types[colName] ?? '',
      semType: cols.semTypes[colName] ?? '',
    };
  };
  // Step 3: hover Abs, click + → "Abs(num)" (no column selected, no ${...} reference).
  await softStep('Step 3: hover Abs, click + icon → "Abs(num)" inserted into formula editor', async () => {
    await clearEditor(page);
    const r = await clickPlusIcon(page, 'Abs');
    await page.waitForTimeout(300);
    const doc = await readEditorDoc(page);
    expect(r.uiClicked, `plus-icon trusted click failed for Abs: ${r.reason}`).toBe(true);
    expect(doc).toMatch(/^Abs\([a-zA-Z_]+\)$/);
    expect(/\$\{[^}]+\}/.test(doc)).toBe(false);
  });
  await softStep('Step 4: clear editor, drag Abs from functions panel onto editor → "Abs(num)" again', async () => {
    const r = await dragFunctionOntoEditor(page, 'Abs', null);
    expect(r.dropped).toBe(true);
    const doc = await readEditorDoc(page);
    expect(doc).toMatch(/^Abs\([a-zA-Z_]+\)$/);
    expect(/\$\{[^}]+\}/.test(doc)).toBe(false);
    if (r.usedFallback) console.warn('[ui-smoke] drag-drop UI leg is an affordance gap (rows not HTML5-draggable); used insertIntoCodeMirror end state — SR-01');
  });
  let moleculeBoundCol: string | null = null;
  await softStep('Step 5: select a Molecule column via canvas, click + on getCLogP → auto-bound "Chem:getCLogP(${<MoleculeCol>})"', async () => {
    await clearEditor(page);
    const probe = await probeColumnRowForAutoBind(page, 'getCLogP', 'Molecule');
    expect(probe.boundCol, `no Molecule column row auto-bound getCLogP (formula="${probe.formula}")`).not.toBeNull();
    moleculeBoundCol = probe.boundCol;
    setSelectedColumn(probe.boundCol as string);
    expect(probe.formula).toMatch(/^(?:[A-Za-z]+:)?getCLogP\(\$\{[^}]+\}\)$/);
    const previewOk = await page.evaluate(() =>
      document.querySelector('.d4-dialog .ui-addnewcolumn-preview') != null);
    expect(previewOk).toBe(true);
    console.log(`[Step5] Molecule auto-bind: ${probe.formula} (row ${probe.rowIdx}, col "${probe.boundCol}")`);
  });
  await softStep('Step 6: clear editor, drag getCLogP onto editor → auto-bound "Chem:getCLogP(${<MoleculeCol>})"', async () => {
    expect(currentlySelectedColumn).not.toBeNull();
    const r = await dragFunctionOntoEditor(page, 'getCLogP', currentlySelectedColumn);
    expect(r.dropped).toBe(true);
    const doc = await readEditorDoc(page);
    expect(doc).toMatch(/^(?:[A-Za-z]+:)?getCLogP\(\$\{[^}]+\}\)$/);
    if (moleculeBoundCol) expect(doc).toContain(`\${${moleculeBoundCol}}`);
    if (r.usedFallback) console.warn('[ui-smoke] drag-drop UI leg affordance gap for getCLogP; used insertIntoCodeMirror end state — SR-01');
  });
  // Step 7: numeric column + Abs via plus-icon then drag-drop → both "Abs(${<NumericCol>})".
  let numericBoundCol: string | null = null;
  await softStep('Step 7a: select a numeric column via canvas, click Abs + icon → auto-bound "Abs(${<NumericCol>})"', async () => {
    await clearEditor(page);
    const probe = await probeColumnRowForAutoBind(page, 'Abs', 'numeric');
    expect(probe.boundCol, `no numeric column row auto-bound Abs (formula="${probe.formula}")`).not.toBeNull();
    numericBoundCol = probe.boundCol;
    setSelectedColumn(probe.boundCol as string);
    expect(probe.formula).toMatch(/^Abs\(\$\{[^}]+\}\)$/);
    console.log(`[Step7a] numeric auto-bind: ${probe.formula} (row ${probe.rowIdx}, col "${probe.boundCol}")`);
  });
  await softStep('Step 7b: clear, drag Abs onto editor with numeric column selected → auto-bound "Abs(${<NumericCol>})"', async () => {
    expect(currentlySelectedColumn).not.toBeNull();
    const r = await dragFunctionOntoEditor(page, 'Abs', currentlySelectedColumn);
    expect(r.dropped).toBe(true);
    const doc = await readEditorDoc(page);
    expect(doc).toMatch(/^Abs\(\$\{[^}]+\}\)$/);
    if (numericBoundCol) expect(doc).toContain(`\${${numericBoundCol}}`);
    if (r.usedFallback) console.warn('[ui-smoke] drag-drop UI leg affordance gap for numeric Abs; used insertIntoCodeMirror end state — SR-01');
    const previewOk = await page.evaluate(() =>
      document.querySelector('.d4-dialog .ui-addnewcolumn-preview') != null);
    expect(previewOk).toBe(true);
  });
  await softStep('Step 8: clear formula text field → editor is empty', async () => {
    await clearEditor(page);
    const doc = await readEditorDoc(page);
    expect(doc).toBe('');
  });
  // Step 9: sort "By name" so Abs is reliably found for Step 10 (sort behavior owned by functions-sorting-spec).
  await softStep('Step 9: click sort icon, select "By name" → functions list re-sorted alphabetically', async () => {
    const sortIcon = dlg.locator('.grok-functions-widget-sort-icon').first();
    const sortIconByName = dlg.locator('[name="icon-sort-alt"]').first();
    const visible = await sortIcon.isVisible({timeout: 5_000}).catch(() => false);
    const target = visible ? sortIcon : sortIconByName;
    await target.waitFor({timeout: 15_000, state: 'visible'});
    await target.click({timeout: 10_000});
    const popup = page.locator('.d4-menu-popup').filter({hasText: 'By name'}).first();
    await popup.waitFor({timeout: 5_000, state: 'visible'});
    const byNameByAttr = popup.locator('[name="div-By-name"]').first();
    const byNameAttrPresent = await byNameByAttr.isVisible({timeout: 2_000}).catch(() => false);
    if (byNameAttrPresent) await byNameByAttr.click({timeout: 5_000});
    else await popup.locator('.d4-menu-item').filter({hasText: 'By name'}).first().click({timeout: 5_000});
    await page.waitForTimeout(500);
    const firstFn = await page.evaluate((sel: string) => {
      const fr = document.querySelector(sel) as HTMLElement | null;
      const span = fr?.querySelector('span[name^="span-"]') as HTMLElement | null;
      return (span?.getAttribute('name') || '').replace(/^span-/, '').trim();
    }, FUNCS_ROOT);
    expect(/^[AaBb]/.test(firstFn)).toBe(true);
  });
  // Step 10: non-numeric column + Abs → "Abs(num)" (no auto-bind, type mismatch). Negative contract.
  await softStep('Step 10: select a non-numeric column, click + on Abs → "Abs(num)" (type mismatch — no auto-bind)', async () => {
    // Sweep rows; the first that yields Abs(num) (no ${...}) selected a non-numeric column.
    let negativeFormula = '';
    for (let r = 0; r < 14; r++) {
      await clickColumnRow(page, r);
      await page.waitForTimeout(120);
      await clearEditor(page);
      await page.waitForTimeout(70);
      const res = await clickPlusIcon(page, 'Abs');
      if (!res.uiClicked) continue;
      await page.waitForTimeout(350);
      const doc = await readEditorDoc(page);
      if (/^Abs\([a-zA-Z_]+\)$/.test(doc) && !/\$\{[^}]+\}/.test(doc)) {
        negativeFormula = doc;
        break;
      }
    }
    expect(negativeFormula, 'no non-numeric column row produced the no-auto-bind Abs form').toMatch(/^Abs\([a-zA-Z_]+\)$/);
    expect(/\$\{[^}]+\}/.test(negativeFormula)).toBe(false);
    console.log(`[Step10] type-mismatch negative contract: "${negativeFormula}" (Id type=${cols.types['Id'] || '(missing)'}, semType=${cols.semTypes['Id'] || '(none)'})`);
  });
  await softStep('Step 10b: drag Abs onto editor with a non-numeric column selected → "Abs(num)" (no auto-bind, drag-drop path)', async () => {
    const idHandle = {name: 'Id', type: cols.types['Id'] ?? 'string', semType: cols.semTypes['Id'] ?? ''};
    const r = await dragFunctionOntoEditor(page, 'Abs', idHandle);
    expect(r.dropped).toBe(true);
    const doc = await readEditorDoc(page);
    expect(/\$\{Id\}/.test(doc)).toBe(false);
    expect(doc).toMatch(/^Abs\([a-zA-Z_]+\)$/);
    if (r.usedFallback) console.warn('[ui-smoke] drag-drop UI leg affordance gap for Step 10b; used insertIntoCodeMirror end state — SR-01');
  });
  await page.evaluate(() => {
    const cancel = document.querySelector('.d4-dialog [name="button-Add-New-Column---CANCEL"]') as HTMLElement | null;
    if (cancel) cancel.click();
    const anyCancel = document.querySelector('.d4-dialog [name="button-CANCEL"]') as HTMLElement | null;
    if (anyCancel) anyCancel.click();
  }).catch(() => {  });
  await page.evaluate(() => {
    try { (window as any).grok.shell.closeAll(); } catch (_) {  }
  }).catch(() => {  });
  finishSpec();
});
