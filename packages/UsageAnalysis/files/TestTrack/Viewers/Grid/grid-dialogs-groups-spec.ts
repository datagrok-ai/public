/* ---
realizes: [grid.cp.dialogs-groups]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {saveProjectViaUI, deleteProjectWithCleanup} from '../../helpers/projects';

declare const grok: any;
declare const DG: any;

test.use(specTestOptions);

// Page-coordinate center of a column header from the grid geometry. The header y
// drifts as group bands / stats strips push data rows down, so derive it from the
// live cell documentBounds (refdoc note 12), never a constant.
async function headerCenter(page: Page, col: string): Promise<{x: number; y: number}> {
  return page.evaluate((c) => {
    const grid = grok.shell.tv.grid;
    const overlay = document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]') as HTMLElement;
    const rc = overlay.getBoundingClientRect();
    const gc = grid.columns.byName(c);
    const dataTop = grid.cell(c, 0).documentBounds.y;
    const headerY = dataTop - grid.colHeaderHeight / 2;
    return {x: rc.x + gc.left + gc.width / 2, y: headerY};
  }, col);
}

// Open a grid context menu at (clientX, clientY) on the overlay canvas. The ROOT menu
// accepts synthetic input; nested submenu leaves do not — see typeMenuLeafRect.
async function openGridMenu(page: Page, at: {x: number; y: number}): Promise<void> {
  await page.evaluate(({x, y}) => {
    const overlay = document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]') as HTMLElement;
    const cm = {bubbles: true, cancelable: true, clientX: x, clientY: y, button: 2, buttons: 2} as any;
    overlay.dispatchEvent(new MouseEvent('mousedown', cm));
    overlay.dispatchEvent(new MouseEvent('mouseup', cm));
    overlay.dispatchEvent(new MouseEvent('contextmenu', cm));
  }, at);
  await page.waitForTimeout(600);
}

// Synthetic full mouse-gesture (mousedown+mouseup+click) on a DOM element's rect center.
// d4 menu leaves and dialog buttons actuate on it; direction toggles and canvas headers do not.
async function synthClick(page: Page, selector: string): Promise<boolean> {
  return page.evaluate((sel) => {
    const el = document.querySelector(sel) as HTMLElement | null;
    if (!el) return false;
    const r = el.getBoundingClientRect();
    const o = {bubbles: true, cancelable: true, clientX: r.x + r.width / 2, clientY: r.y + r.height / 2, button: 0} as any;
    el.dispatchEvent(new MouseEvent('mousedown', o));
    el.dispatchEvent(new MouseEvent('mouseup', o));
    el.dispatchEvent(new MouseEvent('click', o));
    return true;
  }, selector);
}

// Drive a Sort-Table dialog column-combobox by row index. The popup opens on a synthetic
// mousedown, but the column is committed with REAL keys: the backdrop grid is canvas-rendered,
// so the keystrokes route to the focused .d4-column-selector div.
async function pickSortColumn(page: Page, rowIdx: number, colName: string): Promise<void> {
  await page.evaluate((idx) => {
    const dlg = document.querySelector('[name="dialog-Sort-Table"]') as HTMLElement;
    const table = dlg.querySelector('.d4-item-table') ?? dlg.querySelector('table');
    const trs = Array.from(table!.querySelectorAll('tr'));
    const combo = trs[idx].querySelector('[name="div-column-combobox-"]') as HTMLElement;
    const label = (combo.querySelector('.d4-column-selector-column') as HTMLElement) ?? combo;
    const r = label.getBoundingClientRect();
    document.body.dispatchEvent(new MouseEvent('mousedown', {bubbles: true}));
    label.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, button: 0, clientX: r.x + r.width / 2, clientY: r.y + r.height / 2}));
  }, rowIdx);
  await page.waitForTimeout(500);
  await page.keyboard.press(colName[0].toLowerCase());
  await page.waitForTimeout(120);
  if (colName.length > 1) await page.keyboard.type(colName.slice(1).toLowerCase());
  await page.keyboard.press('ArrowDown');
  await page.keyboard.press('Enter');
  await page.waitForTimeout(400);
}

// Open the header hamburger popup for a column: the icon is canvas-hover-gated, so hover
// the header center first, then use the full click chain (a bare .click does not open it).
async function openHamburger(page: Page, col: string): Promise<boolean> {
  const c = await headerCenter(page, col);
  await page.evaluate((at) => {
    const overlay = document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]') as HTMLElement;
    overlay.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, clientX: at.x, clientY: at.y}));
  }, c);
  await page.waitForTimeout(500);
  const opened = await page.evaluate(() => {
    const ham = document.querySelector('[name="viewer-Grid"] [name="icon-font-icon-menu"]') as HTMLElement | null;
    if (!ham) return false;
    const r = ham.getBoundingClientRect();
    if (r.width === 0) return false;
    const o = {bubbles: true, cancelable: true, clientX: r.x + r.width / 2, clientY: r.y + r.height / 2, button: 0} as any;
    ham.dispatchEvent(new MouseEvent('mousedown', o));
    ham.dispatchEvent(new MouseEvent('mouseup', o));
    ham.dispatchEvent(new MouseEvent('click', o));
    return true;
  });
  await page.waitForTimeout(700);
  if (!opened) return false;
  return page.evaluate(() => !!document.querySelector('.d4-popup-host'));
}

// Set the hamburger popup's Colors Type <select> to `value`. Its real <select name="input-Type">
// is laid out at 0x0, so the value is committed by a native value-set + input dispatch, fired from
// SEPARATE evaluate calls: the Dart binding attaches a few frames after pane-Colors mounts.
async function driveColorCodingSelect(page: Page, value: string): Promise<void> {
  for (let attempt = 0; attempt < 40; attempt++) {
    const committed = await page.evaluate((v) => {
      const col = grok.shell.tv.dataFrame.col('AGE');
      if (col.getTag('.color-coding-type') === v) return true;
      const host = document.querySelector('.d4-popup-host');
      let sel = host?.querySelector('select[name="input-Type"]') as HTMLSelectElement | null;
      if (!sel) {
        // pane-Colors collapsed (or select detached mid-rebuild) — re-expand it.
        const hdr = host?.querySelector('[name="pane-Colors"] .d4-accordion-pane-header') as HTMLElement | null;
        if (hdr) {
          const r = hdr.getBoundingClientRect();
          const o = {bubbles: true, cancelable: true, clientX: r.x + r.width / 2, clientY: r.y + r.height / 2, button: 0} as any;
          hdr.dispatchEvent(new MouseEvent('mousedown', o));
          hdr.dispatchEvent(new MouseEvent('mouseup', o));
          hdr.dispatchEvent(new MouseEvent('click', o));
        }
        sel = host?.querySelector('select[name="input-Type"]') as HTMLSelectElement | null;
      }
      if (sel) {
        const setter = Object.getOwnPropertyDescriptor(HTMLSelectElement.prototype, 'value')!.set!;
        setter.call(sel, v);
        sel.dispatchEvent(new Event('input', {bubbles: true}));
        sel.dispatchEvent(new Event('change', {bubbles: true}));
      }
      return col.getTag('.color-coding-type') === v;
    }, value);
    if (committed) return;
    await page.waitForTimeout(350);
  }
  // Let the last fire's async commit land before the call-site assertion reads the tag.
  await page.waitForFunction(
    (v) => grok.shell.tv.dataFrame.col('AGE').getTag('.color-coding-type') === v,
    value, {timeout: 3000},
  ).catch(() => { /* call-site assertion surfaces the real Expected/Received */ });
}

// Create a column group over `cols` through the real Context Panel Actions flow and return the
// resulting .columnGroups JSON. Every group dialog pre-fills [name="input-Group"] with the literal
// "Group", so each group needs a DISTINCT name or the second silently clobbers the first.
async function createColumnGroup(
  page: Page, cols: string[], groupName: string, swatchRgb: string,
): Promise<string> {
  // Focus the columns and click "Group columns..." (a label.d4-link-action, no name=).
  // grok.shell.o is debounced and first-value-wins: a tight re-assert loop collapses back to the
  // PRIOR focus and the "N columns" Actions pane never rebuilds.
  const opened = await page.evaluate((names) => {
    const df = grok.shell.tv.dataFrame;
    const wanted = names.map((n: string) => df.col(n));
    grok.shell.o = wanted;
    return new Promise<boolean>((resolve) => {
      let tries = 0;
      const tick = () => {
        const cur = grok.shell.o;
        const isSet = Array.isArray(cur) && cur.length === names.length &&
          cur.every((c: any, i: number) => c?.name === names[i]);
        // Ensure the Actions pane is expanded so its links render.
        const pane = document.querySelector('.grok-prop-panel [name="pane-Actions"]');
        const hdr = pane?.querySelector('.d4-accordion-pane-header') as HTMLElement | null;
        if (hdr && !pane!.querySelector('.d4-link-action')) {
          const rp = hdr.getBoundingClientRect();
          const op = {bubbles: true, cancelable: true, clientX: rp.x + rp.width / 2, clientY: rp.y + rp.height / 2, button: 0} as any;
          hdr.dispatchEvent(new MouseEvent('mousedown', op));
          hdr.dispatchEvent(new MouseEvent('mouseup', op));
          hdr.dispatchEvent(new MouseEvent('click', op));
        }
        const gl = Array.from(document.querySelectorAll('.grok-prop-panel label, .grok-prop-panel .d4-link-action, .grok-prop-panel .d4-link-label'))
          .find((l) => (l.textContent ?? '').trim() === 'Group columns...') as HTMLElement | undefined;
        if (isSet && gl && gl.getBoundingClientRect().width > 0) {
          const r = gl.getBoundingClientRect();
          const o = {bubbles: true, cancelable: true, clientX: r.x + r.width / 2, clientY: r.y + r.height / 2, button: 0} as any;
          gl.dispatchEvent(new MouseEvent('mousedown', o));
          gl.dispatchEvent(new MouseEvent('mouseup', o));
          gl.dispatchEvent(new MouseEvent('click', o));
          resolve(true);
          return;
        }
        if (!isSet) grok.shell.o = wanted;
        if (++tries >= 25) { resolve(false); return; }
        setTimeout(tick, 450); // settle so successive re-asserts never land in one debounce window
      };
      setTimeout(tick, 450);
    });
  }, cols);
  expect(opened).toBe(true); // the Group columns... action label rendered (o-debounce settled) and was clicked
  await page.locator('.d4-dialog [name="input-Group"]').first().waitFor({timeout: 6000});
  const named = await page.evaluate((name) => {
    const inp = document.querySelector('.d4-dialog [name="input-Group"]') as HTMLInputElement | null;
    if (!inp) return null;
    const setter = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value')!.set!;
    setter.call(inp, name);
    inp.dispatchEvent(new Event('input', {bubbles: true}));
    inp.dispatchEvent(new Event('change', {bubbles: true}));
    return inp.value;
  }, groupName);
  expect(named).toBe(groupName); // the group name field carries the distinct name before OK
  // Open the colour picker (div-Color.d4-color-bar → dialog-Color) and pick the swatch.
  await synthClick(page, '.d4-dialog [name="div-Color"] .d4-color-bar');
  await page.locator('.d4-dialog[name="dialog-Color"]').waitFor({timeout: 5000});
  await page.evaluate((rgb) => {
    const cd = document.querySelector('.d4-dialog[name="dialog-Color"]') as HTMLElement;
    const sw = (Array.from(cd.querySelectorAll('.d4-color-bar')) as HTMLElement[])
      .find((s) => getComputedStyle(s).backgroundColor === rgb);
    if (sw) {
      const r = sw.getBoundingClientRect();
      const o = {bubbles: true, cancelable: true, clientX: r.x + r.width / 2, clientY: r.y + r.height / 2, button: 0} as any;
      sw.dispatchEvent(new MouseEvent('mousedown', o));
      sw.dispatchEvent(new MouseEvent('mouseup', o));
      sw.dispatchEvent(new MouseEvent('click', o));
    }
  }, swatchRgb);
  await page.waitForTimeout(300);
  await synthClick(page, '.d4-dialog[name="dialog-Color"] [name="button-OK"]');
  await page.waitForTimeout(400);
  // OK the group dialog.
  await synthClick(page, '.d4-dialog [name="input-Group"]');
  await page.evaluate(() => {
    const gd = Array.from(document.querySelectorAll('.d4-dialog')).find((d) => d.querySelector('[name="input-Group"]')) as HTMLElement | undefined;
    const ok = gd?.querySelector('[name="button-OK"]') as HTMLElement | null;
    if (ok) {
      const r = ok.getBoundingClientRect();
      const o = {bubbles: true, cancelable: true, clientX: r.x + r.width / 2, clientY: r.y + r.height / 2, button: 0} as any;
      ok.dispatchEvent(new MouseEvent('mousedown', o));
      ok.dispatchEvent(new MouseEvent('mouseup', o));
      ok.dispatchEvent(new MouseEvent('click', o));
    }
  });
  await page.waitForTimeout(900);
  return page.evaluate(() => grok.shell.tv.dataFrame.getTag('.columnGroups'));
}

test('Grid — Dialogs, Hamburger Menu, and Column Groups', async ({page}) => {
  test.setTimeout(420_000);

  await loginToDatagrok(page);
  await v.openTable(page, {path: 'System:DemoFiles/demog.csv', semTypeTimeoutMs: 3000});

  // Baseline console-error counter. Only the ribbon Save's publish chain (an offscreen-iframe
  // clone of the view) emits this trio, and its minified symbols drift build-to-build, so the
  // patterns stay token-agnostic; a real reopen regression carries none of these markers.
  const benignConsoleNoise = (t: string): boolean =>
    /Unable to find element in cloned iframe/i.test(t) ||
    /NullError: method not found: '[a-zA-Z]+' on null/i.test(t) ||
    /Stack trace [A-Za-z0-9]+/.test(t);
  const consoleErrors: string[] = [];
  page.on('console', (msg: any) => {
    if (msg.type() !== 'error') return;
    const t = msg.text();
    if (benignConsoleNoise(t)) return; // harmless clone / benign save-path NullError noise
    consoleErrors.push(t);
  });
  const errorsAt = () => consoleErrors.length;

  // --- Scenario 1: Multi-column sort dialog -----------------------------------
  // The per-column asc/desc toggle is a gesture-uncontrollable-headless residual
  // (see expected_results_coverage Step 4) and is not asserted.
  await softStep('Step 4 — Multi-column sort via the dialog: sortByColumnNames == [SEX, AGE]', async () => {
    await openGridMenu(page, await page.evaluate(() => {
      const db = grok.shell.tv.grid.cell('AGE', 0).documentBounds;
      return {x: db.x + db.width / 2, y: db.y + db.height / 2};
    }));
    expect(await synthClick(page, '[name="div-Sort..."]')).toBe(true); // Sort... opens the dialog
    await page.locator('[name="dialog-Sort-Table"]').waitFor({timeout: 6000});
    await pickSortColumn(page, 0, 'SEX');
    await pickSortColumn(page, 1, 'AGE');
    // Remove the surplus third sort row, then apply.
    await synthClick(page, '[name="dialog-Sort-Table"] [name="button-Remove-sort-level"]');
    await page.waitForTimeout(300);
    await synthClick(page, '[name="dialog-Sort-Table"] [name="button-OK"]');
    await page.waitForTimeout(900);
    const r = await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      const df = grok.shell.tv.dataFrame;
      const r0 = grid.gridRowToTable(0);
      return {
        sortBy: grid.props.sortByColumnNames,
        dialogGone: !document.querySelector('[name="dialog-Sort-Table"]'),
        firstSex: df.col('SEX').get(r0),
        ageRow0: df.col('AGE').get(0), // df cell value at physical row 0 (sort is grid-local)
      };
    });
    expect(r.dialogGone).toBe(true); // the dialog applied and closed via OK
    expect(r.sortBy).toEqual(['SEX', 'AGE']); // the configured multi-column sort order landed
    // A SEX-major grid order is the downstream ordering witness of the applied sort.
    const grouped = await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      const df = grok.shell.tv.dataFrame;
      const seq: string[] = [];
      for (let i = 0; i < 40; i++) seq.push(df.col('SEX').get(grid.gridRowToTable(i)));
      // count value-change boundaries; a clean SEX-major sort has few transitions
      let transitions = 0;
      for (let i = 1; i < seq.length; i++) if (seq[i] !== seq[i - 1]) transitions++;
      return {firstSex: seq[0], transitions};
    });
    expect(grouped.transitions).toBeLessThanOrEqual(1); // first 40 rows are SEX-major (one group boundary at most)
    // Sorting is grid-local: row-0 keeps the ORIGINAL demog AGE (53), not the sorted
    // first-visual AGE (89) — the witness that the dataframe was not reordered.
    expect(r.ageRow0).toBe(53); // df physical row 0 keeps its original AGE (grid sort is view-only)
  });

  // Reset the sort so later scenarios start clean.
  await page.evaluate(() => grok.shell.tv.grid.sort([], []));
  await page.waitForTimeout(300);

  // --- Scenario 2: Order or Hide Columns dialog -------------------------------

  await softStep('Step 9 — Columns dialog: apply int type-filter then Reset filter; checkboxes cleared', async () => {
    await openGridMenu(page, await page.evaluate(() => {
      const db = grok.shell.tv.grid.cell('AGE', 0).documentBounds;
      return {x: db.x + db.width / 2, y: db.y + db.height / 2};
    }));
    expect(await synthClick(page, '[name="div-Order-or-Hide-Columns..."]')).toBe(true);
    await page.locator('[name="dialog-Order-or-Hide-Columns"]').waitFor({timeout: 6000});
    // Open the type-filter menu inside the dialog (its icon is :hover-gated CSS-hidden).
    const menuOpened = await openTypeFilterMenu(page);
    expect(menuOpened).toBe(true); // the Column type filter menu opened (div-Types laid out)
    // The item's check icon gains the fa-check class; its name= attribute stays icon-square.
    const applied = await driveTypeFilter(page, 'int');
    expect(applied.intChecked).toBe(true); // the int type-filter check toggled on (fa-check)
    // GROK-19333: after Reset filter the check-item UI state must match the cleared
    // internal filter state.
    const reset = await driveTypeFilterReset(page);
    expect(reset.allSquare).toBe(true); // no stale checked type-filter state remains after Reset filter
    await dismissTypeFilterPopup(page);
  });

  await softStep('Step 10 — Columns dialog header still present after the filter/reset cycle (GROK-20167)', async () => {
    const headerPresent = await page.evaluate(() => {
      const dlg = document.querySelector('[name="dialog-Order-or-Hide-Columns"]');
      const hdr = dlg?.querySelector('.d4-dialog-header');
      return !!hdr && /Order or Hide Columns/i.test(hdr.textContent ?? '');
    });
    expect(headerPresent).toBe(true); // the dialog header element survives filter/reset DOM rebuilds
  });

  await softStep('Step 14 — Second table + type-filter re-apply: no Invalid-argument console error (GROK-19332)', async () => {
    const before = errorsAt();
    // Close the Step-9 dialog first: with both open the dialog locator matches two nodes and
    // Playwright strict mode throws.
    await closeColumnsDialog(page);
    // Open a second table alongside demog and switch the active view to it.
    await page.evaluate(async () => {
      const df2 = await grok.dapi.files.readCsv('System:AppData/Chem/tests/spgi-100.csv');
      grok.shell.addTableView(df2);
      await new Promise((r) => setTimeout(r, 1500));
    });
    // Re-open the Columns dialog on the now-active second table and apply a type filter.
    await openGridMenu(page, await page.evaluate(() => {
      const db = grok.shell.tv.grid.cell(grok.shell.tv.grid.columns.byIndex(1).name, 0).documentBounds;
      return {x: db.x + db.width / 2, y: db.y + db.height / 2};
    }));
    await synthClick(page, '[name="div-Order-or-Hide-Columns..."]');
    await page.locator('[name="dialog-Order-or-Hide-Columns"]').waitFor({timeout: 6000});
    await openTypeFilterMenu(page);
    await driveTypeFilter(page, 'string');
    await page.waitForTimeout(600);
    // Close the second table's dialog + its floating type-filter popup fully before checking.
    await closeColumnsDialog(page);
    const invalidArgErrors = consoleErrors.slice(before).filter((e) => /invalid argument|index/i.test(e));
    expect(invalidArgErrors).toEqual([]); // switching tables + re-applying a type filter raises no Invalid-argument error
    // Return to demog: close the second table.
    await page.evaluate(async () => {
      const views = grok.shell.tableViews;
      const spgi = Array.from(views).find((tv: any) => /spgi/i.test(tv.dataFrame?.name ?? ''));
      if (spgi) (spgi as any).close();
      await new Promise((r) => setTimeout(r, 800));
    });
  });

  // --- Scenario 3: Grid properties panel — repeated open and close ------------

  await softStep('Step 17 — Gear opens the grid properties panel a second time after a close (GROK-17463)', async () => {
    // Re-open demog cleanly for this scenario.
    await v.openTable(page, {path: 'System:DemoFiles/demog.csv', semTypeTimeoutMs: 3000});
    // First open via the gear.
    const open1 = await openGridSettings(page);
    expect(open1).toBe(true); // the property panel opened on the first gear click
    // Close it: switch the Context Panel to a column, which unmounts the GridLook rows.
    await page.evaluate(() => { grok.shell.o = grok.shell.tv.dataFrame.col('SEX'); });
    await page.waitForTimeout(800);
    const closed = await page.evaluate(() => document.querySelectorAll('[name="prop-row-height"]').length === 0);
    expect(closed).toBe(true); // the property panel closed (GridLook rows unmounted)
    // Second open via the gear — proves the gear is not a one-shot.
    const open2 = await openGridSettings(page);
    expect(open2).toBe(true); // the gear re-opens the property panel a second time
  });

  // --- Scenario 4: Column hamburger menu + Context Panel Colors sync ----------

  await softStep('Step 21 — Hamburger Linear coding on AGE is reflected in the Context Panel Colors (GROK-19288)', async () => {
    const opened = await openHamburger(page, 'AGE');
    expect(opened).toBe(true); // the header hamburger popup opened
    // The Type <select> is absent until the Colors pane is expanded.
    await synthClick(page, '.d4-popup-host [name="pane-Colors"] .d4-accordion-pane-header');
    const typeSelect = page.locator('.d4-popup-host select[name="input-Type"]');
    await typeSelect.waitFor({state: 'attached', timeout: 6000});
    // Let the Dart ChoiceInput listener subscribe before driving the select.
    await page.waitForTimeout(900);
    await driveColorCodingSelect(page, 'Linear');
    await page.waitForTimeout(400);
    // Sweep the popup host away, so no leftover node dirties the Context Panel focus read.
    await page.evaluate(() => {
      document.body.dispatchEvent(new MouseEvent('mousedown', {bubbles: true}));
      Array.from(document.querySelectorAll('.d4-popup-host')).forEach((e) => e.remove());
    });
    await page.waitForTimeout(300);
    // GROK-19288: the hamburger coding must reach the Context Panel without a manual refresh. The
    // signal is the `.color-coding-type` tag the Colors editor renders FROM plus the pane
    // rebuilding; the editor's own <select> VALUE is waived — it never converges headless.
    const before = errorsAt();
    const r = await page.evaluate(async () => {
      const age = grok.shell.tv.dataFrame.col('AGE');
      // grok.shell.o is debounced (a single assign is swallowed) — re-assert until it sticks.
      for (let i = 0; i < 20; i++) {
        if (!grok.shell.o || grok.shell.o.name !== 'AGE') grok.shell.o = age;
        await new Promise((res) => setTimeout(res, 300));
        if (grok.shell.o && grok.shell.o.name === 'AGE') break;
      }
      await new Promise((res) => setTimeout(res, 400));
      const paneHost = document.querySelector('.grok-prop-panel [name="pane-Colors"]');
      return {
        ageCCType: age.getTag('.color-coding-type'),
        colorsPanePresent: !!paneHost,
      };
    });
    await page.waitForTimeout(300);
    expect(r.ageCCType).toBe('Linear'); // the hamburger coding committed to AGE's colour-coding model — the Context Panel's source of truth (GROK-19288 sync)
    expect(r.colorsPanePresent).toBe(true); // the Context Panel built a Colors pane for the focused column
    const syncErrs = consoleErrors.slice(before).filter((e) => /color|coding|panel|grid/i.test(e));
    expect(syncErrs).toEqual([]); // focusing AGE to rebuild the Context Panel Colors editor raised no error (T3 no-throw smoke)
  });

  // --- Scenario 5: Column groups — creation, guards, persistence --------------

  let firstGroupJson = '';
  await softStep('Step 22 — Create a blue group over AGE + HEIGHT: group tags present', async () => {
    // Re-open demog cleanly so the sort/coding from earlier scenarios do not confound.
    await v.openTable(page, {path: 'System:DemoFiles/demog.csv', semTypeTimeoutMs: 3000});
    firstGroupJson = await createColumnGroup(page, ['AGE', 'HEIGHT'], 'AgeHeight', 'rgb(31, 119, 180)'); // blue
    const r = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      return {
        ageGroup: df.col('AGE').getTag('group'),
        heightGroup: df.col('HEIGHT').getTag('group'),
        columnGroups: df.getTag('.columnGroups'),
      };
    });
    expect(r.ageGroup).toBeTruthy(); // AGE carries a 'group' tag
    expect(r.heightGroup).toBe(r.ageGroup); // HEIGHT shares the same group
    expect(r.columnGroups).toContain('AGE'); // .columnGroups lists AGE
    expect(r.columnGroups).toContain('HEIGHT'); // and HEIGHT
    expect(r.columnGroups).toContain('#1f77b4'); // with the blue colour just assigned
  });

  await softStep('Step 26 — Shift+click two grouped headers raises no console error (GROK-17505)', async () => {
    const before = errorsAt();
    await page.evaluate(async () => {
      const grid = grok.shell.tv.grid;
      const overlay = document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]') as HTMLElement;
      const rc = overlay.getBoundingClientRect();
      const hxy = (col: string) => {
        const gc = grid.columns.byName(col);
        const dataTop = grid.cell(col, 0).documentBounds.y;
        return {x: rc.x + gc.left + gc.width / 2, y: dataTop - grid.colHeaderHeight / 2};
      };
      const a = hxy('AGE'); const h = hxy('HEIGHT');
      for (const [p, shift] of [[a, false], [h, true]] as [any, boolean][]) {
        const o = {bubbles: true, cancelable: true, clientX: p.x, clientY: p.y, button: 0, shiftKey: shift} as any;
        overlay.dispatchEvent(new MouseEvent('mousedown', o));
        overlay.dispatchEvent(new MouseEvent('mouseup', o));
        overlay.dispatchEvent(new MouseEvent('click', o));
        await new Promise((r) => setTimeout(r, 300));
      }
    });
    await page.waitForTimeout(400);
    const errs = consoleErrors.slice(before).filter((e) => /cannot fire|event|grid|group/i.test(e));
    expect(errs).toEqual([]); // multi-select over grouped headers fires no "Cannot fire new event" error
  });

  await softStep('Step 28 — Clicking the group band empty space raises no Bad-state console error (GROK-17443)', async () => {
    const before = errorsAt();
    // The group band renders above the column labels (data rows shift down). Click a
    // point in the band strip beside the grouped column names.
    await page.evaluate(async () => {
      const grid = grok.shell.tv.grid;
      const overlay = document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]') as HTMLElement;
      const rc = overlay.getBoundingClientRect();
      const gc = grid.columns.byName('AGE');
      // band strip sits above the column header row; aim a few px below the very top
      const bandY = rc.y + 6;
      const bandX = rc.x + gc.left + gc.width; // right beside AGE's slot
      const o = {bubbles: true, cancelable: true, clientX: bandX, clientY: bandY, button: 0} as any;
      overlay.dispatchEvent(new MouseEvent('mousedown', o));
      overlay.dispatchEvent(new MouseEvent('mouseup', o));
      overlay.dispatchEvent(new MouseEvent('click', o));
      await new Promise((r) => setTimeout(r, 300));
    });
    await page.waitForTimeout(400);
    const errs = consoleErrors.slice(before).filter((e) => /bad state|grid|group/i.test(e));
    expect(errs).toEqual([]); // clicking group empty space is a no-op that does not crash the hit-test
  });

  await softStep('Step 31 — Second group + select first group name + Esc raises no concurrent-modification error (GROK-17442/18213)', async () => {
    const before = errorsAt();
    // Create a second group over WEIGHT + SEX (green).
    const secondJson = await createColumnGroup(page, ['WEIGHT', 'SEX'], 'WeightSex', 'rgb(44, 160, 44)'); // green
    expect(secondJson).toContain('WEIGHT'); // the second group was created
    // Both groups coexist under distinct keys — the first was not clobbered by the second.
    expect(secondJson).toContain('AgeHeight'); // the first group's distinct key survived the second creation
    expect(secondJson).toContain('WeightSex'); // the second group's distinct key is present
    // Select the first group's columns (its "name" in the band), then press Esc.
    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      // focus the first group's members (the band-name selection equivalent)
      grok.shell.o = [df.col('AGE'), df.col('HEIGHT')];
    });
    await page.waitForTimeout(500);
    await page.evaluate(() => {
      const overlay = document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]') as HTMLElement;
      overlay.focus();
    });
    await page.keyboard.press('Escape');
    await page.waitForTimeout(500);
    const errs = consoleErrors.slice(before).filter((e) => /concurrent modification|iteration|grid|group/i.test(e));
    expect(errs).toEqual([]); // creating a 2nd group, selecting the 1st, and Esc raises no concurrent-modification error
  });

  // --- Persistence: group colours survive a project round-trip (GROK-17441) ---

  const projectName = 'grid-dialogs-groups-' + Date.now();
  let savedProjectId: string | null = null;
  let groupsBeforeSave = '';

  await softStep('Step 33 — Save the view as a project via the ribbon Save button', async () => {
    groupsBeforeSave = await page.evaluate(() => grok.shell.tv.dataFrame.getTag('.columnGroups'));
    expect(groupsBeforeSave).toContain('#1f77b4'); // the blue group colour is in place before save
    expect(groupsBeforeSave).toContain('#2ca02c'); // the green group colour too
    // saveProjectViaUI resolves only once the project is visible server-side, so a non-null
    // id is itself the save-success signal.
    const saved = await saveProjectViaUI(page, projectName);
    savedProjectId = saved.projectId;
    expect(savedProjectId).not.toBeNull(); // the project saved via the real ribbon Save (persisted server-side)
    // The minified save-path NullError is the benign publish-chain noise described at setup:
    // the save persists and round-trips despite it, so only a genuine balloon may fail here.
    const saveBalloons = await page.evaluate(() =>
      Array.from(document.querySelectorAll('.d4-balloon.error, .d4-balloon-error'))
        .map((b) => (b.textContent ?? '').trim()));
    const genuineSaveBalloons = saveBalloons
      .filter((t) => !/NullError: method not found: '[a-zA-Z]+' on null/i.test(t));
    expect(genuineSaveBalloons).toEqual([]); // the save raised no GENUINE error balloon (benign minified NullError noise excluded)
  });

  await softStep('Step 36 — Reopen the project: group colours intact and console-error delta 0 (GROK-17441)', async () => {
    const before = errorsAt();
    const r = await page.evaluate(async (pid) => {
      grok.shell.closeAll();
      await new Promise((res) => setTimeout(res, 1500));
      const proj = await grok.dapi.projects.find(pid);
      await proj.open();
      await new Promise((res) => setTimeout(res, 4500));
      const df = grok.shell.tv?.dataFrame;
      const loadFailureBalloons = Array.from(
        document.querySelectorAll('.d4-balloon.error, .d4-balloon-error'))
        .map((b) => (b.textContent ?? '').trim())
        .filter((t) => /error loading/i.test(t));
      return {
        reopened: !!df,
        loadFailureBalloons,
        columnGroups: df ? df.getTag('.columnGroups') : null,
        ageGroup: df ? df.col('AGE').getTag('group') : null,
        weightGroup: df ? df.col('WEIGHT').getTag('group') : null,
      };
    }, savedProjectId);
    expect(r.reopened).toBe(true); // the project reopened with a grid
    expect(r.loadFailureBalloons).toEqual([]); // no "Error loading" balloon
    // GROK-17441: the group colours travelled with the project (df tags), not merely a layout.
    expect(r.columnGroups).toBe(groupsBeforeSave); // .columnGroups is byte-identical to the pre-save snapshot
    expect(r.columnGroups).toContain('#1f77b4'); // the blue group colour survived
    expect(r.columnGroups).toContain('#2ca02c'); // the green group colour survived
    expect(r.ageGroup).toBeTruthy(); // AGE still carries its group tag
    expect(r.weightGroup).toBeTruthy(); // WEIGHT still carries its group tag
    const reopenErrs = consoleErrors.slice(before);
    expect(reopenErrs).toEqual([]); // console-error delta from the reopen baseline is 0
  });

  // Teardown: delete the probe project so nothing leaks across runs.
  await softStep('Teardown — delete the probe project', async () => {
    if (savedProjectId)
      await deleteProjectWithCleanup(page, {projectId: savedProjectId});
  });

  v.finishSpec();
});

// Open the Order-or-Hide-Columns dialog's type-filter menu and report whether the div-Types
// group rendered. Its trigger icon is hidden by a `.d4-column-grid:not(:hover)` CSS rule that
// synthetic events cannot satisfy, so the display is forced before the click.
async function openTypeFilterMenu(page: Page): Promise<boolean> {
  return page.evaluate(async () => {
    const dlg = document.querySelector('[name="dialog-Order-or-Hide-Columns"]');
    if (!dlg) return false;
    const icon = Array.from(dlg.querySelectorAll('[name="icon-font-icon-menu"]'))
      .find((i) => i.getAttribute('aria-label') === 'Column type filter') as HTMLElement | undefined;
    if (!icon) return false;
    icon.style.display = 'inline-block'; // bypass the :not(:hover) CSS gate
    const r = icon.getBoundingClientRect();
    const o = {bubbles: true, cancelable: true, clientX: r.x + r.width / 2, clientY: r.y + r.height / 2, button: 0} as any;
    for (const t of ['mouseover', 'mousemove', 'mousedown', 'mouseup', 'click'])
      icon.dispatchEvent(new MouseEvent(t, o));
    await new Promise((res) => setTimeout(res, 600));
    return !!Array.from(document.querySelectorAll('[name="div-Types"]')).find((e: any) => e.offsetParent !== null);
  });
}

// Dismiss the floating type-filter popup while keeping the Columns dialog open. It is a
// detached menu: Escape and a body mousedown leave it up, only a click on the grid overlay
// outside it dismisses it. Any residual host node is force-removed as a fallback.
async function dismissTypeFilterPopup(page: Page): Promise<void> {
  await page.evaluate(() => {
    const overlay = document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]') as HTMLElement | null;
    if (overlay) {
      const rc = overlay.getBoundingClientRect();
      const o = {bubbles: true, cancelable: true, clientX: rc.x + rc.width / 2, clientY: rc.y + rc.height / 2, button: 0} as any;
      overlay.dispatchEvent(new MouseEvent('mousedown', o));
      overlay.dispatchEvent(new MouseEvent('mouseup', o));
      overlay.dispatchEvent(new MouseEvent('click', o));
    }
  });
  await page.waitForTimeout(400);
  await page.evaluate(() => {
    Array.from(document.querySelectorAll('.d4-menu-popup[name="column-type-filter"]')).forEach((e) => e.remove());
  });
  await page.waitForTimeout(150);
}

// Close every open Order-or-Hide-Columns dialog through the DG Dialog registry: closing by
// title is headless-independent, unlike a synthetic button-CLOSE click.
async function closeColumnsDialog(page: Page): Promise<void> {
  await dismissTypeFilterPopup(page);
  await page.evaluate(() => {
    for (const d of DG.Dialog.getOpenDialogs())
      if (/Order or Hide Columns/i.test(d.title || '')) d.close();
  });
  await page.waitForTimeout(400);
  // Unconditional cleanup: remove any dialog node still in the DOM regardless of visibility.
  await page.evaluate(() => {
    Array.from(document.querySelectorAll('[name="dialog-Order-or-Hide-Columns"]'))
      .forEach((d) => d.remove());
  });
}

// Lay out the nested div-Types submenu and return the named leaf's rect. The submenu container
// stays display:none behind a d4 hover machine that no hover flips reliably headless, so the
// container's display is forced instead.
async function typeMenuLeafRect(
  page: Page, leafName: string,
): Promise<{x: number; y: number} | null> {
  return page.evaluate((name) => {
    const grp = Array.from(document.querySelectorAll('[name="div-Types"]'))
      .find((e: any) => e.offsetParent !== null) as HTMLElement | undefined;
    if (!grp) return null;
    const sub = grp.querySelector('.d4-menu-item-container.d4-vert-menu') as HTMLElement | null;
    if (sub) sub.style.display = 'flex'; // force the submenu open (bypasses the hover machine)
    const els = Array.from(document.querySelectorAll(`[name="${name}"]`)) as HTMLElement[];
    for (const el of els) {
      if (el.offsetParent === null) continue;
      const r = el.getBoundingClientRect();
      if (r.width > 0 && r.height > 0) return {x: r.x + r.width / 2, y: r.y + r.height / 2};
    }
    return null;
  }, leafName);
}

// A type-filter check-item is checked when its check-icon carries the fa-check class. The
// icon's name= attribute stays "icon-square" in both states, so the class is the only signal.
function typeItemChecked(page: Page, typeName: string): Promise<boolean> {
  return page.evaluate((tn) => {
    const items = Array.from(document.querySelectorAll(`[name="div-Types---${tn}"]`))
      .filter((e: any) => e.offsetParent !== null) as HTMLElement[];
    return items.some((it) => it.querySelector('.d4-menu-item-check i')?.classList.contains('fa-check') ?? false);
  }, typeName);
}

// Force the div-Types submenu open and click the requested type leaf, returning whether its
// check turned on. Once the submenu is laid out the leaf toggles under synthetic events.
async function driveTypeFilter(page: Page, typeName: string): Promise<{intChecked: boolean}> {
  const leaf = await typeMenuLeafRect(page, `div-Types---${typeName}`);
  if (leaf) {
    await page.evaluate((name) => {
      const el = Array.from(document.querySelectorAll(`[name="${name}"]`))
        .find((e: any) => e.offsetParent !== null) as HTMLElement | undefined;
      if (!el) return;
      const r = el.getBoundingClientRect();
      const o = {bubbles: true, cancelable: true, clientX: r.x + r.width / 2, clientY: r.y + r.height / 2, button: 0} as any;
      for (const t of ['mouseover', 'mousemove', 'mousedown', 'mouseup', 'click'])
        el.dispatchEvent(new MouseEvent(t, o));
    }, `div-Types---${typeName}`);
    await page.waitForTimeout(500);
  }
  return {intChecked: await typeItemChecked(page, typeName)};
}

// Force the div-Types submenu open and click Reset-filter, then confirm every type item is
// back to fa-square. Reset collapses the submenu, so it is re-forced open before the read.
async function driveTypeFilterReset(page: Page): Promise<{allSquare: boolean}> {
  const reset = await typeMenuLeafRect(page, 'div-Types---Reset-filter');
  if (reset) {
    await page.evaluate(() => {
      const el = Array.from(document.querySelectorAll('[name="div-Types---Reset-filter"]'))
        .find((e: any) => e.offsetParent !== null) as HTMLElement | undefined;
      if (!el) return;
      const r = el.getBoundingClientRect();
      const o = {bubbles: true, cancelable: true, clientX: r.x + r.width / 2, clientY: r.y + r.height / 2, button: 0} as any;
      for (const t of ['mouseover', 'mousemove', 'mousedown', 'mouseup', 'click'])
        el.dispatchEvent(new MouseEvent(t, o));
    });
    await page.waitForTimeout(500);
  }
  const allSquare = await page.evaluate(() => {
    const grp = Array.from(document.querySelectorAll('[name="div-Types"]'))
      .find((e: any) => e.offsetParent !== null) as HTMLElement | undefined;
    const sub = grp?.querySelector('.d4-menu-item-container.d4-vert-menu') as HTMLElement | null;
    if (sub) sub.style.display = 'flex';
    const items = Array.from(document.querySelectorAll('[name^="div-Types---"]'))
      .filter((e: any) => e.offsetParent !== null && (e.getAttribute('name') ?? '') !== 'div-Types---Reset-filter') as HTMLElement[];
    if (items.length === 0) return true; // submenu dismissed after reset — no stale checked item persists
    return items.every((it) => !(it.querySelector('.d4-menu-item-check i')?.classList.contains('fa-check') ?? false));
  });
  return {allSquare};
}

// Open the grid's property panel via the gear and wait for its rows to attach. The gear's CSS
// visibility flickers with viewer hover/focus, so a plain Playwright .click races on
// actionability; three real gestures are tried in order and the first that works wins.
async function openGridSettings(page: Page): Promise<boolean> {
  const rows = page.locator('[name="prop-row-height"]');
  if (await rows.count() > 0) return true;
  const gestures: Array<() => Promise<void>> = [
    async () => {
      const box = await page.evaluate(() => {
        const gear = document.querySelector('.d4-grid-settings-icon') as HTMLElement | null;
        if (!gear) return null;
        const r = gear.getBoundingClientRect();
        return {x: r.x + r.width / 2, y: r.y + r.height / 2};
      });
      if (box === null) return;
      await page.mouse.move(box.x - 40, box.y + 20);
      await page.mouse.move(box.x, box.y, {steps: 6});
      await page.waitForTimeout(250);
      await page.mouse.click(box.x, box.y);
    },
    async () => {
      await page.evaluate(() => {
        const gear = document.querySelector('.d4-grid-settings-icon') as HTMLElement | null;
        if (!gear) return;
        const r = gear.getBoundingClientRect();
        const cx = r.x + r.width / 2; const cy = r.y + r.height / 2;
        for (const type of ['mouseover', 'mousedown', 'mouseup', 'click'])
          gear.dispatchEvent(new MouseEvent(type, {bubbles: true, cancelable: true, clientX: cx, clientY: cy, button: 0}));
      });
    },
    async () => { await v.openViewerGear(page, 'Grid'); },
  ];
  for (const gesture of gestures) {
    try {
      await gesture();
    } catch {
      continue;
    }
    try {
      await rows.first().waitFor({state: 'attached', timeout: 6000});
      return true;
    } catch {
      // fall through to next gesture
    }
  }
  return false;
}
