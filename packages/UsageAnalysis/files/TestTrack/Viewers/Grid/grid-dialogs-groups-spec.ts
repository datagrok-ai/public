/* ---
realizes: [grid.cp.dialogs-groups]
--- */
import {expect, Page} from '@playwright/test';
import {localTest as test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import * as g from './grid-helpers';

declare const grok: any;
declare const DG: any;

// Steps 33 and 36 (the project round-trip of the column groups) run on the server lane in
// grid-server-spec.ts.
test.use(specTestOptions);

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

// The dialog opens with three fixed levels (the first column plus two empty ones); a pick never
// appends a level, and the single remove button sits on the last row.
function sortDialogPick(page: Page, rowIdx: number, colName: string): Promise<boolean> {
  return g.pickColumnInCombo(page, `[name="dialog-Sort-Table"] tr:nth-child(${rowIdx + 1}) [name^="div-column-combobox"]`, colName);
}

async function openHamburger(page: Page, col: string): Promise<boolean> {
  const c = await g.headerCenter(page, col);
  await page.evaluate(({at, sel}) => {
    const overlay = document.querySelector(sel) as HTMLElement;
    overlay.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, clientX: at.x, clientY: at.y}));
  }, {at: c, sel: g.OVERLAY});

  const opened = await v.pollValue(() => page.evaluate(() => {
    const ham = document.querySelector('[name="viewer-Grid"] [name="icon-font-icon-menu"]') as HTMLElement | null;
    if (!ham) return false;
    const r = ham.getBoundingClientRect();
    if (r.width === 0) return false;
    const o = {bubbles: true, cancelable: true, clientX: r.x + r.width / 2, clientY: r.y + r.height / 2, button: 0} as any;
    ham.dispatchEvent(new MouseEvent('mousedown', o));
    ham.dispatchEvent(new MouseEvent('mouseup', o));
    ham.dispatchEvent(new MouseEvent('click', o));
    return true;
  }), (clicked) => clicked, 500, 50);
  if (!opened) return false;
  return v.pollValue(() => page.evaluate(() => !!document.querySelector('.d4-popup-host')), (p) => p, 700, 50);
}

async function driveColorCodingSelect(page: Page, value: string): Promise<void> {
  for (let attempt = 0; attempt < 40; attempt++) {
    const committed = await page.evaluate((v) => {
      const col = grok.shell.tv.dataFrame.col('AGE');
      if (col.getTag('.color-coding-type') === v) return true;
      const host = document.querySelector('.d4-popup-host');
      let sel = host?.querySelector('select[name="input-Type"]') as HTMLSelectElement | null;
      if (!sel) {
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
    await page.waitForTimeout(150);
  }
}

async function createColumnGroup(
  page: Page, cols: string[], groupName: string, swatchRgb: string,
): Promise<string> {
  const opened = await page.evaluate(async (names) => {
    const w = window as any;
    const df = grok.shell.tv.dataFrame;
    const wanted = names.map((n: string) => df.col(n));
    const click = (el: HTMLElement) => {
      const r = el.getBoundingClientRect();
      const o = {bubbles: true, cancelable: true, clientX: r.x + r.width / 2, clientY: r.y + r.height / 2, button: 0} as any;
      el.dispatchEvent(new MouseEvent('mousedown', o));
      el.dispatchEvent(new MouseEvent('mouseup', o));
      el.dispatchEvent(new MouseEvent('click', o));
    };
    grok.shell.o = wanted;
    let lastExpand = 0;
    const link = await w.__poll(() => {
      const cur = grok.shell.o;
      const isSet = Array.isArray(cur) && cur.length === names.length && cur.every((c: any, i: number) => c?.name === names[i]);
      if (!isSet) { grok.shell.o = wanted; return null; }
      const gl = Array.from(document.querySelectorAll('.grok-prop-panel label, .grok-prop-panel .d4-link-action, .grok-prop-panel .d4-link-label'))
        .find((l) => (l.textContent ?? '').trim() === 'Group columns...') as HTMLElement | undefined;
      if (gl && gl.getBoundingClientRect().width > 0) return gl;
      // the Actions pane may come collapsed with its links already rendered; expand it, once per settle
      const hdr = document.querySelector('.grok-prop-panel [name="pane-Actions"] .d4-accordion-pane-header') as HTMLElement | null;
      if (hdr && Date.now() - lastExpand > 1500) { lastExpand = Date.now(); click(hdr); }
      return null;
    }, (e: any) => !!e, 11_000, 450);
    if (!link) {
      const cur = grok.shell.o;
      return 'no Group columns link; shell.o=' + (Array.isArray(cur) ? cur.map((c: any) => c?.name).join('+') : String(cur?.name ?? cur)) +
        '; panes=' + Array.from(document.querySelectorAll('.grok-prop-panel [name^="pane-"]')).map((p) => p.getAttribute('name') +
          ((p as HTMLElement).offsetParent === null ? '(hidden)' : '')).join(',') +
        '; actions=' + Array.from(document.querySelectorAll('.grok-prop-panel .d4-link-action')).map((l) => (l.textContent ?? '').trim()).join('|');
    }
    click(link);
    return true;
  }, cols);
  expect(opened, String(opened)).toBe(true);
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
  expect(named).toBe(groupName);

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
  await page.locator('.d4-dialog[name="dialog-Color"]').waitFor({state: 'detached', timeout: 2000}).catch(() => {});

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
  return v.pollValue(() => page.evaluate(() => grok.shell.tv.dataFrame.getTag('.columnGroups') ?? ''),
    (t: string) => t.includes(groupName), 1500, 50);
}

test('Grid — Dialogs, Hamburger Menu, and Column Groups', async ({page}) => {
  test.setTimeout(300_000);

  await openDatagrok(page);
  const flags = await g.readShellFlags(page);
  const errors = g.trackErrors(page, g.BENIGN_NOISE);
  try {
    await v.openTable(page, {path: g.DEMOG, semTypeTimeoutMs: 3000});

    await softStep('Step 4 — Multi-column sort via the dialog: sortByColumnNames == [SEX, AGE]', async () => {
      expect(await g.clickMenuLeaf(page, await g.cellCenter(page, 'AGE', 0), [], 'div-Sort...')).toBe(true);
      await page.locator('[name="dialog-Sort-Table"]').waitFor({timeout: 6000});
      expect(await sortDialogPick(page, 0, 'SEX')).toBe(true);
      expect(await sortDialogPick(page, 1, 'AGE')).toBe(true);

      const levelsBefore = await page.evaluate(() => document.querySelectorAll('[name="dialog-Sort-Table"] tr').length);
      await synthClick(page, '[name="dialog-Sort-Table"] [name="button-Remove-sort-level"]');
      const levels = await v.pollValue(() => page.evaluate(() => ({
        rows: document.querySelectorAll('[name="dialog-Sort-Table"] tr').length,
        labels: Array.from(document.querySelectorAll('[name="dialog-Sort-Table"] .d4-column-selector-column')).map((e) => e.textContent?.trim()),
      })), (x) => x.rows === levelsBefore - 1, 1000, 50);
      expect(levels.labels).toEqual(['SEX', 'AGE']);

      await page.evaluate(() => new Promise<void>((resolve) => {
        const grid = grok.shell.tv.grid;
        const sub = grid.onRowsSorted.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(() => { sub.unsubscribe(); resolve(); }, 1500);
        const ok = document.querySelector('[name="dialog-Sort-Table"] [name="button-OK"]') as HTMLElement;
        const r = ok.getBoundingClientRect();
        const o = {bubbles: true, cancelable: true, clientX: r.x + r.width / 2, clientY: r.y + r.height / 2, button: 0} as any;
        for (const t of ['mousedown', 'mouseup', 'click']) ok.dispatchEvent(new MouseEvent(t, o));
      }));
      const r = await v.pollValue(() => page.evaluate(() => {
        const grid = grok.shell.tv.grid;
        const df = grok.shell.tv.dataFrame;
        const r0 = grid.gridRowToTable(0);
        return {
          sortBy: grid.props.sortByColumnNames,
          dialogGone: !document.querySelector('[name="dialog-Sort-Table"]'),
          firstSex: df.col('SEX').get(r0),
          ageRow0: df.col('AGE').get(0),
        };
      }), (x) => x.dialogGone && JSON.stringify(x.sortBy) === JSON.stringify(['SEX', 'AGE']), 1500, 50);
      expect(r.dialogGone).toBe(true);
      expect(r.sortBy).toEqual(['SEX', 'AGE']);

      const grouped = await page.evaluate(() => {
        const grid = grok.shell.tv.grid;
        const df = grok.shell.tv.dataFrame;
        const seq: string[] = [];
        for (let i = 0; i < 40; i++) seq.push(df.col('SEX').get(grid.gridRowToTable(i)));
        let transitions = 0;
        for (let i = 1; i < seq.length; i++) if (seq[i] !== seq[i - 1]) transitions++;
        return {firstSex: seq[0], transitions};
      });
      expect(grouped.transitions).toBeLessThanOrEqual(1);
      expect(r.ageRow0).toBe(53);
    });

    await page.evaluate(() => grok.shell.tv.grid.sort([], []));
    await v.waitForViewerRendered(page, 'Grid', 300);

    await softStep('Step 9 — Columns dialog: apply int type-filter then Reset filter; checkboxes cleared', async () => {
      expect(await g.clickMenuLeaf(page, await g.cellCenter(page, 'AGE', 0), [], 'div-Order-or-Hide-Columns...')).toBe(true);
      await page.locator('[name="dialog-Order-or-Hide-Columns"]').waitFor({timeout: 6000});

      expect(await openTypeFilterMenu(page)).toBe(true);
      const applied = await driveTypeFilter(page, 'int');
      expect(applied.intChecked).toBe(true);
      const reset = await driveTypeFilterReset(page);
      expect(reset.allSquare).toBe(true);
      await dismissTypeFilterPopup(page);
    });

    await softStep('Step 10 — Columns dialog header still present after the filter/reset cycle (GROK-20167)', async () => {
      const headerPresent = await page.evaluate(() => {
        const dlg = document.querySelector('[name="dialog-Order-or-Hide-Columns"]');
        const hdr = dlg?.querySelector('.d4-dialog-header');
        return !!hdr && /Order or Hide Columns/i.test(hdr.textContent ?? '');
      });
      expect(headerPresent).toBe(true);
    });

    await softStep('Step 14 — Second table + type-filter re-apply: no Invalid-argument console error (GROK-19332)', async () => {
      const before = errors.count();
      await closeColumnsDialog(page);

      await page.evaluate(async () => {
        const w = window as any;
        const df2 = await w.__readCsv('System:AppData/Chem/tests/spgi-100.csv');
        grok.shell.addTableView(df2);
        await w.__tableReady(3000);
      });
      await page.locator('.d4-grid[name="viewer-Grid"]').first().waitFor({timeout: 10000});

      const at = await page.evaluate(() => {
        const db = grok.shell.tv.grid.cell(grok.shell.tv.grid.columns.byIndex(1).name, 0).documentBounds;
        return {x: db.x + db.width / 2, y: db.y + db.height / 2};
      });
      expect(await g.clickMenuLeaf(page, at, [], 'div-Order-or-Hide-Columns...')).toBe(true);
      await page.locator('[name="dialog-Order-or-Hide-Columns"]').waitFor({timeout: 6000});
      await openTypeFilterMenu(page);
      await driveTypeFilter(page, 'string');
      // the error window after the type filter is the assertion
      await page.waitForTimeout(600);

      await closeColumnsDialog(page);
      expect(errors.list.slice(before)).toEqual([]);

      await page.evaluate(async () => {
        const w = window as any;
        const spgi = Array.from(grok.shell.tableViews).find((tv: any) => /spgi/i.test(tv.dataFrame?.name ?? ''));
        if (spgi) (spgi as any).close();
        await w.__poll(() => Array.from(grok.shell.tableViews).length, (n: number) => n === 1, 800, 50);
      });
    });

    await softStep('Step 17 — Gear opens the grid properties panel a second time after a close (GROK-17463)', async () => {
      await v.openTable(page, {path: g.DEMOG, semTypeTimeoutMs: 3000});

      expect(await g.openGridSettings(page)).toBe(true);
      await page.evaluate(() => { grok.shell.o = grok.shell.tv.dataFrame.col('SEX'); });
      const closed = await v.pollValue(
        () => page.evaluate(() => document.querySelectorAll('[name="prop-row-height"]').length === 0),
        (c) => c, 800, 50);
      expect(closed).toBe(true);
      expect(await g.openGridSettings(page)).toBe(true);
    });

    await softStep('Step 21 — Hamburger Linear coding on AGE is reflected in the Context Panel Colors (GROK-19288)', async () => {
      expect(await openHamburger(page, 'AGE')).toBe(true);

      await synthClick(page, '.d4-popup-host [name="pane-Colors"] .d4-accordion-pane-header');
      await page.locator('.d4-popup-host select[name="input-Type"]').waitFor({state: 'attached', timeout: 6000});
      await driveColorCodingSelect(page, 'Linear');

      await page.evaluate(() => {
        document.body.dispatchEvent(new MouseEvent('mousedown', {bubbles: true}));
        Array.from(document.querySelectorAll('.d4-popup-host')).forEach((e) => e.remove());
      });

      const before = errors.count();
      const r = await page.evaluate(async () => {
        const w = window as any;
        const age = grok.shell.tv.dataFrame.col('AGE');
        grok.shell.o = age;
        const paneHost = await w.__poll(() => {
          if (!grok.shell.o || grok.shell.o.name !== 'AGE') grok.shell.o = age;
          return document.querySelector('.grok-prop-panel [name="pane-Colors"]');
        }, (e: any) => !!e, 6000, 100);
        return {ageCCType: age.getTag('.color-coding-type'), colorsPanePresent: !!paneHost};
      });
      // the error window after the panel sync is the assertion
      await page.waitForTimeout(300);
      expect(r.ageCCType).toBe('Linear');
      expect(r.colorsPanePresent).toBe(true);
      const syncErrs = errors.list.slice(before).filter((e) => /color|coding|panel|grid/i.test(e));
      expect(syncErrs).toEqual([]);
    });

    await softStep('Step 22 — Create a blue group over AGE + HEIGHT: group tags present', async () => {
      await v.openTable(page, {path: g.DEMOG, semTypeTimeoutMs: 3000});
      await createColumnGroup(page, ['AGE', 'HEIGHT'], 'AgeHeight', 'rgb(31, 119, 180)');
      const r = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        return {
          ageGroup: df.col('AGE').getTag('group'),
          heightGroup: df.col('HEIGHT').getTag('group'),
          columnGroups: df.getTag('.columnGroups'),
        };
      });
      expect(r.ageGroup).toBeTruthy();
      expect(r.heightGroup).toBe(r.ageGroup);
      expect(r.columnGroups).toContain('AGE');
      expect(r.columnGroups).toContain('HEIGHT');
      expect(r.columnGroups).toContain('#1f77b4');
    });

    await softStep('Step 26 — Shift+click two grouped headers raises no console error (GROK-17505)', async () => {
      const before = errors.count();
      await page.evaluate(async (sel) => {
        const grid = grok.shell.tv.grid;
        const overlay = document.querySelector(sel) as HTMLElement;
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
      }, g.OVERLAY);
      // the error window after the clicks is the assertion
      await page.waitForTimeout(400);
      expect(errors.list.slice(before)).toEqual([]);
    });

    await softStep('Step 28 — Clicking the group band empty space raises no Bad-state console error (GROK-17443)', async () => {
      const before = errors.count();
      await page.evaluate(async (sel) => {
        const grid = grok.shell.tv.grid;
        const overlay = document.querySelector(sel) as HTMLElement;
        const rc = overlay.getBoundingClientRect();
        const gc = grid.columns.byName('AGE');
        const bandY = rc.y + 6;
        const bandX = rc.x + gc.left + gc.width;
        const o = {bubbles: true, cancelable: true, clientX: bandX, clientY: bandY, button: 0} as any;
        overlay.dispatchEvent(new MouseEvent('mousedown', o));
        overlay.dispatchEvent(new MouseEvent('mouseup', o));
        overlay.dispatchEvent(new MouseEvent('click', o));
        await new Promise((r) => setTimeout(r, 300));
      }, g.OVERLAY);
      // the error window after the click is the assertion
      await page.waitForTimeout(400);
      expect(errors.list.slice(before)).toEqual([]);
    });

    await softStep('Step 31 — Second group + select first group name + Esc raises no concurrent-modification error (GROK-17442/18213)', async () => {
      const before = errors.count();

      const secondJson = await createColumnGroup(page, ['WEIGHT', 'SEX'], 'WeightSex', 'rgb(44, 160, 44)');
      expect(secondJson).toContain('WEIGHT');
      expect(secondJson).toContain('AgeHeight');
      expect(secondJson).toContain('WeightSex');

      await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        grok.shell.o = [df.col('AGE'), df.col('HEIGHT')];
      });
      await page.waitForTimeout(500);
      await g.focusGrid(page);
      await page.keyboard.press('Escape');
      // the error window after Esc is the assertion
      await page.waitForTimeout(500);
      expect(errors.list.slice(before)).toEqual([]);
    });
  } finally {
    errors.stop();
    await g.leaveShellClean(page, flags);
  }
  v.finishSpec();
});

async function openTypeFilterMenu(page: Page): Promise<boolean> {
  const clicked = await page.evaluate(() => {
    const dlg = document.querySelector('[name="dialog-Order-or-Hide-Columns"]');
    if (!dlg) return false;
    const icon = Array.from(dlg.querySelectorAll('[name="icon-font-icon-menu"]'))
      .find((i) => i.getAttribute('aria-label') === 'Column type filter') as HTMLElement | undefined;
    if (!icon) return false;
    icon.style.display = 'inline-block';
    const r = icon.getBoundingClientRect();
    const o = {bubbles: true, cancelable: true, clientX: r.x + r.width / 2, clientY: r.y + r.height / 2, button: 0} as any;
    for (const t of ['mouseover', 'mousemove', 'mousedown', 'mouseup', 'click'])
      icon.dispatchEvent(new MouseEvent(t, o));
    return true;
  });
  if (!clicked) return false;
  return v.pollValue(() => page.evaluate(() =>
    !!Array.from(document.querySelectorAll('[name="div-Types"]')).find((e: any) => e.offsetParent !== null)),
  (open) => open, 600, 50);
}

async function dismissTypeFilterPopup(page: Page): Promise<void> {
  await page.evaluate((sel) => {
    const overlay = document.querySelector(sel) as HTMLElement | null;
    if (overlay) {
      const rc = overlay.getBoundingClientRect();
      const o = {bubbles: true, cancelable: true, clientX: rc.x + rc.width / 2, clientY: rc.y + rc.height / 2, button: 0} as any;
      overlay.dispatchEvent(new MouseEvent('mousedown', o));
      overlay.dispatchEvent(new MouseEvent('mouseup', o));
      overlay.dispatchEvent(new MouseEvent('click', o));
    }
  }, g.OVERLAY);
  await v.pollValue(() => page.evaluate(() =>
    !Array.from(document.querySelectorAll('[name="div-Types"]')).find((e: any) => e.offsetParent !== null)),
  (gone) => gone, 400, 50);
  await page.evaluate(() => {
    Array.from(document.querySelectorAll('.d4-menu-popup[name="column-type-filter"]')).forEach((e) => e.remove());
  });
}

async function closeColumnsDialog(page: Page): Promise<void> {
  await dismissTypeFilterPopup(page);
  await page.evaluate(() => {
    for (const d of DG.Dialog.getOpenDialogs())
      if (/Order or Hide Columns/i.test(d.title || '')) d.close();
  });
  await page.locator('[name="dialog-Order-or-Hide-Columns"]').waitFor({state: 'detached', timeout: 2000}).catch(() => {});
  await page.evaluate(() => {
    Array.from(document.querySelectorAll('[name="dialog-Order-or-Hide-Columns"]')).forEach((d) => d.remove());
  });
}

async function typeMenuLeafRect(page: Page, leafName: string): Promise<{x: number; y: number} | null> {
  return page.evaluate((name) => {
    const grp = Array.from(document.querySelectorAll('[name="div-Types"]'))
      .find((e: any) => e.offsetParent !== null) as HTMLElement | undefined;
    if (!grp) return null;
    const sub = grp.querySelector('.d4-menu-item-container.d4-vert-menu') as HTMLElement | null;
    if (sub) sub.style.display = 'flex';
    const els = Array.from(document.querySelectorAll(`[name="${name}"]`)) as HTMLElement[];
    for (const el of els) {
      if (el.offsetParent === null) continue;
      const r = el.getBoundingClientRect();
      if (r.width > 0 && r.height > 0) return {x: r.x + r.width / 2, y: r.y + r.height / 2};
    }
    return null;
  }, leafName);
}

function typeItemChecked(page: Page, typeName: string): Promise<boolean> {
  return page.evaluate((tn) => {
    const items = Array.from(document.querySelectorAll(`[name="div-Types---${tn}"]`))
      .filter((e: any) => e.offsetParent !== null) as HTMLElement[];
    return items.some((it) => it.querySelector('.d4-menu-item-check i')?.classList.contains('fa-check') ?? false);
  }, typeName);
}

function clickVisibleLeaf(page: Page, name: string): Promise<void> {
  return page.evaluate((name) => {
    const el = Array.from(document.querySelectorAll(`[name="${name}"]`))
      .find((e: any) => e.offsetParent !== null) as HTMLElement | undefined;
    if (!el) return;
    const r = el.getBoundingClientRect();
    const o = {bubbles: true, cancelable: true, clientX: r.x + r.width / 2, clientY: r.y + r.height / 2, button: 0} as any;
    for (const t of ['mouseover', 'mousemove', 'mousedown', 'mouseup', 'click'])
      el.dispatchEvent(new MouseEvent(t, o));
  }, name);
}

async function driveTypeFilter(page: Page, typeName: string): Promise<{intChecked: boolean}> {
  const leaf = await typeMenuLeafRect(page, `div-Types---${typeName}`);
  if (leaf) {
    await clickVisibleLeaf(page, `div-Types---${typeName}`);
    return {intChecked: await v.pollValue(() => typeItemChecked(page, typeName), (c) => c, 500, 50)};
  }
  return {intChecked: await typeItemChecked(page, typeName)};
}

function allTypesSquare(page: Page): Promise<boolean> {
  return page.evaluate(() => {
    const grp = Array.from(document.querySelectorAll('[name="div-Types"]'))
      .find((e: any) => e.offsetParent !== null) as HTMLElement | undefined;
    const sub = grp?.querySelector('.d4-menu-item-container.d4-vert-menu') as HTMLElement | null;
    if (sub) sub.style.display = 'flex';
    const items = Array.from(document.querySelectorAll('[name^="div-Types---"]'))
      .filter((e: any) => e.offsetParent !== null && (e.getAttribute('name') ?? '') !== 'div-Types---Reset-filter') as HTMLElement[];
    if (items.length === 0) return true;
    return items.every((it) => !(it.querySelector('.d4-menu-item-check i')?.classList.contains('fa-check') ?? false));
  });
}

async function driveTypeFilterReset(page: Page): Promise<{allSquare: boolean}> {
  const reset = await typeMenuLeafRect(page, 'div-Types---Reset-filter');
  if (reset) {
    await clickVisibleLeaf(page, 'div-Types---Reset-filter');
    return {allSquare: await v.pollValue(() => allTypesSquare(page), (s) => s, 500, 50)};
  }
  return {allSquare: await allTypesSquare(page)};
}
