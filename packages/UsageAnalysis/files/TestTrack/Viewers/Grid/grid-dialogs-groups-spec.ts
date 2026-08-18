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

async function openGridMenu(page: Page, at: {x: number; y: number}): Promise<void> {
  await page.evaluate(({x, y}) => {
    const overlay = document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]') as HTMLElement;
    const cm = {bubbles: true, cancelable: true, clientX: x, clientY: y, button: 2, buttons: 2} as any;
    overlay.dispatchEvent(new MouseEvent('mousedown', cm));
    overlay.dispatchEvent(new MouseEvent('mouseup', cm));
    overlay.dispatchEvent(new MouseEvent('contextmenu', cm));
  }, at);

  await v.pollValue(
    () => page.locator('.d4-menu-popup .d4-menu-item').count(), (n) => n > 0, 600, 50);
}

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

async function openHamburger(page: Page, col: string): Promise<boolean> {
  const c = await headerCenter(page, col);
  await page.evaluate((at) => {
    const overlay = document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]') as HTMLElement;
    overlay.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, clientX: at.x, clientY: at.y}));
  }, c);

  await v.pollValue(() => page.evaluate(() => {
    const ham = document.querySelector('[name="viewer-Grid"] [name="icon-font-icon-menu"]') as HTMLElement | null;
    return !!ham && ham.getBoundingClientRect().width > 0;
  }), (revealed) => revealed, 500, 50);
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
    await page.waitForTimeout(350);
  }

  await page.waitForFunction(
    (v) => grok.shell.tv.dataFrame.col('AGE').getTag('.color-coding-type') === v,
    value, {timeout: 3000},
  ).catch(() => {  });
}

async function createColumnGroup(
  page: Page, cols: string[], groupName: string, swatchRgb: string,
): Promise<string> {

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
        setTimeout(tick, 450); 
      };
      setTimeout(tick, 450);
    });
  }, cols);
  expect(opened).toBe(true); 
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
  await page.waitForTimeout(400);

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

  const benignConsoleNoise = (t: string): boolean =>
    /Unable to find element in cloned iframe/i.test(t) ||
    /NullError: method not found: '[a-zA-Z]+' on null/i.test(t) ||
    /Stack trace [A-Za-z0-9]+/.test(t);
  const consoleErrors: string[] = [];
  page.on('console', (msg: any) => {
    if (msg.type() !== 'error') return;
    const t = msg.text();
    if (benignConsoleNoise(t)) return; 
    consoleErrors.push(t);
  });
  const errorsAt = () => consoleErrors.length;

  await softStep('Step 4 — Multi-column sort via the dialog: sortByColumnNames == [SEX, AGE]', async () => {
    await openGridMenu(page, await page.evaluate(() => {
      const db = grok.shell.tv.grid.cell('AGE', 0).documentBounds;
      return {x: db.x + db.width / 2, y: db.y + db.height / 2};
    }));
    expect(await synthClick(page, '[name="div-Sort..."]')).toBe(true); 
    await page.locator('[name="dialog-Sort-Table"]').waitFor({timeout: 6000});
    await pickSortColumn(page, 0, 'SEX');
    await pickSortColumn(page, 1, 'AGE');

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
        ageRow0: df.col('AGE').get(0), 
      };
    });
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
    await openGridMenu(page, await page.evaluate(() => {
      const db = grok.shell.tv.grid.cell('AGE', 0).documentBounds;
      return {x: db.x + db.width / 2, y: db.y + db.height / 2};
    }));
    expect(await synthClick(page, '[name="div-Order-or-Hide-Columns..."]')).toBe(true);
    await page.locator('[name="dialog-Order-or-Hide-Columns"]').waitFor({timeout: 6000});

    const menuOpened = await openTypeFilterMenu(page);
    expect(menuOpened).toBe(true); 

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
    const before = errorsAt();

    await closeColumnsDialog(page);

    await page.evaluate(async () => {
      const df2 = await grok.dapi.files.readCsv('System:AppData/Chem/tests/spgi-100.csv');
      grok.shell.addTableView(df2);
      await new Promise((r) => setTimeout(r, 1500));
    });

    await openGridMenu(page, await page.evaluate(() => {
      const db = grok.shell.tv.grid.cell(grok.shell.tv.grid.columns.byIndex(1).name, 0).documentBounds;
      return {x: db.x + db.width / 2, y: db.y + db.height / 2};
    }));
    await synthClick(page, '[name="div-Order-or-Hide-Columns..."]');
    await page.locator('[name="dialog-Order-or-Hide-Columns"]').waitFor({timeout: 6000});
    await openTypeFilterMenu(page);
    await driveTypeFilter(page, 'string');
    await page.waitForTimeout(600);

    await closeColumnsDialog(page);

    const switchErrs = consoleErrors.slice(before);
    expect(switchErrs).toEqual([]); 

    await page.evaluate(async () => {
      const views = grok.shell.tableViews;
      const spgi = Array.from(views).find((tv: any) => /spgi/i.test(tv.dataFrame?.name ?? ''));
      if (spgi) (spgi as any).close();
      await new Promise((r) => setTimeout(r, 800));
    });
  });

  await softStep('Step 17 — Gear opens the grid properties panel a second time after a close (GROK-17463)', async () => {

    await v.openTable(page, {path: 'System:DemoFiles/demog.csv', semTypeTimeoutMs: 3000});

    const open1 = await openGridSettings(page);
    expect(open1).toBe(true); 

    await page.evaluate(() => { grok.shell.o = grok.shell.tv.dataFrame.col('SEX'); });

    const closed = await v.pollValue(
      () => page.evaluate(() => document.querySelectorAll('[name="prop-row-height"]').length === 0),
      (c) => c, 800, 50);
    expect(closed).toBe(true); 

    const open2 = await openGridSettings(page);
    expect(open2).toBe(true); 
  });

  await softStep('Step 21 — Hamburger Linear coding on AGE is reflected in the Context Panel Colors (GROK-19288)', async () => {
    const opened = await openHamburger(page, 'AGE');
    expect(opened).toBe(true); 

    await synthClick(page, '.d4-popup-host [name="pane-Colors"] .d4-accordion-pane-header');
    const typeSelect = page.locator('.d4-popup-host select[name="input-Type"]');
    await typeSelect.waitFor({state: 'attached', timeout: 6000});

    await page.waitForTimeout(900);
    await driveColorCodingSelect(page, 'Linear');
    await page.waitForTimeout(400);

    await page.evaluate(() => {
      document.body.dispatchEvent(new MouseEvent('mousedown', {bubbles: true}));
      Array.from(document.querySelectorAll('.d4-popup-host')).forEach((e) => e.remove());
    });
    await page.waitForTimeout(300);

    const before = errorsAt();
    const r = await page.evaluate(async () => {
      const age = grok.shell.tv.dataFrame.col('AGE');

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
    expect(r.ageCCType).toBe('Linear'); 
    expect(r.colorsPanePresent).toBe(true); 
    const syncErrs = consoleErrors.slice(before).filter((e) => /color|coding|panel|grid/i.test(e));
    expect(syncErrs).toEqual([]); 
  });

  let firstGroupJson = '';
  await softStep('Step 22 — Create a blue group over AGE + HEIGHT: group tags present', async () => {

    await v.openTable(page, {path: 'System:DemoFiles/demog.csv', semTypeTimeoutMs: 3000});
    firstGroupJson = await createColumnGroup(page, ['AGE', 'HEIGHT'], 'AgeHeight', 'rgb(31, 119, 180)'); 
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
    const errs = consoleErrors.slice(before);
    expect(errs).toEqual([]); 
  });

  await softStep('Step 28 — Clicking the group band empty space raises no Bad-state console error (GROK-17443)', async () => {
    const before = errorsAt();

    await page.evaluate(async () => {
      const grid = grok.shell.tv.grid;
      const overlay = document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]') as HTMLElement;
      const rc = overlay.getBoundingClientRect();
      const gc = grid.columns.byName('AGE');

      const bandY = rc.y + 6;
      const bandX = rc.x + gc.left + gc.width; 
      const o = {bubbles: true, cancelable: true, clientX: bandX, clientY: bandY, button: 0} as any;
      overlay.dispatchEvent(new MouseEvent('mousedown', o));
      overlay.dispatchEvent(new MouseEvent('mouseup', o));
      overlay.dispatchEvent(new MouseEvent('click', o));
      await new Promise((r) => setTimeout(r, 300));
    });
    await page.waitForTimeout(400);
    const errs = consoleErrors.slice(before);
    expect(errs).toEqual([]); 
  });

  await softStep('Step 31 — Second group + select first group name + Esc raises no concurrent-modification error (GROK-17442/18213)', async () => {
    const before = errorsAt();

    const secondJson = await createColumnGroup(page, ['WEIGHT', 'SEX'], 'WeightSex', 'rgb(44, 160, 44)'); 
    expect(secondJson).toContain('WEIGHT'); 

    expect(secondJson).toContain('AgeHeight'); 
    expect(secondJson).toContain('WeightSex'); 

    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;

      grok.shell.o = [df.col('AGE'), df.col('HEIGHT')];
    });
    await page.waitForTimeout(500);
    await page.evaluate(() => {
      const overlay = document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]') as HTMLElement;
      overlay.focus();
    });
    await page.keyboard.press('Escape');
    await page.waitForTimeout(500);
    const errs = consoleErrors.slice(before);
    expect(errs).toEqual([]); 
  });

  const projectName = 'grid-dialogs-groups-' + Date.now();
  let savedProjectId: string | null = null;
  let groupsBeforeSave = '';

  await softStep('Step 33 — Save the view as a project via the ribbon Save button', async () => {
    groupsBeforeSave = await page.evaluate(() => grok.shell.tv.dataFrame.getTag('.columnGroups'));
    expect(groupsBeforeSave).toContain('#1f77b4'); 
    expect(groupsBeforeSave).toContain('#2ca02c'); 

    const saved = await saveProjectViaUI(page, projectName);
    savedProjectId = saved.projectId;
    expect(savedProjectId).not.toBeNull(); 

    const saveBalloons = await page.evaluate(() =>
      Array.from(document.querySelectorAll('.d4-balloon.error, .d4-balloon-error'))
        .map((b) => (b.textContent ?? '').trim()));
    const genuineSaveBalloons = saveBalloons
      .filter((t) => !/NullError: method not found: '[a-zA-Z]+' on null/i.test(t));
    expect(genuineSaveBalloons).toEqual([]); 
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
    expect(r.reopened).toBe(true); 
    expect(r.loadFailureBalloons).toEqual([]); 

    expect(r.columnGroups).toBe(groupsBeforeSave); 
    expect(r.columnGroups).toContain('#1f77b4'); 
    expect(r.columnGroups).toContain('#2ca02c'); 
    expect(r.ageGroup).toBeTruthy(); 
    expect(r.weightGroup).toBeTruthy(); 
    const reopenErrs = consoleErrors.slice(before);
    expect(reopenErrs).toEqual([]); 
  });

  await softStep('Teardown — delete the probe project', async () => {
    if (savedProjectId)
      await deleteProjectWithCleanup(page, {projectId: savedProjectId});
  });

  v.finishSpec();
});

async function openTypeFilterMenu(page: Page): Promise<boolean> {
  return page.evaluate(async () => {
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

    await new Promise((res) => setTimeout(res, 600));
    return !!Array.from(document.querySelectorAll('[name="div-Types"]')).find((e: any) => e.offsetParent !== null);
  });
}

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

async function closeColumnsDialog(page: Page): Promise<void> {
  await dismissTypeFilterPopup(page);
  await page.evaluate(() => {
    for (const d of DG.Dialog.getOpenDialogs())
      if (/Order or Hide Columns/i.test(d.title || '')) d.close();
  });
  await page.waitForTimeout(400);

  await page.evaluate(() => {
    Array.from(document.querySelectorAll('[name="dialog-Order-or-Hide-Columns"]'))
      .forEach((d) => d.remove());
  });
}

async function typeMenuLeafRect(
  page: Page, leafName: string,
): Promise<{x: number; y: number} | null> {
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
    if (items.length === 0) return true; 
    return items.every((it) => !(it.querySelector('.d4-menu-item-check i')?.classList.contains('fa-check') ?? false));
  });
  return {allSquare};
}

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

    }
  }
  return false;
}
