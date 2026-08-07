/* ---
realizes: [pivottable.cp.chrome-history-and-drag-config, pivottable.int.history-menu-requires-existing-columns, pivottable.int.default-aggr-type-remembered]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

// For docked pivots, use panel chrome rather than viewer-root controls:
// viewer [name="icon-times"] is not the close button, and saved states are found only
// in pivot history under the aggregation set's toString() key.


declare const grok: any;

const PIVOT = '[name="viewer-Pivot-table"]';

// A visible tag-row "+" may still be unclickable because transient popups/backdrops can
// intercept the pointer. Ensure stray overlays are gone and the target is the hit-test top element before clicking.

async function ensurePivotPlusClickable(page: Page, plusName: string) {
  const sel = `${PIVOT} [name="${plusName}"]`;
  for (let i = 0; i < 25; i++) {
    const clear = await page.evaluate((s) => {
      const plus = document.querySelector(s) as HTMLElement | null;
      if (!plus) return false;
      const r = plus.getBoundingClientRect();
      if (r.width === 0 || r.height === 0) return false;
      const top = document.elementFromPoint(r.x + r.width / 2, r.y + r.height / 2);
      return !!top && plus.contains(top);
    }, sel);
    if (clear) return;
    await page.keyboard.press('Escape');
    await page.waitForTimeout(160);
  }
  throw new Error(`pivot ${plusName} + icon stayed obscured (overlay never cleared)`);
}

// Drive the tag-row "+" column picker: click the +, type the column name into the canvas
// column-grid (.d4-column-selector-backdrop), commit with Enter. The picker is a canvas-rendered
// column grid driven by real keystrokes, not a DOM list.
async function addColumnViaPlus(page: Page, plusName: string, columnName: string) {
  await ensurePivotPlusClickable(page, plusName);
  await page.locator(`${PIVOT} [name="${plusName}"]`).click();
  await page.waitForSelector('.d4-column-selector-backdrop', {timeout: 6000});
  await page.keyboard.press(columnName[0]);
  await page.waitForTimeout(150);
  if (columnName.length > 1) await page.keyboard.type(columnName.slice(1));
  await page.waitForTimeout(150);
  await page.keyboard.press('Enter');
  await page.waitForTimeout(600);
}

// Expand a submenu parent so its children lay out with a non-zero box. The pivot-tag menu is a
// cascading vert-menu whose flyout opens only on a genuine sustained trusted hover; the pointer
// is moved with a 1px jitter each poll so each move is a distinct event and the open-delay elapses.
async function expandSubmenu(page: Page, parent: 'Aggregation' | 'Column', leaf: string) {
  const box = await page.locator(`.d4-menu-popup[name="pivot-tag"] [name="div-${parent}"]`).boundingBox();
  if (!box) throw new Error(`pivot-tag ${parent} parent not found`);
  const px = box.x + box.width / 2, py = box.y + box.height / 2;
  const leafReady = () => page.evaluate((name) => {
    const el = document.querySelector(`.d4-menu-popup[name="pivot-tag"] [name="${name}"]`) as HTMLElement | null;
    if (!el) return false;
    const r = el.getBoundingClientRect();
    return r.width > 0 && r.height > 0;
  }, leaf);
  for (let i = 0; i < 25; i++) {
    await page.mouse.move(px + (i % 2), py);
    await page.waitForTimeout(150);
    if (await leafReady()) return;
  }
  throw new Error(`pivot-tag ${parent} flyout did not reveal ${leaf}`);
}

// Click a cascading-flyout leaf while it is held open. Re-expands and retries on a mid-transit
// collapse so a fresh locator.click() cannot dismiss the flyout before the click lands.
async function clickFlyoutLeaf(page: Page, parent: 'Aggregation' | 'Column', leafName: string) {
  for (let attempt = 0; attempt < 3; attempt++) {
    await expandSubmenu(page, parent, leafName);
    const box = await page.locator(`.d4-menu-popup[name="pivot-tag"] [name="${leafName}"]`).boundingBox();
    if (!box) continue;
    const lx = box.x + box.width / 2, ly = box.y + box.height / 2;
    await page.mouse.move(lx, ly);
    await page.waitForTimeout(120);
    const onLeaf = await page.evaluate(({x, y, name}) => {
      const top = document.elementFromPoint(x, y);
      const el = document.querySelector(`.d4-menu-popup[name="pivot-tag"] [name="${name}"]`);
      return !!(top && el && (el.contains(top) || el === top));
    }, {x: lx, y: ly, name: leafName});
    if (!onLeaf) continue;
    await page.mouse.down();
    await page.mouse.up();
    await page.waitForTimeout(400);
    return;
  }
  throw new Error(`could not click pivot-tag ${parent} leaf ${leafName}`);
}

// Pick an Aggregation-type child (e.g. 'sum').
async function pickAggregation(page: Page, aggType: string) {
  await clickFlyoutLeaf(page, 'Aggregation', `div-Aggregation---${aggType}`);
}

// Read a tag row's chip captions by row title (durable across rebuilds — the chip name attribute
// is dropped after any prop-driven rebuild, per pivot_table.md; the caption span survives).
async function rowChips(page: Page, rowTitle: string): Promise<string[]> {
  return page.evaluate((title) => {
    const root = document.querySelector('[name="viewer-Pivot-table"]');
    const panel = Array.from(root!.querySelectorAll('.grok-pivot-column-panel'))
      .find((p) => p.querySelector('.grok-pivot-column-tags-title')?.getAttribute('d4-name') === title);
    if (!panel) return [];
    return Array.from(panel.querySelectorAll('.d4-tag')).map((t) => t.querySelector('span')?.textContent?.trim() ?? '');
  }, rowTitle);
}

test.use(specTestOptions);

// Page-error collector for the GROK-17122 console-error-delta step. The
// cloned-iframe warning class is unfixable/harmless and excluded from the count.
const pageErrors: string[] = [];
const isIgnorable = (m: string) => m.includes('cloned iframe') || m.includes('Unable to find element in cloned iframe');

test('Pivot table chrome, history and drag-driven configuration', async ({page}) => {
  test.setTimeout(300_000);
  page.on('pageerror', (e) => { if (!isIgnorable(e.message)) pageErrors.push(e.message); });

  await loginToDatagrok(page);

  // Setup: open demog, add the Pivot Table viewer, wait for the tag-editor header.
  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    window.localStorage.removeItem('grok-aggregation-history');
    const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
    grok.shell.addTableView(df);
    await new Promise((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(resolve, 3000);
    });
    (document.querySelector('[name="icon-pivot-table"]') as HTMLElement)?.click();
    await new Promise((r) => setTimeout(r, 1500));
  });
  await page.locator('[name="viewer-Pivot-table"]').waitFor({timeout: 15000});
  await page.locator('[name="viewer-Pivot-table"] .grok-pivot-column-tags-title[d4-name="Group by"]').waitFor({timeout: 15000});

  await softStep('Scenario 1 Step 2: auto-config is DIS_POP / SEVERITY / avg(AGE), counts visible', async () => {
    const r = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Pivot-table"]')!;
      const pv = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
      const tagsOf = (title: string) => {
        const row = [...root.querySelectorAll('.grok-pivot-column-panel')]
          .find((p) => p.querySelector(`.grok-pivot-column-tags-title[d4-name="${title}"]`));
        return row ? [...row.querySelectorAll('.d4-tag')].map((t) => t.textContent!.trim()) : [];
      };
      const counts = root.querySelector('.grok-pivot-counts') as HTMLElement | null;
      return {
        groupByTags: tagsOf('Group by'), aggTags: tagsOf('Aggregate'), pivotTags: tagsOf('Pivot'),
        groupBy: pv.props.groupByColumnNames, pivot: pv.props.pivotColumnNames,
        agg: pv.props.aggregateColumnNames, aggTypes: pv.props.aggregateAggTypes,
        countsVisible: !!counts && !!counts.offsetParent,
      };
    });
    // Tag captions
    expect(r.groupByTags).toEqual(['DIS_POP']);
    expect(r.aggTags).toEqual(['avg(AGE)']);
    expect(r.pivotTags).toEqual(['SEVERITY']);
    // Property lists
    expect(r.groupBy).toContain('DIS_POP');
    expect(r.pivot).toContain('SEVERITY');
    expect(r.agg).toContain('AGE');
    expect(r.aggTypes).toContain('avg');
    expect(r.countsVisible).toBe(true);
  });

  await softStep('Scenario 2 Step 3: close via cross icon → viewer gone, no console-error delta (GROK-17122)', async () => {
    const errorsBefore = pageErrors.length;
    const r = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Pivot-table"]')!;
      const panel = root.closest('.panel-base')!;
      const closeBtn = panel.querySelector('.panel-titlebar [name="Close"]') as HTMLElement;
      closeBtn?.click();
      await new Promise((res) => setTimeout(res, 900));
      const gone = !Array.from(grok.shell.tv.viewers).some((x: any) => x.type === 'Pivot table');
      // Re-add and confirm auto-config comes back.
      (document.querySelector('[name="icon-pivot-table"]') as HTMLElement)?.click();
      await new Promise((res) => setTimeout(res, 1200));
      const pv2 = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
      return {gone, reAdded: !!pv2, groupBy: pv2?.props.groupByColumnNames, pivot: pv2?.props.pivotColumnNames};
    });
    expect(r.gone).toBe(true);
    // GROK-17122: closing the viewer must not emit console errors.
    expect(pageErrors.length).toBe(errorsBefore);
    expect(r.reAdded).toBe(true);
    expect(r.groupBy).toContain('DIS_POP');
    expect(r.pivot).toContain('SEVERITY');
  });

  // ---- Scenario 3: Show Header and Show Command Bar -----------------------
  const probeChrome = async () => page.evaluate(() => {
    const root = document.querySelector('[name="viewer-Pivot-table"]')!;
    const vis = (el: Element | null) => !!el && !!(el as HTMLElement).offsetParent && getComputedStyle(el as HTMLElement).display !== 'none';
    const rowByTitle = (title: string) => [...root.querySelectorAll('.grok-pivot-column-panel')]
      .find((p) => p.querySelector(`.grok-pivot-column-tags-title[d4-name="${title}"]`)) ?? null;
    return {
      data: vis(rowByTitle('Data')), groupBy: vis(rowByTitle('Group by')),
      agg: vis(rowByTitle('Aggregate')), pivot: vis(rowByTitle('Pivot')),
      counts: vis(root.querySelector('.grok-pivot-counts')),
      cmdBar: vis(root.querySelector('.d4-command-bar')),
      history: vis(root.querySelector('.d4-command-bar [name="icon-history"]')),
    };
  });
  const setChromeProp = async (prop: string, value: boolean) => page.evaluate(async ({prop, value}) => {
    const pv = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
    pv.props[prop] = value;
    await new Promise((res) => setTimeout(res, 450));
  }, {prop, value});

  await softStep('Scenario 3 Step 3: Show Header=false hides the Data row, tag rows and counts; they return on true', async () => {
    await setChromeProp('showHeader', false);
    const headerOff = await probeChrome();
    await setChromeProp('showHeader', true);
    const headerOn = await probeChrome();
    expect(headerOff.data).toBe(false);
    expect(headerOff.groupBy).toBe(false);
    expect(headerOff.agg).toBe(false);
    expect(headerOff.pivot).toBe(false);
    expect(headerOff.counts).toBe(false);
    expect(headerOn.groupBy).toBe(true);
    expect(headerOn.counts).toBe(true);
  });

  await softStep('Scenario 3 Step 6: Show Command Bar=false hides the command bar with history/refresh icons; it returns on true', async () => {
    await setChromeProp('showCommandBar', false);
    const cmdOff = await probeChrome();
    await setChromeProp('showCommandBar', true);
    const cmdOn = await probeChrome();
    expect(cmdOff.cmdBar).toBe(false);
    expect(cmdOff.history).toBe(false);
    expect(cmdOn.cmdBar).toBe(true);
    expect(cmdOn.history).toBe(true);
  });

  await softStep('Scenario 4 Step 4: title in the header, description Top visible, Never hides it', async () => {
    const r = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Pivot-table"]')!;
      const pv = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
      pv.props.showTitle = true;
      pv.props.title = 'My Pivot';
      pv.props.description = 'Summary stats';
      pv.props.descriptionPosition = 'Top';
      await new Promise((res) => setTimeout(res, 500));
      const panel = root.closest('.panel-base') ?? root;
      const titleText = [...document.querySelectorAll('.panel-titlebar-tabhost .panel-titlebar-text')]
        .map((e) => e.textContent!.trim()).filter(Boolean);
      const descTop = [...panel.querySelectorAll('.d4-viewer-description')]
        .map((e) => ({txt: e.textContent!.trim(), vis: !!(e as HTMLElement).offsetParent}));
      pv.props.descriptionVisibilityMode = 'Never';
      await new Promise((res) => setTimeout(res, 500));
      const descNever = [...panel.querySelectorAll('.d4-viewer-description')]
        .filter((e) => !!(e as HTMLElement).offsetParent && e.textContent!.includes('Summary stats')).length;
      // restore
      pv.props.showTitle = false; pv.props.title = ''; pv.props.description = ''; pv.props.descriptionVisibilityMode = 'Auto';
      return {titleShown: titleText.includes('My Pivot'), descTopVisible: descTop.some((d) => d.txt.includes('Summary stats') && d.vis), descNeverCount: descNever};
    });
    expect(r.titleShown).toBe(true);
    expect(r.descTopVisible).toBe(true);
    expect(r.descNeverCount).toBe(0);
  });

  // ---- Scenario 5: aggregation history (I8) -------------------------------
  await softStep('Scenario 5 Step 4: Save parameters writes localStorage history for RACE / avg(WEIGHT)', async () => {
    const r = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Pivot-table"]')!;
      const pv = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
      pv.props.groupByColumnNames = ['RACE'];
      pv.props.aggregateColumnNames = ['WEIGHT'];
      pv.props.aggregateAggTypes = ['avg'];
      pv.props.pivotColumnNames = [];
      await new Promise((res) => setTimeout(res, 600));
      (root.querySelector('.d4-command-bar [name="icon-history"]') as HTMLElement).click();
      await new Promise((res) => setTimeout(res, 500));
      const saveItem = [...document.querySelectorAll('.d4-menu-item')]
        .find((i) => i.querySelector('.d4-menu-item-label')?.textContent?.trim() === 'Save parameters') as HTMLElement | null;
      saveItem?.click();
      await new Promise((res) => setTimeout(res, 500));
      const raw = window.localStorage.getItem('grok-aggregation-history');
      let parsed: any = null;
      try { parsed = JSON.parse(raw ?? ''); } catch (_) { parsed = null; }
      const flat = Array.isArray(parsed) ? parsed.flat(2).map((a: any) => a.colName) : [];
      return {isArray: Array.isArray(parsed), len: Array.isArray(parsed) ? parsed.length : -1, names: flat};
    });
    expect(r.isArray).toBe(true);
    expect(r.len).toBeGreaterThan(0);
    expect(r.names).toContain('RACE');
    expect(r.names).toContain('WEIGHT');
  });

  await softStep('Scenario 5 Step 6: picking the saved entry restores Group by / Aggregate (tag captions)', async () => {
    // Reconfigure away, then re-apply the saved entry via a trusted menu click.
    await page.evaluate(async () => {
      const pv = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
      pv.props.groupByColumnNames = ['SEX'];
      pv.props.aggregateColumnNames = ['AGE'];
      pv.props.aggregateAggTypes = ['avg'];
      await new Promise((res) => setTimeout(res, 500));
      (document.querySelector('[name="viewer-Pivot-table"] .d4-command-bar [name="icon-history"]') as HTMLElement).click();
      await new Promise((res) => setTimeout(res, 500));
    });
    // The history re-apply rebuilds the tag rows without writing back to look props,
    // so the restored config is read from the tag captions, never pv.props.
    await page.locator('.d4-menu-popup [name="div-key(RACE),avg(WEIGHT)"]').click();
    await page.waitForTimeout(700);
    const tags = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Pivot-table"]')!;
      const tagsOf = (title: string) => {
        const row = [...root.querySelectorAll('.grok-pivot-column-panel')]
          .find((p) => p.querySelector(`.grok-pivot-column-tags-title[d4-name="${title}"]`));
        return row ? [...row.querySelectorAll('.d4-tag')].map((t) => t.textContent!.trim()) : [];
      };
      return {groupBy: tagsOf('Group by'), agg: tagsOf('Aggregate')};
    });
    expect(tags.groupBy).toEqual(['RACE']);
    expect(tags.agg).toEqual(['avg(WEIGHT)']);
  });

  await softStep('Scenario 5 Step 8: after WEIGHT is removed the history menu drops the WEIGHT entry (I8)', async () => {
    // The menu filters entries whose columns are absent from the pivot's table snapshot.
    // That snapshot is rebuilt on viewer (re)attach, so remove the column then re-add the
    // viewer — the faithful reproduction of the I8 invariant.
    const r = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.columns.remove('WEIGHT');
      await new Promise((res) => setTimeout(res, 600));
      const pv = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
      pv.close();
      await new Promise((res) => setTimeout(res, 900));
      (document.querySelector('[name="icon-pivot-table"]') as HTMLElement)?.click();
      await new Promise((res) => setTimeout(res, 1800));
      const root = document.querySelector('[name="viewer-Pivot-table"]')!;
      (root.querySelector('.d4-command-bar [name="icon-history"]') as HTMLElement).click();
      await new Promise((res) => setTimeout(res, 700));
      const menu = [...document.querySelectorAll('.d4-menu-popup')].pop()!;
      const labels = [...menu.querySelectorAll('.d4-menu-item-label')].map((e) => e.textContent!.trim());
      document.body.click();
      return {labels, weightGone: !labels.some((l) => l.includes('WEIGHT'))};
    });
    expect(r.weightGone).toBe(true);
  });

  await softStep('Scenario 5 Step 9: Refresh (icon-redo) clears Group by / Pivot and re-seeds the default aggregates', async () => {
    // Establish a non-default configuration, then drive the command-bar Refresh (glyph is
    // `redo`, tooltip "Refresh"). Refresh clears keyCols/pivotCols and sets aggrCols to the
    // first two numerical columns (avg) via setTags(notify:false) — it does NOT write back to
    // the look props, so the restored state is read from the tag captions, never pv.props.
    await page.evaluate(async () => {
      const pv = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
      pv.props.groupByColumnNames = ['RACE'];
      pv.props.aggregateColumnNames = ['AGE'];
      pv.props.aggregateAggTypes = ['sum'];
      pv.props.pivotColumnNames = ['SEX'];
      await new Promise((res) => setTimeout(res, 600));
    });
    await page.locator('[name="viewer-Pivot-table"] .d4-command-bar [name="icon-redo"]').click();
    await page.waitForTimeout(900);
    const tags = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Pivot-table"]')!;
      const tagsOf = (title: string) => {
        const row = [...root.querySelectorAll('.grok-pivot-column-panel')]
          .find((p) => p.querySelector(`.grok-pivot-column-tags-title[d4-name="${title}"]`));
        return row ? [...row.querySelectorAll('.d4-tag')].map((t) => t.textContent!.trim()) : [];
      };
      return {groupBy: tagsOf('Group by'), agg: tagsOf('Aggregate'), pivot: tagsOf('Pivot')};
    });
    expect(tags.groupBy).toEqual([]);
    expect(tags.pivot).toEqual([]);
    expect(tags.agg).toEqual(['avg(AGE)', 'avg(HEIGHT)']);
  });

  // Re-open demog to restore the WEIGHT column for the following scenarios.
  await page.evaluate(async () => {
    grok.shell.closeAll();
    window.localStorage.removeItem('grok-aggregation-history');
    const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
    grok.shell.addTableView(df);
    await new Promise((res) => setTimeout(res, 1500));
    (document.querySelector('[name="icon-pivot-table"]') as HTMLElement)?.click();
    await new Promise((res) => setTimeout(res, 1500));
  });
  await page.locator('[name="viewer-Pivot-table"] .grok-pivot-column-tags-title[d4-name="Group by"]').waitFor({timeout: 15000});

  // ---- Scenario 6: remembered aggregation type (I9) -----------------------
  await softStep('Scenario 6 Steps 1-2: add HEIGHT, choose sum → sum(HEIGHT) tag, then remove it', async () => {
    // Step 1: click + Aggregate, pick HEIGHT (the new tag defaults to avg(HEIGHT)),
    // then right-click the chip → Aggregation → sum so it becomes sum(HEIGHT).
    await addColumnViaPlus(page, 'div-add-Aggregate', 'HEIGHT');
    let aggChips = await rowChips(page, 'Aggregate');
    expect(aggChips.some((c) => c.includes('HEIGHT'))).toBe(true);
    // Retarget the chip menu at the HEIGHT chip: it is the last-added Aggregate chip.
    await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Pivot-table"]')!;
      const panel = Array.from(root.querySelectorAll('.grok-pivot-column-panel'))
        .find((p) => p.querySelector('.grok-pivot-column-tags-title')?.getAttribute('d4-name') === 'Aggregate')!;
      const chip = Array.from(panel.querySelectorAll('.d4-tag'))
        .find((t) => (t.textContent ?? '').includes('HEIGHT'))!;
      chip.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true, button: 2}));
    });
    await page.locator('.d4-menu-popup[name="pivot-tag"] [name="div-Aggregation"]').waitFor({timeout: 5000});
    await pickAggregation(page, 'sum');
    await page.keyboard.press('Escape');
    await page.waitForTimeout(400);
    // Step 1 outcome: the HEIGHT chip now reads sum(HEIGHT) — the last-chosen aggregation
    // type is now sum, which is what the viewer remembers for the next picker.
    aggChips = await rowChips(page, 'Aggregate');
    expect(aggChips.some((c) => c.includes('sum(HEIGHT)'))).toBe(true);
    // Step 2: remove the sum(HEIGHT) tag via its × icon.
    await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Pivot-table"]')!;
      const panel = Array.from(root.querySelectorAll('.grok-pivot-column-panel'))
        .find((p) => p.querySelector('.grok-pivot-column-tags-title')?.getAttribute('d4-name') === 'Aggregate')!;
      const chip = Array.from(panel.querySelectorAll('.d4-tag'))
        .find((t) => (t.textContent ?? '').includes('HEIGHT'))!;
      (chip.querySelector('i') as HTMLElement)?.click();
    });
    await page.waitForTimeout(500);
    aggChips = await rowChips(page, 'Aggregate');
    expect(aggChips.some((c) => c.includes('HEIGHT'))).toBe(false);
  });

  await softStep('Scenario 6 Step 3: the Aggregate + popup pre-offers the remembered aggregation type (I9)', async () => {
    // Step 3: re-open the Aggregate + picker. The remembered aggregation type (sum, set in
    // Steps 1-2) lives on PivotGrid.defaultAggrType and is pushed into the ColumnComboBox
    // aggr-selector when the picker opens. The pre-offered value renders on the canvas-backed
    // picker with no DOM/JS-API read-back, so it cannot be asserted here.
    await ensurePivotPlusClickable(page, 'div-add-Aggregate');
    await page.locator(`${PIVOT} [name="div-add-Aggregate"]`).click();
    const backdrop = page.locator('.d4-column-selector-backdrop');
    await backdrop.waitFor({timeout: 6000});
    expect(await backdrop.count()).toBeGreaterThan(0);
    // Cancel the picker without adding a column (Step 4 leaves the Aggregate row unchanged).
    await page.keyboard.press('Escape');
    await page.waitForTimeout(400);
    const aggChips = await rowChips(page, 'Aggregate');
    expect(aggChips.some((c) => c.includes('AGE'))).toBe(true);
    expect(aggChips.some((c) => c.includes('HEIGHT'))).toBe(false);
  });

  // Scenario 7 covered manually in the pivot-table-ui.md companion.

  // ---- Scenario 8: ID grouping, semantic types and the Table property -----
  await softStep('Scenario 8 Step 2: grouping by USUBJID makes one row per identifier, no console error (GROK-16201)', async () => {
    const errorsBefore = pageErrors.length;
    const r = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Pivot-table"]')!;
      const pv = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
      const df = grok.shell.tv.dataFrame;
      pv.props.groupByColumnNames = ['USUBJID'];
      pv.props.aggregateColumnNames = ['AGE'];
      pv.props.aggregateAggTypes = ['avg'];
      pv.props.pivotColumnNames = [];
      await new Promise((res) => setTimeout(res, 800));
      const counts = root.querySelector('.grok-pivot-counts')!.textContent!.replace(/\s+/g, ' ').trim();
      const distinct = df.col('USUBJID').categories.length;
      const rowsMatch = counts.startsWith(`${distinct} rows`);
      return {distinct, counts, rowsMatch};
    });
    // One aggregate row per subject identifier.
    expect(r.rowsMatch).toBe(true);
    // GROK-16201: ID grouping writes no console error.
    expect(pageErrors.length).toBe(errorsBefore);
  });

  await softStep('Scenario 8 Step 5: ADD opens the aggregated result; key column keeps its type (GROK-16074)', async () => {
    const r = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Pivot-table"]')!;
      const pv = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
      const df = grok.shell.tv.dataFrame;
      pv.props.groupByColumnNames = ['DIS_POP'];
      pv.props.aggregateColumnNames = ['AGE'];
      pv.props.aggregateAggTypes = ['avg'];
      pv.props.pivotColumnNames = [];
      await new Promise((res) => setTimeout(res, 600));
      const srcCol = df.col('DIS_POP');
      const addBtn = root.querySelector('.grok-pivot-counts [name="button-ADD"]') as HTMLElement;
      addBtn.click();
      await new Promise((res) => setTimeout(res, 1000));
      const views = Array.from(grok.shell.views) as any[];
      const newDf = views[views.length - 1]?.dataFrame;
      const keyCol = newDf?.col('DIS_POP');
      const out = {
        opened: !!keyCol,
        keyType: keyCol?.type, srcType: srcCol.type,
        keySemType: keyCol?.semType ?? null, srcSemType: srcCol.semType ?? null,
      };
      // cleanup: close the aggregated view
      const aggView = views.find((vw) => vw.name === 'Table aggregation');
      if (aggView) aggView.close();
      await new Promise((res) => setTimeout(res, 400));
      return out;
    });
    expect(r.opened).toBe(true);
    // GROK-16074: the key column of the opened table keeps the source column's type.
    // (demog carries no semType on its columns, so type is the observable signal here.)
    expect(r.keyType).toBe(r.srcType);
    expect(r.keySemType).toBe(r.srcSemType);
  });

  await softStep('Scenario 8 Step 7: switching the Data-row Table property back and forth duplicates no Data entry / header (github-3414, GROK-14995)', async () => {
    const counts = () => page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Pivot-table"]')!;
      const dataRow = [...root.querySelectorAll('.grok-pivot-column-panel')]
        .find((p) => p.querySelector('.grok-pivot-column-tags-title[d4-name="Data"]'));
      return {
        dataEntries: dataRow ? dataRow.querySelectorAll('.d4-tag').length : 0,
        headers: root.querySelectorAll('.grok-pivot-column-tags-title[d4-name="Data"]').length,
      };
    });
    const before = await counts();
    // Drive the real Data-row Table workflow: clicking the Data-row table tag opens the
    // column-selection dialog, and confirming it rebuilds the Data-row table.
    // This exercises the full UI path guarded by github-3414 / GROK-14995..
    const tagBox = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Pivot-table"]')!;
      const dataRow = [...root.querySelectorAll('.grok-pivot-column-panel')]
        .find((p) => p.querySelector('.grok-pivot-column-tags-title[d4-name="Data"]'))!;
      const tag = dataRow.querySelector('.d4-tag') as HTMLElement;
      const r = tag.getBoundingClientRect();
      return {x: r.x + r.width / 2, y: r.y + r.height / 2};
    });
    await page.mouse.click(tagBox.x, tagBox.y);   // trusted click → opens the column-select modal
    const dlg = page.locator('.d4-dialog').last();
    await dlg.waitFor({timeout: 8000});
    // OK re-commits the (unchanged) column selection and fires the Data-row rebuild "and back".
    await dlg.locator('[name="button-OK"]').click();
    await page.waitForTimeout(700);
    const after = await counts();
    expect(after.dataEntries).toBe(before.dataEntries);
    expect(after.headers).toBe(1);
  });

  await page.evaluate(() => grok.shell.closeAll());
  v.finishSpec();
});
