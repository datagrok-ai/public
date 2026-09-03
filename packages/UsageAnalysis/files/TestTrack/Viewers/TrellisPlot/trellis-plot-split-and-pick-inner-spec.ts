/* ---
realizes: [trellisplot.cp.split-and-pick-inner, trellisplot.int.split-columns-drive-inner-viewer-grid, trellisplot.int.viewer-type-change-control-panel-axes]
--- */
import {expect, Page} from '@playwright/test';
import {localTest as test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep, isLocalBootNoise} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

// Scenarios 1-2 of the split-and-pick scenario on the local lane; Scenario 3 (layout and project
// round-trips) is trellis-plot-split-and-pick-inner-server-spec.ts.
test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

const isBenignError = (text: string) =>
  /Failed to load resource/.test(text) || /404 \(\)/.test(text) || /favicon/.test(text) ||
  /Unable to find element in cloned iframe/.test(text) ||
  /NullError: method not found: '\w+' on null/.test(text) || isLocalBootNoise(text);

// The column popup is a canvas grid with a 20px header row: row i is centred at grid.top + 20 + 16 * i + 8
// (measured 2026-09-03 on dev). The grid canvas moves down when the search box opens, so rows are
// located from the canvas rect, and the hover preview only follows a pointer that MOVES inside the
// popup, so the scan enters through the header rather than jumping onto a row.
const POPUP_HEADER_H = 20;
const POPUP_ROW_H = 16;

const readXState = (page: Page) => page.evaluate(() => {
  const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
  return {x: [...tp.props.xColumnNames] as string[],
    xCat: tp.xCategoriesCount as number, yCat: tp.yCategoriesCount as number,
    cells: document.querySelectorAll('[name="viewer-Trellis-plot"] .d4-trellis-plot-cell').length};
});

const popupRect = (page: Page) => page.evaluate(() => {
  const b = document.querySelector('.d4-column-selector-backdrop')!.getBoundingClientRect();
  return {left: b.left, top: b.top, width: b.width, height: b.height};
});

// the grid canvas: the sized canvas inside the popup that fits its width
const popupGridRect = (page: Page) => page.evaluate(() => {
  const popup = document.querySelector('.d4-column-selector-backdrop')!;
  const pb = popup.getBoundingClientRect();
  const fits = Array.from(popup.querySelectorAll('canvas'))
    .map((c) => c.getBoundingClientRect())
    .filter((b) => b.width > 0 && b.height > 0 && b.width <= pb.width && b.top >= pb.top)
    .sort((a, b) => b.width * b.height - a.width * a.height);
  const g = fits[0] ?? pb;
  return {left: g.left, top: g.top, width: g.width, height: g.height};
});

async function hoverPopupRow(page: Page, row: number): Promise<{hx: number; hy: number}> {
  const rect = await popupRect(page);
  const grid = await popupGridRect(page);
  const hx = rect.left + Math.min(60, rect.width / 2);
  const hy = grid.top + POPUP_HEADER_H + POPUP_ROW_H * row + POPUP_ROW_H / 2;
  await page.mouse.move(hx, rect.top + 4);
  await page.mouse.move(hx, hy);
  return {hx, hy};
}

async function closePopup(page: Page): Promise<void> {
  for (let i = 0; i < 3; i++) {
    if (!await page.evaluate(() => !!document.querySelector('.d4-column-selector-backdrop'))) return;
    await page.keyboard.press('Escape').catch(() => {});
    await v.pollValue(() => page.evaluate(() => !document.querySelector('.d4-column-selector-backdrop')), (ok) => ok, 800, 25);
  }
  await page.evaluate(() => document.body.dispatchEvent(new MouseEvent('mousedown', {bubbles: true}))).catch(() => {});
}

// Finds the popup row whose hover preview adds `column` to the X split, by the preview itself.
async function findRowByPreview(page: Page, column: string): Promise<{hx: number; hy: number} | null> {
  const grid = await popupGridRect(page);
  const rows = Math.floor((grid.height - POPUP_HEADER_H) / POPUP_ROW_H);
  for (let row = 0; row < rows; row++) {
    const pt = await hoverPopupRow(page, row);
    const state = await v.pollValue(() => readXState(page), (st) => st.x[st.x.length - 1] === column, 700, 25);
    if (state.x[state.x.length - 1] === column) return pt;
  }
  return null;
}

// Types `column` into the popup's search box (which appears on the first letter and takes the
// keyboard focus), so the list is filtered to that column and its row is row 0.
async function searchPopup(page: Page, column: string): Promise<void> {
  const searchFocused = () => page.evaluate(() =>
    document.activeElement?.classList.contains('d4-column-selector-search-input') ?? false);
  for (let i = 0; i < 2 && !await searchFocused(); i++) {
    await page.keyboard.press(column[0].toLowerCase());
    await v.pollValue(searchFocused, (ok) => ok, 1000, 25);
  }
  expect(await searchFocused(), 'the column popup search box did not open').toBe(true);
  await page.keyboard.press('Control+A');
  await page.keyboard.type(column.toLowerCase());
  const typed = await v.pollValue(() => page.evaluate(() =>
    (document.activeElement as HTMLInputElement | null)?.value ?? ''), (t) => t === column.toLowerCase(), 1000, 25);
  expect(typed).toBe(column.toLowerCase());
}

test('Trellis plot: split columns, inner-type switching', async ({page}) => {
  test.setTimeout(240_000);

  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  const onPageError = (e: Error) => { pageErrors.push(String(e)); };
  const onConsole = (m: any) => { if (m.type() === 'error' && !isBenignError(m.text())) consoleErrors.push(m.text()); };
  page.on('pageerror', onPageError);
  page.on('console', onConsole);

  await openDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  const addPageErrBefore = pageErrors.length;
  const addConsoleErrBefore = consoleErrors.length;
  await v.addViewerByIcon(page, 'trellis-plot', 'Trellis-plot', 15000);
  await v.waitForViewerRendered(page, 'Trellis plot', 900);

  await softStep('Scenario 1 Step 1', async () => {
    await expect(page.locator('[name="viewer-Trellis-plot"]')).toHaveCount(1);

    expect(pageErrors.slice(addPageErrBefore),
      'adding the Trellis Plot viewer raised an uncaught page error (github-964 smoke guard)').toEqual([]);
    expect(consoleErrors.slice(addConsoleErrBefore),
      'adding the Trellis Plot viewer raised a non-benign console error (github-964 smoke guard)').toEqual([]);
  });

  const cardinalities = await page.evaluate(() => {
    const df = grok.shell.tv.dataFrame;
    return {
      sex: df.col('SEX').categories.length,
      race: df.col('RACE').categories.length,
      control: df.col('CONTROL').categories.length,
    };
  });
  expect(cardinalities).toEqual({sex: 2, race: 4, control: 2});

  const cellLocator = page.locator('[name="viewer-Trellis-plot"] .d4-trellis-plot-cell');

  await softStep('Scenario 1 Step 4', async () => {
    const counts = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      await (window as any).__settled('viewer:Trellis plot.onViewerRendered', () => {
        tp.props.xColumnNames = ['SEX'];
        tp.props.yColumnNames = ['RACE'];
      }, 1500);
      return {xCat: tp.xCategoriesCount, yCat: tp.yCategoriesCount};
    });

    await expect(cellLocator).toHaveCount(counts.xCat * counts.yCat);
    await expect(cellLocator).toHaveCount(8);
    expect(counts).toEqual({xCat: 2, yCat: 4});
  });

  await softStep('Scenario 1 Step 6', async () => {
    const backdrop = page.locator('.d4-column-selector-backdrop');
    const openAddXPopup = async (): Promise<boolean> => {
      const plus = page.locator('[name="viewer-Trellis-plot"] [name="add-x-column"]').first();
      const box = await plus.boundingBox({timeout: 8000}).catch(() => null);
      if (!box) return false;
      await page.mouse.move(box.x + box.width / 2, box.y + box.height / 2);
      await page.mouse.down();
      await page.mouse.up();
      const opened = await backdrop.waitFor({timeout: 6000}).then(() => true).catch(() => false);
      if (opened) await page.evaluate(() => { (window as any).__popupFocus = document.activeElement; });
      return opened;
    };

    try {
      expect(await openAddXPopup(),
        'the (+) add-X-column control did not open the column selector popup').toBe(true);
      const first = await findRowByPreview(page, 'CONTROL');
      expect(first,
        'CONTROL row was not found by hover-preview scan of the popup — no blind-position fallback is taken').not.toBeNull();
      const previewed = await v.pollValue(() => readXState(page), (st) => st.cells !== 8 || st.x.length !== 1, 900, 50);
      expect(previewed.cells !== 8 || previewed.x.length !== 1,
        'hover preview never rebuilt the grid during the dwell — preview-regression signal (GROK-19673 class)').toBe(true);
      // the dwell: the preview lands on the popup's own timer, and an Escape that beats it leaves the
      // previewed column committed after the popup is gone (seen 2026-09-03), so hold until it is still
      await page.evaluate(() => (window as any).__settledFor(() => {
        const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
        return tp.props.xColumnNames.join(',') + '|' + document.querySelectorAll('[name="viewer-Trellis-plot"] .d4-trellis-plot-cell').length;
      }, 900, 1500, 25));
      // the preview rebuild can move keyboard focus off the popup, and an Escape that lands on the
      // viewer instead leaves the popup open with the preview committed (seen 2026-09-03, 1 run in 4)
      const focused = await v.pollValue(() => page.evaluate(() => {
        const a = document.activeElement as HTMLElement | null;
        return !!a && (a.classList.contains('d4-column-selector-search-input') || a === (window as any).__popupFocus);
      }), (ok) => ok, 500, 25);
      if (!focused) await page.evaluate(() => (window as any).__popupFocus?.focus?.());

      await page.keyboard.press('Escape');
      const closedOnFirstEsc = await backdrop
        .waitFor({state: 'detached', timeout: 2000}).then(() => true).catch(() => false);
      if (!closedOnFirstEsc) {
        await page.keyboard.press('Escape');
        await backdrop.waitFor({state: 'detached', timeout: 5000});
      }
      const afterEsc = await v.pollValue(() => readXState(page), (st) => st.cells === 8, 3000, 50);
      console.log(`[Step 6] popupFocused=${focused} closedOnFirstEsc=${closedOnFirstEsc} afterEsc=${JSON.stringify(afterEsc)}`);
      await expect(cellLocator).toHaveCount(8);
      expect((await readXState(page)).x).toEqual(['SEX']);

      expect(await openAddXPopup(),
        'the (+) add-X-column control did not reopen the column selector popup').toBe(true);
      await searchPopup(page, 'CONTROL');
      await page.keyboard.press('Enter');
      await backdrop.waitFor({state: 'detached', timeout: 6000});

      const committed = await v.pollValue(() => readXState(page), (s) => s.xCat === 4 && s.cells === 16, 1500, 50);
      expect(committed.x).toEqual(['SEX', 'CONTROL']);
      expect({xCat: committed.xCat, yCat: committed.yCat}).toEqual({xCat: 4, yCat: 4});
      await expect(cellLocator).toHaveCount(16);
    } finally {
      await closePopup(page);
    }
  });

  await softStep('Scenario 1 Step 7', async () => {
    const BORDER = 1;
    const HEADER_H = 20;
    const ROW_H = 16;
    const AUTOSIZE_PAD = 10;

    const backdrop = page.locator('.d4-column-selector-backdrop');

    const slot = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]');
      const plus = root ? root.querySelector('[name="add-x-column"]') : null;
      const host = plus ? plus.parentElement : null;
      if (!host) return null;
      const slots = Array.from(host.children).filter((el) =>
        el !== plus && el.classList.contains('d4-column-selector')) as HTMLElement[];
      const labels = slots.map((el) =>
        (el.querySelector('.d4-column-selector-column') as HTMLElement | null)?.innerText?.trim() ?? '');
      const hits = labels.filter((l) => l === 'CONTROL').length;
      if (hits !== 1) return {labels, hits, cx: -1, cy: -1};
      const b = slots[labels.indexOf('CONTROL')].getBoundingClientRect();
      return {labels, hits, cx: b.left + b.width / 2, cy: b.top + b.height / 2};
    });
    expect(slot,
      'the X selectors host ([name="add-x-column"] parent) was not found — Step 6 must have committed CONTROL first').not.toBeNull();
    expect(slot!.hits,
      `expected exactly one CONTROL slot among the X selectors, found labels ${JSON.stringify(slot?.labels)} — nothing is clicked`).toBe(1);

    await page.mouse.move(slot!.cx, slot!.cy);
    await page.mouse.down();
    await page.mouse.up();
    const opened = await backdrop.waitFor({timeout: 6000}).then(() => true).catch(() => false);
    expect(opened, 'the CONTROL X-slot did not open the column selector popup').toBe(true);

    try {
      const rect = await popupRect(page);

      const modelRows = (rect.height - (2 * BORDER + HEADER_H + AUTOSIZE_PAD)) / ROW_H;
      expect(Number.isInteger(modelRows) && modelRows >= 2,
        `column picker measured ${rect.width}x${rect.height}: its height is not 32 + 16*rows, so the layout model that locates the blank first row no longer holds — re-derive the offsets from these numbers`).toBe(true);

      const {hx, hy} = await hoverPopupRow(page, 0);
      const previewed = await v.pollValue(() => readXState(page), (s) => s.x.length === 1, 900, 25);

      expect(previewed.x,
        `hovering the blank first row at (${Math.round(hx)}, ${Math.round(hy)}) of a ${rect.width}x${rect.height} popup left the X split columns as ${JSON.stringify(previewed.x)} — the pointer did not land on the blank row, so nothing is committed`).toEqual(['SEX']);

      await page.mouse.click(hx, hy);
      await backdrop.waitFor({state: 'detached', timeout: 6000});

      const after = await v.pollValue(() => readXState(page), (s) => s.x.length === 1 && s.cells === 8, 1500, 50);
      expect(after.x).toEqual(['SEX']);
      expect({xCat: after.xCat, yCat: after.yCat}).toEqual({xCat: 2, yCat: 4});
      await expect(cellLocator).toHaveCount(8);
    } finally {
      await closePopup(page);
    }
  });

  await softStep('Scenario 2 Step 3', async () => {
    const perType = await page.evaluate(async () => {
      const w = window as any;
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;

      function cellSig(cellIdx: number): number | null {
        const cell = root.querySelectorAll('.d4-trellis-plot-cell')[cellIdx];
        const cv = cell?.querySelector('canvas') as HTMLCanvasElement | null;
        if (!cv) return null;
        try {
          const img = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
          let sig = 0;
          for (let i = 0; i < img.length; i += 4)
            sig = (sig * 31 + ((img[i] << 16) | (img[i + 1] << 8) | img[i + 2])) % 2147483647;
          return sig;
        } catch { return null; }
      }

      const typeIcons: Record<string, string> = {
        'Scatter plot': 'icon-scatter-plot', 'Bar chart': 'icon-bar-chart',
        'Histogram': 'icon-histogram', 'Pie chart': 'icon-pie-chart',
      };
      const out: Record<string, {eventType: string | null; distinct: boolean}> = {};
      for (const t of ['Scatter plot', 'Bar chart', 'Histogram', 'Pie chart']) {
        let eventType: string | null = null;
        const sub = tp.onEvent('d4-trellis-plot-viewer-type-changed').subscribe((arg: any) => {
          eventType = (typeof arg === 'string' ? arg : (arg?.args?.viewerType ?? null));
        });
        const vs = root.querySelector('[name="viewer selector"]') as HTMLElement;
        vs.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, button: 0}));
        await w.__poll(() => document.querySelector('.d4-combo-drop-down'), (e: Element | null) => !!e, 600, 40);
        const item = document.querySelector(`.d4-combo-drop-down [name="${typeIcons[t]}"]`);
        await w.__settled('viewer:Trellis plot.onViewerRendered',
          () => (item?.closest('.d4-list-item') as HTMLElement | null)?.click(), 1500);
        sub?.unsubscribe?.();

        let distinct = false;
        const deadline = Date.now() + 3000;
        do {
          const populated: number[] = [];
          const cells = root.querySelectorAll('.d4-trellis-plot-cell');
          for (let i = 0; i < cells.length && populated.length < 2; i++)
            if (cells[i].querySelector('canvas')) populated.push(i);
          if (populated.length === 2) {
            const a = cellSig(populated[0]);
            const b = cellSig(populated[1]);
            distinct = a !== null && b !== null && a !== b;
          }
          if (!distinct) await new Promise((r) => setTimeout(r, 150));
        } while (!distinct && Date.now() < deadline);
        out[t] = {eventType, distinct};
      }
      return out;
    });
    for (const t of ['Scatter plot', 'Bar chart', 'Histogram', 'Pie chart']) {
      expect(perType[t].eventType).toBe(t);
      expect(perType[t].distinct).toBe(true);
    }
  });

  await softStep('Scenario 2 Step 5', async () => {
    const result = await page.evaluate(async () => {
      const w = window as any;
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      const settle = (capMs: number, act: () => void) => w.__settled('viewer:Trellis plot.onViewerRendered', act, capMs);

      function populatedCellIdxs(limit: number): number[] {
        const idxs: number[] = [];
        const cells = root.querySelectorAll('.d4-trellis-plot-cell');
        for (let i = 0; i < cells.length && idxs.length < limit; i++)
          if (cells[i].querySelector('canvas')) idxs.push(i);
        return idxs;
      }
      function cellHash(cellIdx: number): number | null {
        const cell = root.querySelectorAll('.d4-trellis-plot-cell')[cellIdx];
        const cv = cell?.querySelector('canvas') as HTMLCanvasElement | null;
        if (!cv) return null;
        try {
          const img = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
          let h = 0;
          for (let i = 0; i < img.length; i += 4)
            h = (h * 31 + ((img[i] << 16) | (img[i + 1] << 8) | img[i + 2])) % 2147483647;
          return h;
        } catch { return null; }
      }

      await settle(1500, () => {
        tp.props.viewerType = 'Pie chart';
        tp.props.xColumnNames = ['SEX'];
        tp.props.yColumnNames = ['RACE'];
      });
      await settle(1200, () => {
        tp.props.yColumnNames = [];
        tp.props.viewerType = 'Bar chart';
      });
      await settle(1500, () => { tp.props.yColumnNames = ['RACE']; });
      const idxs = populatedCellIdxs(2);
      const barA = cellHash(idxs[0]);
      const barB = cellHash(idxs[1]);

      await settle(1800, () => { tp.props.viewerType = 'Pie chart'; });

      let pieA = cellHash(idxs[0]);
      let pieB = cellHash(idxs[1]);
      const deadline = Date.now() + 3000;
      while ((pieA === barA || pieB === barB || pieA === pieB) && Date.now() < deadline) {
        await new Promise((r) => setTimeout(r, 150));
        pieA = cellHash(idxs[0]);
        pieB = cellHash(idxs[1]);
      }

      return {
        hashesRead: barA !== null && barB !== null && pieA !== null && pieB !== null,
        distinctAfter: pieA !== null && pieB !== null && pieA !== pieB,
        reRendered: barA !== pieA && barB !== pieB,
      };
    });
    expect(result.hashesRead).toBe(true);
    expect(result.distinctAfter).toBe(true);
    expect(result.reRendered).toBe(true);
  });

  await softStep('Scenario 2 Step 7', async () => {
    const result = await page.evaluate(async () => {
      const w = window as any;
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      function selectorVisible(): boolean {
        const sel = root.querySelector('[name="viewer selector"]') as HTMLElement | null;
        if (!sel) return false;
        const r = sel.getBoundingClientRect();
        return r.width > 0 && r.height > 0;
      }
      tp.props.viewerType = 'Pie chart';
      tp.props.showControlPanel = true;
      const visibleBefore = await w.__poll(selectorVisible, (ok: boolean) => ok, 800, 25);

      tp.props.showControlPanel = false;
      const visibleAfterHide = await w.__poll(selectorVisible, (ok: boolean) => !ok, 1000, 25);
      const typeAfterHide = tp.props.viewerType;

      tp.props.showControlPanel = true;
      await w.__poll(selectorVisible, (ok: boolean) => ok, 800, 25);
      return {visibleBefore, visibleAfterHide, typeAfterHide};
    });
    expect(result.visibleBefore).toBe(true);
    expect(result.visibleAfterHide).toBe(false);
    expect(result.typeAfterHide).toBe('Pie chart');
  });

  page.off('pageerror', onPageError);
  page.off('console', onConsole);
  await v.closeAllAndWait(page);
  v.finishSpec();
});
