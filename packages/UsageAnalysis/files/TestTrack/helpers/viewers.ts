import {Page, expect} from '@playwright/test';
import {createHash} from 'crypto';
import {stepErrors, StepError} from '../spec-login';

export interface OpenTableOptions {

  path?: string;

  withFilterPanel?: boolean;

  semType?: 'Molecule' | 'Macromolecule';

  sdf?: boolean;

  settleMs?: number;

  semTypeTimeoutMs?: number;
}

export async function openTable(page: Page, options?: OpenTableOptions): Promise<void> {
  const p = options?.path ?? 'System:AppData/Chem/tests/spgi-100.csv';
  const useOpenFile = options?.sdf === true || /\.(sdf|nwk|pdb)$/i.test(p);
  const semTypeTimeoutMs = options?.semTypeTimeoutMs ?? 5000;
  await page.evaluate(async ({path, openFile, semTypeTimeoutMs}) => {
    document.body.classList.add('selenium');
    (window as any).grok.shell.settings.showFiltersIconsConstantly = true;
    (window as any).grok.shell.windows.simpleMode = true;
    (window as any).grok.shell.closeAll();
    let df: any;
    if (openFile) {
      const DG = (window as any).DG;
      const grok = (window as any).grok;
      const fns = DG.Func.find({name: 'OpenFile'});
      if (!fns?.length) throw new Error('OpenFile function not registered');
      const fullPathForOpen = path.replace(/^([^:/]+):/, '$1.');
      await fns[0].prepare({fullPath: fullPathForOpen}).call(undefined, undefined, {processed: false});
      for (let i = 0; i < 24; i++) {
        const tv = grok.shell.tv;
        if (tv?.dataFrame && typeof tv.addViewer === 'function') {
          df = tv.dataFrame;
          break;
        }
        await new Promise((r) => setTimeout(r, 500));
      }
      if (!df) throw new Error(`OpenFile("${path}") did not produce a TableView (12s settle)`);
    } else {
      df = await (window as any).grok.dapi.files.readCsv(path);
      (window as any).grok.shell.addTableView(df);
    }
    await new Promise((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(null); });
      setTimeout(resolve, semTypeTimeoutMs);
    });
    const hasBioChem = Array.from({length: df.columns.length}, (_, i: number) => df.columns.byIndex(i))
      .some((c: any) => c.semType === 'Molecule' || c.semType === 'Macromolecule');
    if (hasBioChem) {
      for (let i = 0; i < 50; i++) {
        const grid = document.querySelector('[name="viewer-Grid"]');
        if (grid?.querySelector('canvas')) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      await new Promise((r) => setTimeout(r, 5000));
    }
  }, {path: p, openFile: useOpenFile, semTypeTimeoutMs});
  await page.locator('.d4-grid[name="viewer-Grid"]').first().waitFor({timeout: 30000});
  if (options?.settleMs) await page.waitForTimeout(options.settleMs);
  if (options?.withFilterPanel) await openFilterPanel(page);
}

export async function openFilterPanel(page: Page): Promise<void> {
  await page.evaluate(() => (window as any).grok.shell.tv.getFiltersGroup());
  // a molecule column builds a substructure filter here, which has been seen to take
  // longer than 15s on dev — the wait exits on the element, so the higher cap is free
  await page.locator('[name="viewer-Filters"] .d4-filter').first().waitFor({timeout: 30000});

  await page.locator('[name="viewer-Filters"] .d4-filter').first().hover();
}

export async function addViewerByIcon(
  page: Page, iconName: string, viewerName: string, timeoutMs = 5000, expectedType?: string,
): Promise<void> {
  await page.evaluate((n) => {
    (document.querySelector('[name="icon-' + n + '"]') as HTMLElement).click();
  }, iconName);
  await page.locator('[name="viewer-' + viewerName + '"]').waitFor({timeout: timeoutMs});

  await page.waitForFunction((vn) => {
    const tv = (window as any).grok?.shell?.tv;
    if (!tv) return false;
    const norm = (s: string) => s.replace(/[\s-]+/g, ' ').toLowerCase();
    const v = Array.from(tv.viewers).find((x: any) => norm(x.type) === norm(vn));
    return !!(v && (v as any).props);
  }, expectedType ?? viewerName, {timeout: timeoutMs});
}

export const LEGEND_COLUMN_PROP: Record<string, {prop: string; array: boolean}> = {
  'Scatter plot':  {prop: 'colorColumnName', array: false},
  'Histogram':     {prop: 'splitColumnName', array: false},
  'Line chart':    {prop: 'splitColumnName', array: false},
  'Bar chart':     {prop: 'splitColumnName', array: false},
  'Pie chart':     {prop: 'categoryColumnName', array: false},
  'Trellis plot':  {prop: 'xColumnNames', array: true},
  'Box plot':      {prop: 'categoryColumnNames', array: true},
};

export interface ViewerSpec {

  type: string;

  column?: string;

  prop?: string;

  array?: boolean;
}

export async function addLegendViewers(
  page: Page,
  options: {column: string; viewers: (string | ViewerSpec)[]; settleMs?: number},
): Promise<void> {
  const specs: ViewerSpec[] = options.viewers.map((v) => typeof v === 'string' ? {type: v} : v);
  const settleMs = options.settleMs ?? 1500;
  await page.evaluate(async ({s, defaultCol, settle, map}) => {
    const tv = (window as any).grok.shell.tv;
    for (const spec of s) {
      tv.addViewer(spec.type);
      await new Promise((r) => setTimeout(r, 300));
    }
    for (const v of tv.viewers) {
      if (v.type === 'Grid') continue;
      try {
        const spec = s.find((x: ViewerSpec) => x.type === v.type);
        if (spec) {
          const col = spec.column ?? defaultCol;
          const entry = (map as any)[v.type];
          const prop = spec.prop ?? entry?.prop;
          const array = spec.array ?? entry?.array ?? false;
          if (prop) (v.props as any)[prop] = array ? [col] : col;
        }
        try { v.props.legendVisibility = 'Always'; } catch (_) {}
      } catch (_) {}
    }
    await new Promise((r) => setTimeout(r, settle));
  }, {s: specs, defaultCol: options.column, settle: settleMs, map: LEGEND_COLUMN_PROP});
}

export interface PickColumnOptions {

  comboboxSuffix: string;

  columnName: string;

  viewerType?: string;

  propName?: string;

  allowFallback?: boolean;

  popupWaitStrategy?: 'sleep' | 'backdrop' | 'either';

  scopeSelector?: string;
}

export async function pickColumnViaSelector(
  page: Page, opts: PickColumnOptions,
): Promise<{usedFallback: boolean}> {
  const selectorName = `div-column-combobox-${opts.comboboxSuffix}`;
  const scope = opts.scopeSelector ?? null;
  const strategy = opts.popupWaitStrategy ?? 'sleep';

  await page.evaluate(({name, sc}) => {
    document.body.dispatchEvent(new MouseEvent('mousedown', {bubbles: true}));
    const root: Document | Element = sc
      ? (document.querySelector(sc) as Element | null) ?? document
      : document;
    const sel = root.querySelector(`[name="${name}"]`);
    if (!sel) return;
    const colLabel = sel.querySelector('.d4-column-selector-column');
    (colLabel || sel).dispatchEvent(new MouseEvent('mousedown', {bubbles: true, button: 0}));
  }, {name: selectorName, sc: scope});
  if (strategy === 'sleep') {
    await page.waitForTimeout(500);
  } else if (strategy === 'backdrop') {
    await page.waitForFunction(() => !!document.querySelector('.d4-column-selector-backdrop'),
      null, {timeout: 3000}).catch(() => {});
  } else { 
    await Promise.race([
      page.waitForFunction(() => !!document.querySelector('.d4-column-selector-backdrop'),
        null, {timeout: 3000}).catch(() => {}),
      page.waitForTimeout(500),
    ]);
  }

  await page.keyboard.press(opts.columnName[0].toLowerCase());
  await page.waitForTimeout(100);
  if (opts.columnName.length > 1)
    await page.keyboard.type(opts.columnName.slice(1).toLowerCase());
  await page.keyboard.press('ArrowDown');
  await page.keyboard.press('Enter');
  await page.waitForTimeout(300);

  let usedFallback = false;
  if (opts.viewerType && opts.propName) {
    const applied = await page.evaluate(({vt, prop, col}) => {
      const view = (window as any).grok.shell.tv?.viewers?.find((x: any) => x.type === vt);
      return !!view && (view.props as any)[prop] === col;
    }, {vt: opts.viewerType, prop: opts.propName, col: opts.columnName});
    if (!applied && opts.allowFallback === true) {
      await page.evaluate(({vt, prop, col}) => {
        const view = (window as any).grok.shell.tv?.viewers?.find((x: any) => x.type === vt);
        if (view) (view.props as any)[prop] = col;
      }, {vt: opts.viewerType, prop: opts.propName, col: opts.columnName});
      usedFallback = true;
    }
  }
  return {usedFallback};
}

export interface TrustedPickColumnOptions {

  role: string;

  columnName: string;

  target?: 'auto' | 'column';

  viewerType?: string;

  propName?: string;

  scopeSelector?: string;

  backdropTimeoutMs?: number;

  commitSettleMs?: number;

  requirePopup?: boolean;
}

export async function pickColumnViaSelectorTrusted(
  page: Page, opts: TrustedPickColumnOptions,
): Promise<{popupOpened: boolean}> {
  const viewerType = opts.viewerType ?? 'Scatter plot';
  const rootName = `viewer-${viewerType.replace(/\s+/g, '-')}`;
  const target = opts.target ?? 'auto';
  const propName = opts.propName ?? `${opts.role}ColumnName`;

  const canvas = await page.evaluate((rn: string) => {
    const root = [...document.querySelectorAll(`[name="${rn}"]`)].find((e) => !e.closest('.d4-dialog'));
    const el = root?.querySelector('canvas[name="canvas"]') ?? root;
    if (!el) return null;
    const r = el.getBoundingClientRect();
    return {x: r.x + r.width / 2, y: r.y + r.height / 2};
  }, rootName);
  if (!canvas) throw new Error(`pickColumnViaSelectorTrusted: no ${viewerType} root on the page`);
  await page.mouse.move(canvas.x, canvas.y);

  const locate = ({rn, r, t, sc}: {rn: string; r: string; t: string; sc: string | null}) => {
    const root = sc
      ? document.querySelector(sc)
      : [...document.querySelectorAll(`[name="${rn}"]`)].find((e) => !e.closest('.d4-dialog'));
    const sel = root?.querySelector(`[name="div-column-combobox-${r}"]`) as HTMLElement | null;
    if (!sel) return null;
    const candidates: Element[] = [];
    const queries = t === 'column'
      ? ['.d4-column-selector-column']
      : ['.d4-column-selector-column', '.d4-column-selector-caption'];
    for (const q of queries) {
      const el = sel.querySelector(q);
      if (el) candidates.push(el);
    }
    if (t !== 'column') candidates.push(sel);
    for (const el of candidates) {
      const b = el.getBoundingClientRect();
      if (b.width > 0 && b.height > 0) return {x: b.x + b.width / 2, y: b.y + b.height / 2};
    }
    return null;
  };
  // the on-canvas selectors are hover-revealed, so the settle after mouse.move IS
  // "the selector became clickable" — which this lookup already reports
  const point = await pollValue(
    () => page.evaluate(locate, {rn: rootName, r: opts.role, t: target, sc: opts.scopeSelector ?? null}),
    (p) => p !== null, 400, 50);
  if (!point)
    throw new Error(`pickColumnViaSelectorTrusted: the ${opts.role} selector exposes no clickable text`);

  let popupOpened = false;
  for (let attempt = 0; attempt < 2 && !popupOpened; attempt++) {
    if (attempt > 0) await page.waitForTimeout(500);
    await page.mouse.click(point.x, point.y);
    popupOpened = await page.waitForFunction(
      () => !!document.querySelector('.d4-column-selector-backdrop'),
      null, {timeout: opts.backdropTimeoutMs ?? 6000}).then(() => true).catch(() => false);
  }
  if (!popupOpened) {
    if (opts.requirePopup === false) return {popupOpened: false};
    throw new Error(`pickColumnViaSelectorTrusted: the ${opts.role} column popup did not open`);
  }

  const comboItems = () => page.evaluate(() => document.querySelectorAll('.d4-combo-popup li').length);
  const text = opts.columnName.toLowerCase();
  await page.keyboard.press(text[0]);
  await pollValue(comboItems, (n) => n > 0, 150, 50);
  if (text.length > 1) await page.keyboard.type(text.slice(1));
  await pollValue(comboItems, (n) => n > 0, 200, 50);
  await page.keyboard.press('Enter');

  const readApplied = () => page.evaluate(({vt, prop}: {vt: string; prop: string}) => {
    const view = (window as any).grok.shell.tv?.viewers?.find((x: any) => x.type === vt);
    return view ? ((view.props as any)[prop] ?? null) : null;
  }, {vt: viewerType, prop: propName});
  const applied = await pollValue(readApplied, (a) => a === opts.columnName,
    opts.commitSettleMs ?? 900, 50);
  if (applied !== opts.columnName) {
    throw new Error(`pickColumnViaSelectorTrusted: ${opts.role} did not take — ` +
      `${propName} expected "${opts.columnName}", got "${applied}"`);
  }
  return {popupOpened: true};
}

export async function openViewerGear(page: Page, viewerType: string): Promise<void> {
  await page.evaluate((vt) => {
    const candidates = [
      `[name="viewer-${vt}"]`,
      `[name="viewer-${vt.replace(/\s+/g, '-')}"]`,
    ];
    // A closed table view leaves its viewers in the DOM at zero size, ahead of the
    // current view's. querySelector would pick that one and click a gear nobody can
    // see, so every property edit that followed landed on the wrong viewer.
    const sized = (e: Element) => {
      const b = e.getBoundingClientRect();
      return b.width > 0 && b.height > 0;
    };
    let vEl: HTMLElement | null = null;
    for (const c of candidates) {
      const all = [...document.querySelectorAll(c)] as HTMLElement[];
      vEl = all.find(sized) ?? all[0] ?? null;
      if (vEl) break;
    }
    if (!vEl) return;
    const panel = vEl.closest('.panel-base') as HTMLElement | null;
    const gear = panel?.querySelector('.panel-titlebar [name="icon-font-icon-settings"]') as HTMLElement | null;
    gear?.click();
  }, viewerType);
  // The property panel is a shared surface: switching viewers keeps .property-grid
  // in the DOM and even the same leading prop rows, so its presence says nothing
  // about WHICH viewer is showing. shell.o is what actually changes.
  await pollValue(() => page.evaluate((vt: string) => {
    const norm = (s: string) => s.replace(/[\s-]+/g, ' ').toLowerCase();
    const o = (window as any).grok?.shell?.o;
    return !!document.querySelector('.property-grid') && !!o?.type && norm(o.type) === norm(vt);
  }, viewerType), (ready) => ready, 1000, 50);
}

export async function clickViewerTitlebarIcon(
  page: Page, viewerName: string, iconName: string,
): Promise<void> {
  const point = await page.evaluate(({vn, icon}) => {
    const el = document.querySelector(`[name="viewer-${vn}"]`) as HTMLElement | null;
    const panel = el?.closest('.panel-base') as HTMLElement | null;
    const target = panel?.querySelector(`.panel-titlebar [name="${icon}"]`) as HTMLElement | null;
    if (!target) return null;
    const r = target.getBoundingClientRect();
    return {x: r.x + r.width / 2, y: r.y + r.height / 2};
  }, {vn: viewerName, icon: iconName});
  if (!point) throw new Error(`no ${iconName} icon on the ${viewerName} title bar`);
  await page.mouse.click(point.x, point.y);
}

export async function openViewerProperties(
  page: Page, viewerName: string, probeSelector = '.property-grid',
): Promise<void> {
  if (await page.locator(probeSelector).count() > 0) return;
  await clickViewerTitlebarIcon(page, viewerName, 'icon-font-icon-settings');
  await page.locator(probeSelector).first().waitFor({timeout: 10_000});
  await pollValue(() => page.locator(`${probeSelector} tr[name^="prop-"]`).count(), (n) => n > 0, 500, 50);
}

export async function ensurePropertyCategory(
  page: Page, viewerName: string, category: string, probeProp: string, timeoutMs = 20_000,
): Promise<void> {
  const headerSelector = `[name="prop-category-${category}"]`;
  const header = page.locator(headerSelector);

  const deadline = Date.now() + timeoutMs;
  while (Date.now() < deadline) {

    const index = await propertyRowIndex(page, probeProp, category);
    if (index >= 0 &&
        await page.locator(`.property-grid tr[name="prop-${probeProp}"]`).nth(index).isVisible())
      return;
    if (await header.count() === 0) {
      await clickViewerTitlebarIcon(page, viewerName, 'icon-font-icon-settings').catch(() => {});
      await page.locator(headerSelector).first().waitFor({timeout: 3000}).catch(() => {});
    } else
      await header.first().click().catch(() => {});
    await page.waitForTimeout(400);
  }
  throw new Error(`property-grid category "${category}" never exposed prop-${probeProp}`);
}

async function propertyRowIndex(page: Page, prop: string, category?: string): Promise<number> {
  if (!category) return 0;
  return page.evaluate(({p, c}) => {
    const rows = Array.from(document.querySelectorAll('.property-grid tr'));
    const start = rows.findIndex((r) => r.getAttribute('name') === `prop-category-${c}`);
    if (start < 0) return -1;
    let seen = 0;
    for (let i = 0; i < rows.length; i++) {
      const name = rows[i].getAttribute('name');
      if (name !== `prop-${p}`) continue;
      if (i > start && !rows.slice(start + 1, i).some((r) => r.className.includes('property-grid-category')))
        return seen;
      seen++;
    }
    return -1;
  }, {p: prop, c: category});
}

async function propertyRow(page: Page, prop: string, category?: string) {
  const index = await propertyRowIndex(page, prop, category);
  if (index < 0) throw new Error(`no prop-${prop} row inside category "${category}"`);
  return page.locator(`.property-grid tr[name="prop-${prop}"]`).nth(index);
}

export async function propertyGridValue(page: Page, prop: string, category?: string): Promise<string> {
  const index = await propertyRowIndex(page, prop, category);
  if (index < 0) return '';
  return page.evaluate(({p, i}) => {
    const row = document.querySelectorAll(`.property-grid tr[name="prop-${p}"]`)[i];
    if (!row) return '';
    const check = row.querySelector('input[type="checkbox"]') as HTMLInputElement | null;
    if (check) return String(check.checked);
    const view = row.querySelector(`[name="prop-view-${p}"]`) as HTMLElement | null;
    const input = row.querySelector('input:not([type="checkbox"])') as HTMLInputElement | null;
    const text = view?.innerText?.trim() ?? '';
    return text.length > 0 ? text : (input?.value ?? '');
  }, {p: prop, i: index});
}

export async function setPropertyGridValue(
  page: Page, prop: string, value: string, category?: string,
): Promise<void> {
  const row = await propertyRow(page, prop, category);
  await row.locator('td').last().click();
  await pollValue(() => row.locator('input:not([type="checkbox"])').count(), (n) => n > 0, 400, 50);
  await page.keyboard.press('Control+a');
  await page.keyboard.type(value);
  await page.keyboard.press('Enter');
  await pollValue(() => propertyGridValue(page, prop, category), (v) => v === value, 700, 50);
}

export async function selectPropertyGridChoice(
  page: Page, prop: string, value: string, category?: string,
): Promise<void> {
  const row = await propertyRow(page, prop, category);
  await row.locator('td').last().click();
  await row.locator('select').selectOption(value);
  // the cap, not an assertion: a choice whose display text differs from its option
  // value never settles equal, and this helper has never been the one to fail on it
  await pollValue(() => propertyGridValue(page, prop, category), (v) => v === value, 900, 50);
}

export async function setPropertyGridCheckbox(
  page: Page, prop: string, desired: boolean, category?: string,
): Promise<void> {
  const row = await propertyRow(page, prop, category);
  const box = row.locator('input[type="checkbox"]').first();
  if (await box.isChecked() === desired) return;
  await box.click();
  await pollValue(() => box.isChecked(), (c) => c === desired, 700, 50);
  expect(await box.isChecked()).toBe(desired);
}

export async function togglePropertyGridCheckbox(
  page: Page, prop: string, category?: string,
): Promise<boolean> {
  const row = await propertyRow(page, prop, category);
  const box = row.locator('input[type="checkbox"]').first();
  const before = await box.isChecked();
  await box.click();
  await pollValue(() => box.isChecked(), (c) => c !== before, 700, 50);
  return box.isChecked();
}

export async function driveTopMenuLeaf(
  page: Page, path: string[], opts: {attempts?: number} = {},
): Promise<boolean> {
  const attempts = opts.attempts ?? 4;
  const selAt = (depth: number) =>
    `[name="div-${path.slice(0, depth + 1).map((s) => s.replace(/ /g, '-')).join('---')}"]`;

  const liveCentre = async (sel: string): Promise<{x: number; y: number} | null> =>
    page.evaluate((s) => {
      const e = document.querySelector(s);
      if (!e) return null;
      const r = e.getBoundingClientRect();
      if (r.width === 0 || r.height === 0) return null;
      return {x: r.x + r.width / 2, y: r.y + r.height / 2};
    }, sel);
  const waitCentre = async (sel: string, ms: number): Promise<{x: number; y: number} | null> => {
    const deadline = Date.now() + ms;
    for (;;) {
      const c = await liveCentre(sel);
      if (c || Date.now() > deadline) return c;
      await page.waitForTimeout(150);
    }
  };

  for (let attempt = 0; attempt < attempts; attempt++) {
    const head = await liveCentre(selAt(0));
    if (head) await page.mouse.click(head.x, head.y);
    let ok = head != null;
    for (let depth = 1; ok && depth < path.length; depth++) {

      const c = await waitCentre(selAt(depth), 3000);
      if (!c) { ok = false; break; }
      if (depth < path.length - 1)
        await page.mouse.move(c.x, c.y, {steps: 6}); 
      else
        await page.mouse.click(c.x, c.y);            
    }
    if (ok) return true;
    await page.keyboard.press('Escape');
    await page.waitForTimeout(300);
  }
  return false;
}

export async function viewerSignature(page: Page, viewerName: string): Promise<string> {
  const shot = await page.locator(`[name="viewer-${viewerName}"]`).first().screenshot();
  return createHash('md5').update(shot).digest('hex');
}

export interface LegendInfo {

  itemCount: number;

  labels: string[];

  legendRendered: boolean;

  laidOut: boolean;
}

export async function readLegend(page: Page, viewerType: string): Promise<LegendInfo> {
  return await page.evaluate((vt) => {
    const tv = (window as any).grok.shell.tv;
    const v = tv?.viewers?.find((x: any) => x.type === vt);
    if (!v) return {itemCount: 0, labels: [], legendRendered: false, laidOut: false};
    const legendRoot = v.root.querySelector('[name="legend"]') as HTMLElement | null;
    const items = Array.from(v.root.querySelectorAll('[name="legend"] .d4-legend-item')) as HTMLElement[];
    const rect = legendRoot?.getBoundingClientRect();
    const laidOut = !!legendRoot && getComputedStyle(legendRoot).display !== 'none' &&
      legendRoot.offsetParent !== null && !!rect && rect.width > 0 && rect.height > 0;
    return {
      itemCount: items.length,
      labels: items.map((it) => (it.querySelector('.d4-legend-value')?.textContent ?? '').trim()),
      legendRendered: !!legendRoot,
      laidOut,
    };
  }, viewerType);
}

export interface RgbRange {
  rMin: number; rMax: number;
  gMin: number; gMax: number;
  bMin: number; bMax: number;
}

export interface CanvasPixelCounts {

  total: number;

  matched: number;
}

export async function countCanvasPixels(
  page: Page, viewerType: string, opts?: {rgbRange?: RgbRange; canvasSelector?: string},
): Promise<CanvasPixelCounts> {
  return await page.evaluate(({vt, range, cs}) => {
    try {
      const tv = (window as any).grok?.shell?.tv;
      const norm = (s: string) => s.replace(/[\s-]+/g, ' ').toLowerCase();
      const v = tv ? Array.from(tv.viewers).find((x: any) => norm(x.type) === norm(vt)) as any : null;
      const cv = v?.root?.querySelector(cs) as HTMLCanvasElement | null;
      const ctx = cv?.getContext('2d');
      if (!cv || !ctx) return {total: -1, matched: -1};
      const data = ctx.getImageData(0, 0, cv.width, cv.height).data;
      let total = 0, matched = 0;
      for (let i = 0; i < data.length; i += 4) {
        const r = data[i], g = data[i + 1], b = data[i + 2], a = data[i + 3];
        if (a === 0 || (r >= 250 && g >= 250 && b >= 250)) continue;
        total++;
        if (!range || (r >= range.rMin && r <= range.rMax &&
            g >= range.gMin && g <= range.gMax && b >= range.bMin && b <= range.bMax))
          matched++;
      }
      return {total, matched};
    } catch (_) {
      return {total: -1, matched: -1};
    }
  }, {vt: viewerType, range: opts?.rgbRange ?? null, cs: opts?.canvasSelector ?? 'canvas'});
}

export const SELECTION_HUE_RANGE: RgbRange = {rMin: 150, rMax: 255, gMin: 100, gMax: 200, bMin: 0, bMax: 110};

export async function countSelectionHuePixels(page: Page, viewerType: string): Promise<number> {
  return (await countCanvasPixels(page, viewerType, {rgbRange: SELECTION_HUE_RANGE})).matched;
}

export async function snapshotCanvasColors(
  page: Page, viewerType: string, canvasSelector = 'canvas',
): Promise<boolean> {
  return await page.evaluate(({vt, cs}) => {
    const w = window as any;
    const norm = (s: string) => s.replace(/[\s-]+/g, ' ').toLowerCase();
    const v = Array.from(w.grok?.shell?.tv?.viewers ?? [])
      .find((x: any) => norm(x.type) === norm(vt)) as any;
    const cv = v?.root?.querySelector(cs) as HTMLCanvasElement | null;
    if (!cv) return false;
    try {
      const img = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
      const colors = new Map<number, number>();
      for (let i = 0; i < img.length; i += 4) {
        const key = (img[i] << 16) | (img[i + 1] << 8) | img[i + 2];
        colors.set(key, (colors.get(key) ?? 0) + 1);
      }
      w.__canvasColorSnap = w.__canvasColorSnap || {};
      w.__canvasColorSnap[norm(vt)] = colors;
      return true;
    } catch {
      return false;
    }
  }, {vt: viewerType, cs: canvasSelector});
}

export async function diffCanvasColors(
  page: Page, viewerType: string, canvasSelector = 'canvas',
): Promise<{deltaPx: number}> {
  return await page.evaluate(({vt, cs}) => {
    const w = window as any;
    const norm = (s: string) => s.replace(/[\s-]+/g, ' ').toLowerCase();
    const prev = w.__canvasColorSnap?.[norm(vt)] as Map<number, number> | undefined;
    const v = Array.from(w.grok?.shell?.tv?.viewers ?? [])
      .find((x: any) => norm(x.type) === norm(vt)) as any;
    const cv = v?.root?.querySelector(cs) as HTMLCanvasElement | null;
    if (!cv || !prev) return {deltaPx: -1};
    try {
      const img = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
      const colors = new Map<number, number>();
      for (let i = 0; i < img.length; i += 4) {
        const key = (img[i] << 16) | (img[i + 1] << 8) | img[i + 2];
        colors.set(key, (colors.get(key) ?? 0) + 1);
      }
      let deltaPx = 0;
      for (const [c, n] of colors) deltaPx += Math.abs(n - (prev.get(c) ?? 0));
      for (const [c, n] of prev) if (!colors.has(c)) deltaPx += n;
      w.__canvasColorSnap[norm(vt)] = colors;
      return {deltaPx};
    } catch {
      return {deltaPx: -1};
    }
  }, {vt: viewerType, cs: canvasSelector});
}

export interface CanvasChangeOptions {

  minDelta?: number;

  canvasSelector?: string;

  timeoutMs?: number;
}

export async function waitForCanvasChange(
  page: Page, viewerType: string, opts: CanvasChangeOptions = {},
): Promise<number> {
  const minDelta = opts.minDelta ?? 1;
  const canvasSelector = opts.canvasSelector ?? 'canvas';
  await page.waitForFunction(({vt, cs, min}) => {
    const w = window as any;
    const norm = (s: string) => s.replace(/[\s-]+/g, ' ').toLowerCase();
    const prev = w.__canvasColorSnap?.[norm(vt)] as Map<number, number> | undefined;
    if (!prev) return false;
    const v = Array.from(w.grok?.shell?.tv?.viewers ?? [])
      .find((x: any) => norm(x.type) === norm(vt)) as any;
    const cv = v?.root?.querySelector(cs) as HTMLCanvasElement | null;
    if (!cv) return false;
    try {
      const img = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
      const colors = new Map<number, number>();
      for (let i = 0; i < img.length; i += 4) {
        const key = (img[i] << 16) | (img[i + 1] << 8) | img[i + 2];
        colors.set(key, (colors.get(key) ?? 0) + 1);
      }
      let delta = 0;
      for (const [c, n] of colors) delta += Math.abs(n - (prev.get(c) ?? 0));
      for (const [c, n] of prev) if (!colors.has(c)) delta += n;
      return delta >= min;
    } catch {
      return false;
    }
  }, {vt: viewerType, cs: canvasSelector, min: minDelta},
  {timeout: opts.timeoutMs ?? 15_000, polling: 250});
  return (await diffCanvasColors(page, viewerType, canvasSelector)).deltaPx;
}

/**
 * Waits until the viewer's canvas stops changing.
 *
 * `optional: true` reports a timeout as `false` instead of throwing — for a settle-precheck,
 * "still painting" is what the delta assertion downstream is there to catch, not a reason to
 * abort the step with a waitForFunction error. Default stays throwing: existing callers use it
 * as an assertion that the canvas DID settle.
 *
 * Returns true when the canvas went quiet.
 */
export async function waitForCanvasQuiet(
  page: Page, viewerType: string,
  opts: {canvasSelector?: string; stableReads?: number; timeoutMs?: number; optional?: boolean} = {},
): Promise<boolean> {
  await page.evaluate(() => { (window as any).__canvasQuiet = null; });
  const quiet = page.waitForFunction(({vt, cs, need}) => {
    const w = window as any;
    const norm = (s: string) => s.replace(/[\s-]+/g, ' ').toLowerCase();
    const v = Array.from(w.grok?.shell?.tv?.viewers ?? [])
      .find((x: any) => norm(x.type) === norm(vt)) as any;
    const cv = v?.root?.querySelector(cs) as HTMLCanvasElement | null;
    if (!cv) return false;
    try {
      const img = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
      let sig = 0;
      for (let i = 0; i < img.length; i += 4)
        sig = (sig * 31 + ((img[i] << 16) | (img[i + 1] << 8) | img[i + 2])) % 2147483647;
      const state = w.__canvasQuiet ?? {sig: null, count: 0};
      if (state.sig === sig) state.count++;
      else { state.sig = sig; state.count = 1; }
      w.__canvasQuiet = state;
      return state.count >= need;
    } catch {
      return false;
    }
  }, {vt: viewerType, cs: opts.canvasSelector ?? 'canvas', need: opts.stableReads ?? 2},
  {timeout: opts.timeoutMs ?? 20_000, polling: 300});
  if (!opts.optional) {
    await quiet;
    return true;
  }
  return quiet.then(() => true, () => false);
}

export async function waitForViewerRepaint(
  page: Page, viewerName: string, baseline: string, timeoutMs = 15_000,
): Promise<string> {
  await expect.poll(() => viewerSignature(page, viewerName),
    {timeout: timeoutMs, intervals: [200, 300, 500, 700, 1000]}).not.toBe(baseline);
  return viewerSignature(page, viewerName);
}

export async function waitForPropertyValue(
  page: Page, prop: string, expected: string, category?: string, timeoutMs = 10_000,
): Promise<void> {
  await expect.poll(() => propertyGridValue(page, prop, category),
    {timeout: timeoutMs, intervals: [100, 200, 300, 500]}).toBe(expected);
}

export interface ChangeColorOptions {

  viewerType: string;

  category: string;

  rgb: [number, number, number];

  hex: string;

  column: string;

  additive?: Record<string, string>;

  altContainerNames?: string[];
}

export async function changeLegendItemColor(page: Page, opts: ChangeColorOptions): Promise<void> {

  const dlgName = 'dialog-' + opts.category.replace(/[_\s]/g, '-');
  const containers = [
    `viewer-${opts.viewerType}`,
    `viewer-${opts.viewerType.replace(/\s+/g, '-')}`,
    ...(opts.altContainerNames ?? []),
  ];
  const rgbStr = `rgb(${opts.rgb[0]}, ${opts.rgb[1]}, ${opts.rgb[2]})`;
  const lowerHex = opts.hex.toLowerCase();
  const readTag = () => page.evaluate(({c, col}) => {
    const dfCol = (window as any).grok.shell.tv.dataFrame.col(col);
    const t = JSON.parse(dfCol.tags['.color-coding-categorical'] ?? '{}');
    return String(t[c] ?? '').toLowerCase();
  }, {c: opts.category, col: opts.column});
  const committed = (tag: string) => tag === lowerHex || tag.includes(lowerHex.slice(1));
  let okCommitted = false;
  try {
    let item = page.locator(`[name="${containers[0]}"] [name="legend"] .d4-legend-item`)
      .filter({has: page.locator('.d4-legend-value', {hasText: new RegExp(`^${opts.category}x?$`)})}).first();
    for (let i = 1; i < containers.length && (await item.count()) === 0; i++) {
      item = page.locator(`[name="${containers[i]}"] [name="legend"] .d4-legend-item`)
        .filter({has: page.locator('.d4-legend-value', {hasText: new RegExp(`^${opts.category}x?$`)})}).first();
    }
    await item.scrollIntoViewIfNeeded();
    await item.click({button: 'right', timeout: 5000});
    await page.locator(`.d4-dialog[name="${dlgName}"]`).waitFor({timeout: 5000});
    await page.evaluate(({dn, rgb}) => {
      const dlg = document.querySelector(`.d4-dialog[name="${dn}"]`)!;
      const sw = (Array.from(dlg.querySelectorAll('.d4-color-bar')) as HTMLElement[])
        .find((s) => s.style.backgroundColor === rgb);
      if (sw) {
        const o = {bubbles: true, cancelable: true, view: window, button: 0};
        sw.dispatchEvent(new MouseEvent('mousedown', o));
        sw.dispatchEvent(new MouseEvent('mouseup', o));
        sw.dispatchEvent(new MouseEvent('click', o));
      }
    }, {dn: dlgName, rgb: rgbStr});
    await page.waitForTimeout(200);
    await page.locator(`.d4-dialog[name="${dlgName}"] [name="button-OK"]`).click({timeout: 5000});
    okCommitted = committed(await pollValue(readTag, committed, 700, 50));
  } catch (_) {
    okCommitted = false;
  }
  if (!okCommitted) {
    const additive = opts.additive ?? {[opts.category]: opts.hex};
    await page.evaluate(({col, add}) => {
      const dfCol = (window as any).grok.shell.tv.dataFrame.col(col);
      const current: Record<string, string> = {};
      try { Object.assign(current, JSON.parse(dfCol.tags['.color-coding-categorical'] ?? '{}')); } catch (_) {}
      Object.assign(current, add);
      dfCol.tags['.color-coding-type'] = 'Categorical';
      dfCol.meta.colors.setCategorical(current, {fallbackColor: '#808080'});
      for (const v of (window as any).grok.shell.tv.viewers)
        if (v.type !== 'Grid') try { v.invalidate?.(); } catch (_) {}
    }, {col: opts.column, add: additive});
    await pollValue(readTag, (t) => t === lowerHex, 800, 50);
  }

  expect(await readTag()).toBe(lowerHex);
}

/**
 * Navigates an already-open `.d4-menu-popup` to `group > leaf` and clicks it.
 *
 * Shared by drivePanelMenuLeaf and driveContextMenuLeaf — the two differ only in how the popup
 * is opened. Kept as one function on purpose: the registry records eleven copy-pasted copies of
 * this block spreading through a single section, each inheriting the original's bugs and none of
 * its fixes.
 *
 * Failure modes: throws naming the missing segment and listing the labels actually visible.
 * Cleanup: N/A (leaves the menu closed by the click).
 */
async function navigateMenuPopup(page: Page, group: string | null, leaf: string): Promise<void> {
  await installEventWaits(page);
  await page.evaluate(({g, l}) => (window as any).__menuLeaf(g, l), {g: group, l: leaf});
}

/**
 * Opens a viewer's PANEL menu (the titlebar burger) and clicks `group > leaf`.
 *
 * Failure modes: throws when the burger, the group, or the leaf is absent.
 * Cleanup: N/A.
 */
export async function drivePanelMenuLeaf(
  page: Page, viewerName: string, group: string | null, leaf: string,
): Promise<void> {
  await page.evaluate((vn: string) => {
    const root = document.querySelector(`[name="viewer-${vn}"]`);
    const titlebar = root?.closest('.panel-base')?.querySelector('.panel-titlebar');
    const burger = titlebar?.querySelector('[name="icon-font-icon-menu"]') as HTMLElement | null;
    if (!burger)
      throw new Error(`no [name="icon-font-icon-menu"] in the .panel-titlebar of viewer-${vn}`);
    burger.click();
  }, viewerName);
  await page.locator('.d4-menu-popup').last().waitFor({timeout: 15_000});
  await navigateMenuPopup(page, group, leaf);
}

/**
 * Right-clicks a viewer's canvas and clicks `group > leaf` in the CONTEXT menu.
 *
 * The context menu is a different surface from the panel burger menu — different opener, same
 * popup structure — so drivePanelMenuLeaf cannot drive it. Specs were hand-rolling the
 * dispatch-hover-click sequence with fixed sleeps at every site instead.
 *
 * Failure modes: throws when the canvas, the group, or the leaf is absent.
 * Cleanup: N/A.
 */
export async function driveContextMenuLeaf(
  page: Page, viewerName: string, group: string | null, leaf: string,
  opts: {canvasSelector?: string} = {},
): Promise<void> {
  await page.evaluate(({vn, cs}) => {
    const canvas = (document.querySelector(`[name="viewer-${vn}"]`) as HTMLElement | null)
      ?.querySelector(cs) as HTMLCanvasElement | null;
    if (!canvas)
      throw new Error(`no ${cs} in viewer-${vn}`);
    const r = canvas.getBoundingClientRect();
    canvas.dispatchEvent(new MouseEvent('contextmenu', {
      bubbles: true, cancelable: true, button: 2,
      clientX: r.left + r.width / 2, clientY: r.top + r.height / 2,
    }));
  }, {vn: viewerName, cs: opts.canvasSelector ?? 'canvas'});
  await page.locator('.d4-menu-popup').last().waitFor({timeout: 15_000});
  await navigateMenuPopup(page, group, leaf);
}

const BAR_POSITIONS = (w: number, h: number, nCats: number) => [
  {x: w * 0.5, y: h * 0.2},
  {x: w * 0.5, y: h * 0.4},
  {x: w * 0.5, y: h * 0.6},
  {x: 60 + 0.5 * (w - 80) / Math.max(nCats, 1), y: h * 0.6},
  {x: w * 0.3, y: h * 0.3},
  {x: w * 0.7, y: h * 0.3},
];
const PIE_POSITIONS = (w: number, h: number) => [
  {x: w * 0.45, y: h * 0.4},
  {x: w * 0.5, y: h * 0.35},
  {x: w * 0.6, y: h * 0.5},
  {x: w * 0.4, y: h * 0.5},
  {x: w * 0.5, y: h * 0.6},
];

export interface CanvasFilterOptions {

  viewerType: 'Bar chart' | 'Pie chart';

  column: string;
}

export interface CanvasFilterResult {

  totalFiltered: number;

  survivors: number;

  canvasClickWorked: boolean;
}

export async function clickCanvasFilter(
  page: Page, opts: CanvasFilterOptions,
): Promise<CanvasFilterResult> {
  const setup = await page.evaluate(async ({vt, col}) => {
    const tv = (window as any).grok.shell.tv;
    const df = tv.dataFrame;
    df.filter.setAll(true);
    const fg = tv.getFiltersGroup();
    for (const f of Array.from(fg.filters as any)) { try { fg.remove(f); } catch (_) {} }
    const sp = tv.viewers.find((v: any) => v.type === 'Scatter plot');
    if (sp) try { sp.props.filter = ''; } catch (_) {}
    await new Promise((r) => setTimeout(r, 600));
    const v = tv.viewers.find((x: any) => x.type === vt);
    v.props.onClick = 'Filter';
    await new Promise((r) => setTimeout(r, 500));
    const cv = v.root.querySelector('canvas') as HTMLCanvasElement;
    const r = cv.getBoundingClientRect();
    const nCats = df.col(col)?.categories?.length ?? 0;
    return {w: r.width, h: r.height, before: df.filter.trueCount, nCats};
  }, {vt: opts.viewerType, col: opts.column});

  const positions = opts.viewerType === 'Bar chart'
    ? BAR_POSITIONS(setup.w, setup.h, setup.nCats)
    : PIE_POSITIONS(setup.w, setup.h);

  let canvasClickWorked = false;
  let result: CanvasFilterResult = {totalFiltered: setup.before, survivors: 0, canvasClickWorked: false};
  for (const pos of positions) {
    await page.evaluate(async () => {
      (window as any).grok.shell.tv.dataFrame.filter.setAll(true);
      await new Promise((r) => setTimeout(r, 300));
    });
    try {
      await page.locator(`[name="viewer-${opts.viewerType}"] canvas`).first()
        .click({position: {x: pos.x, y: pos.y}, timeout: 3000});
      await page.waitForTimeout(opts.viewerType === 'Bar chart' ? 800 : 900);
      const probe = await page.evaluate((col) => {
        const df = (window as any).grok.shell.tv.dataFrame;
        const c = df.col(col);
        const cats = c.categories;
        const counts: Record<string, number> = {};
        for (const x of cats) counts[x] = 0;
        for (let i = 0; i < df.rowCount; i++) if (df.filter.get(i)) counts[c.get(i)]++;
        const survivors = Object.entries(counts).filter(([_, n]) => n > 0).length;
        return {survivors, totalFiltered: df.filter.trueCount};
      }, opts.column);
      result = {...probe, canvasClickWorked: false};
      const ok = opts.viewerType === 'Bar chart'
        ? probe.survivors === 1 && probe.totalFiltered > 0
        : probe.totalFiltered !== setup.before && probe.totalFiltered > 0;
      if (ok) { canvasClickWorked = true; break; }
    } catch (_) {  }
  }
  if (!canvasClickWorked) {

    const probe = await page.evaluate(async (col) => {
      const tv = (window as any).grok.shell.tv;
      const df = tv.dataFrame;
      df.filter.setAll(true);
      const fg = tv.getFiltersGroup();
      for (const f of Array.from(fg.filters as any)) { try { fg.remove(f); } catch (_) {} }
      const DG = (window as any).DG;
      const cats = df.col(col).categories;
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: col, selected: [cats[0]]});
      await new Promise((r) => setTimeout(r, 1500));
      const c = df.col(col);
      const counts: Record<string, number> = {};
      for (const x of c.categories) counts[x] = 0;
      for (let i = 0; i < df.rowCount; i++) if (df.filter.get(i)) counts[c.get(i)]++;
      return {
        survivors: Object.entries(counts).filter(([_, n]) => n > 0).length,
        totalFiltered: df.filter.trueCount,
      };
    }, opts.column);
    result = {...probe, canvasClickWorked: false};
  } else {
    result.canvasClickWorked = true;
  }
  return result;
}

export async function applyCategoricalFilter(
  page: Page, column: string, selected: string[], settleMs = 1500,
): Promise<{filteredCount: number}> {
  return await page.evaluate(async ({col, sel, settle}) => {
    const tv = (window as any).grok.shell.tv;
    const df = tv.dataFrame;
    const fg = tv.getFiltersGroup();
    const DG = (window as any).DG;
    const settled = new Promise<void>((resolve) => {
      const sub = df.onRowsFiltered.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, settle);
    });
    fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: col, selected: sel});
    await settled;
    return {filteredCount: df.filter.trueCount};
  }, {col: column, sel: selected, settle: settleMs});
}

export async function applyNumericFilter(
  page: Page, column: string, min: number, max = 1e12, settleMs = 1500,
): Promise<number> {
  return await page.evaluate(async ({col, mn, mx, settle}) => {
    const tv = (window as any).grok.shell.tv;
    const df = tv.dataFrame;
    const fg = tv.getFiltersGroup();
    const settled = new Promise<void>((resolve) => {
      const sub = df.onRowsFiltered.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, settle);
    });
    fg.updateOrAdd({type: 'histogram', column: col, min: mn, max: mx});
    await settled;
    return df.filter.trueCount;
  }, {col: column, mn: min, mx: max, settle: settleMs});
}

export interface LayoutRoundTripResult {

  layoutId: string;

  filteredAfter: number;

  rowCountAfter: number;

  viewersAfter: number;
}

export async function saveAndReloadLayout(
  page: Page, namePrefix: string, applyCapMs = 3500,
): Promise<LayoutRoundTripResult> {
  return await page.evaluate(async ({prefix, cap}) => {
    const grok = (window as any).grok;
    const tv = grok.shell.tv;
    const layout = tv.saveLayout();
    layout.name = prefix + '_' + Date.now();
    const saved = await grok.dapi.layouts.save(layout);
    await new Promise((r) => setTimeout(r, 1000)); 
    const applied = new Promise<void>((resolve) => {
      const sub = grok.events.onViewLayoutApplied.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, cap);
    });
    tv.loadLayout(await grok.dapi.layouts.find(saved.id));
    await applied;

    for (let i = 0; i < 40; i++) {
      const t = grok.shell.tv;
      if (t?.dataFrame?.rowCount > 0 && t?.viewers?.length > 0) break;
      await new Promise((r) => setTimeout(r, 100));
    }
    const now = grok.shell.tv;
    return {
      layoutId: saved.id,
      filteredAfter: now.dataFrame.filter.trueCount,
      rowCountAfter: now.dataFrame.rowCount,
      viewersAfter: now.viewers.length,
    };
  }, {prefix: namePrefix, cap: applyCapMs});
}

export async function resetFilters(page: Page, opts: {clearScatterFilter?: boolean} = {}): Promise<void> {
  await page.evaluate(async (clearSp) => {
    const tv = (window as any).grok.shell.tv;
    if (!tv) return;
    const df = tv.dataFrame;
    const settled = new Promise<void>((resolve) => {
      const sub = df.onRowsFiltered.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 500);
    });
    df.filter.setAll(true);
    const fg = tv.getFiltersGroup();
    for (const f of Array.from(fg.filters as any)) { try { fg.remove(f); } catch (_) {} }
    if (clearSp) {
      const sp = tv.viewers.find((v: any) => v.type === 'Scatter plot');
      if (sp) try { sp.props.filter = ''; } catch (_) {}
    }
    await settled;
  }, opts.clearScatterFilter ?? false);
}

export async function cleanupShell(
  page: Page,
  opts: {clearStereoCategoryColorCoding?: boolean} = {},
): Promise<void> {
  await page.evaluate((clearColors) => {
    if (clearColors) {
      try {
        const col = (window as any).grok.shell.tv?.dataFrame.col('Stereo Category');
        if (col) {
          delete col.tags['.color-coding-categorical'];
          delete col.tags['.color-coding-type'];
        }
      } catch (_) {}
    }
    (window as any).grok.shell.closeAll();
  }, opts.clearStereoCategoryColorCoding ?? false);
  await page.waitForTimeout(500);
}

export interface ViewerPropStep {

  set: Record<string, any>;

  wait?: number;

  read?: string | string[];
}

export async function setViewerProps(
  page: Page, viewerType: string, steps: ViewerPropStep[], delayMs = 300,
): Promise<any[]> {
  return page.evaluate(async ({viewerType, steps, delayMs}) => {
    const h = Array.from((window as any).grok.shell.tv.viewers)
      .find((x: any) => x.type === viewerType) as any;
    const out: any[] = [];
    for (var step of steps) {

      const settled = new Promise<void>((resolve) => {
        let sub: any = null;
        try { sub = h.onViewerRendered.subscribe(() => { sub.unsubscribe(); resolve(); }); }
        catch (_) {  }
        setTimeout(() => { try { sub?.unsubscribe(); } catch (_) {} resolve(); }, step.wait ?? delayMs);
      });
      for (var k of Object.keys(step.set)) h.props[k] = step.set[k];
      await settled;
      if (step.read === undefined) continue;
      if (Array.isArray(step.read)) {
        const obj: Record<string, any> = {};
        for (var rk of step.read) obj[rk] = h.props[rk];
        out.push(obj);
      } else
        out.push(h.props[step.read]);
    }
    return out;
  }, {viewerType, steps, delayMs});
}

export async function waitForViewerRendered(page: Page, viewerType: string, capMs = 1000): Promise<void> {
  await page.evaluate(({type, capMs}) => new Promise<void>((resolve) => {
    const w = window as any;
    const t0 = Date.now();
    const v = Array.from(w.grok.shell.tv.viewers).find((x: any) => x.type === type) as any;
    if (!v) { setTimeout(resolve, 0); return; }

    w.__lastRender = w.__lastRender ?? {};
    if (!v.__renderStamped) {
      try {
        v.onViewerRendered.subscribe(() => { w.__lastRender[type] = Date.now(); });
        v.__renderStamped = true;
      } catch (_) {  }
    }
    const before = w.__lastRender[type] ?? 0;
    if (before && Date.now() - before < 400) { resolve(); return; }
    const tick = () => {
      if ((w.__lastRender[type] ?? 0) > before || Date.now() - t0 >= capMs) { resolve(); return; }
      setTimeout(tick, 25);
    };
    tick();
  }), {type: viewerType, capMs});
}

export async function waitForViewerEvent(
  page: Page, viewerType: string, eventProp: string, capMs = 3000,
): Promise<void> {
  await page.evaluate(({type, eventProp, capMs}) => new Promise<void>((resolve) => {
    const v = Array.from((window as any).grok.shell.tv.viewers).find((x: any) => x.type === type) as any;
    if (!v) { setTimeout(resolve, 0); return; }
    let sub: any = null;
    try { sub = v[eventProp].subscribe(() => { sub.unsubscribe(); resolve(); }); }
    catch (_) {  }
    setTimeout(() => { try { sub?.unsubscribe(); } catch (_) {} resolve(); }, capMs);
  }), {type: viewerType, eventProp, capMs});
}

export async function pollValue<T>(
  read: () => Promise<T>, ok: (value: T) => boolean, timeoutMs = 2500, intervalMs = 100,
): Promise<T> {
  const deadline = Date.now() + timeoutMs;
  let value = await read();
  while (!ok(value) && Date.now() < deadline) {
    await new Promise((res) => setTimeout(res, intervalMs));
    value = await read();
  }
  return value;
}

/**
 * Polls `read` until two consecutive reads compare equal under `same`, and returns that value.
 *
 * The stability idiom the specs kept hand-rolling: read, tick, re-read, stop when nothing moved.
 * Thirty-six copies of it existed across the viewer specs, each with its own iteration count and
 * its own fixed tick — which is what a "fixed sleep" in these files usually was. The tick lives
 * here now, in the helper layer, instead of at every call site.
 *
 * On timeout it returns the last value read rather than throwing: the caller's assertion is what
 * decides whether an unsettled value is a failure, and it gives a far better message than a bare
 * timeout would.
 *
 * Failure modes: none — always resolves. Cleanup: N/A.
 */
/**
 * Arms a JS-API event subscription, then hands back a function that awaits it.
 *
 * The two-phase twin of the in-page `__settled`, for when the ACTION is Playwright-side —
 * `page.mouse.move`, `page.keyboard.press`, `locator.click` — and so cannot be passed into the
 * browser as a callback:
 *
 *   const shown = await armEvent(page, 'grok.events.onTooltipShown', 2000);
 *   await page.mouse.move(x, y);          // the real gesture, driven from Playwright
 *   const args = await shown();           // event payload, or undefined if the cap won
 *
 * Arming BEFORE the action is the whole point: subscribe-after-act races the event and always
 * burns the cap. Without this, every Playwright-driven action had no rung-1 option at all, which
 * is why those sites ended up polling the DOM for a symptom instead of listening for the fact.
 *
 * Failure modes: an unknown channel rejects when awaited. Cleanup: unsubscribes on both paths.
 */
export async function armEvent(
  page: Page, channel: string, capMs = 2000,
): Promise<() => Promise<any>> {
  await installEventWaits(page);
  const token = `__armed_${Date.now()}_${Math.floor(Math.random() * 1e6)}`;
  await page.evaluate(({ch, cap, tok}) => {
    const w = window as any;
    w[tok] = w.__armed(ch, cap);
  }, {ch: channel, cap: capMs, tok: token});
  return async () => page.evaluate(async (tok) => {
    const w = window as any;
    const args = await w[tok];
    delete w[tok];
    return args ?? null;
  }, token);
}

export async function pollStable<T>(
  read: () => Promise<T>,
  same: (a: T, b: T) => boolean,
  timeoutMs = 3000,
  intervalMs = 150,
): Promise<T> {
  const deadline = Date.now() + timeoutMs;
  let prev = await read();
  while (Date.now() < deadline) {
    await new Promise((res) => setTimeout(res, intervalMs));
    const cur = await read();
    if (same(prev, cur)) return cur;
    prev = cur;
  }
  return prev;
}

export async function closeAllAndWait(page: Page): Promise<void> {
  await page.evaluate(() => (window as any).grok.shell.closeAll());
  await page.waitForFunction(() => Array.from((window as any).grok.shell.tableViews).length === 0,
    null, {timeout: 5000}).catch(() => {});
}

/**
 * Installs the in-page channel awaiter used inside `page.evaluate` bodies.
 *
 * Browser-context code cannot reach the Playwright-side wait helpers, which is why fixed sleeps
 * kept reappearing there. After this call a spec evaluates against named channels instead:
 *
 *   await w.__settled('df.onSelectionChanged', () => canvas.dispatchEvent(clickEvent));
 *   const args = await w.__eventFired('grok.events.onContextMenu', 2000);
 *
 * `__settled` subscribes BEFORE running the action and resolves on the first event, so it cannot
 * miss one fired synchronously; the cap is a ceiling, not a delay. Both resolve with the event
 * payload, or `undefined` when the cap wins.
 *
 * Channel grammar: `df.<stream>` (current dataframe), `grok.events.<stream>`,
 * `viewer:<Type>.<stream>` (viewer matched on its type, whitespace/dash-insensitive).
 *
 * Failure modes: an unknown channel throws inside the evaluate; a cap that is routinely consumed
 * means the channel does not fire for that action — a product gap to record, not a cap to raise.
 * Cleanup: N/A (idempotent, re-installed on navigation).
 */
export async function installEventWaits(page: Page): Promise<void> {
  await page.evaluate(() => {
    const w = window as any;
    if (w.__settled) return;

    const norm = (s: string) => s.replace(/[\s-]+/g, ' ').toLowerCase();

    w.__stream = (channel: string) => {
      const dot = channel.indexOf('.');
      const scope = channel.slice(0, dot);
      const rest = channel.slice(dot + 1);
      if (scope === 'df') return w.grok.shell.tv.dataFrame[rest];
      if (scope === 'grok') return w.grok.events[rest.replace(/^events\./, '')];
      if (scope.startsWith('viewer:')) {
        const type = scope.slice('viewer:'.length);
        const v = Array.from(w.grok.shell.tv.viewers).find((x: any) => norm(x.type) === norm(type)) as any;
        if (!v) throw new Error(`__settled: no viewer of type "${type}"`);
        return v[rest];
      }
      throw new Error(`__settled: unknown channel "${channel}"`);
    };

    w.__armed = (channel: string, capMs: number) => {
      const stream = w.__stream(channel);
      let sub: any = null;
      return new Promise<any>((resolve) => {
        sub = stream.subscribe((args: any) => { sub.unsubscribe(); resolve(args); });
        setTimeout(() => { try { sub?.unsubscribe(); } catch (_) {} resolve(undefined); }, capMs);
      });
    };

    w.__settled = async (channel: string, act: () => any, capMs = 2000) => {
      const fired = w.__armed(channel, capMs);
      await act();
      return fired;
    };

    w.__eventFired = (channel: string, capMs = 2000) => w.__armed(channel, capMs);

    // In-page twin of pollValue, for state with no event behind it — a viewer
    // appearing in tv.viewers after a project opens, a legend item rendering.
    // Resolves with the first value that satisfies `ok`, or the last one read.
    // A drag is a sequence of mousemove events PACED over time — the pacing is the gesture, not a
    // wait for something to settle, so there is no channel to race and no event to await mid-drag.
    // That makes it the one shape where a fixed sleep is correct, and it belongs here in the
    // helper rather than repeated at every call site.
    w.__drag = async (el: HTMLElement, from: {x: number; y: number}, to: {x: number; y: number},
      opts: {steps?: number; stepMs?: number; holdMs?: number; hoverFirst?: boolean;
        modifiers?: Record<string, boolean>} = {}) => {
      const {steps = 4, stepMs = 25, holdMs = 30, hoverFirst = false, modifiers = {}} = opts;
      const pause = (ms: number) => new Promise((r) => setTimeout(r, ms));
      const mk = (x: number, y: number) => ({
        bubbles: true, cancelable: true, clientX: x, clientY: y, button: 0, ...modifiers,
      });
      const at = (type: string, x: number, y: number, alsoDocument = false) => {
        el.dispatchEvent(new MouseEvent(type, mk(x, y)));
        if (alsoDocument) document.dispatchEvent(new MouseEvent(type, mk(x, y)));
      };
      if (hoverFirst) {
        at('mousemove', from.x, from.y);
        await pause(holdMs);
      }
      at('mousedown', from.x, from.y);
      await pause(holdMs);
      for (let i = 0; i <= steps; i++) {
        // steps: 0 is the degenerate path a hover-and-click uses; 0/0 would be NaN.
        const t = steps === 0 ? 1 : i / steps;
        at('mousemove', from.x + (to.x - from.x) * t, from.y + (to.y - from.y) * t, true);
        await pause(stepMs);
      }
      at('mouseup', to.x, to.y, true);
    };

    // Menu navigation, in-page so BOTH the Playwright helpers (drivePanelMenuLeaf /
    // driveContextMenuLeaf, which delegate here) and a spec whose surrounding logic sits inside
    // one evaluate use the same implementation. The alternative was a second copy in the browser
    // context, which is how eleven divergent copies of a menu driver got into one section.
    //
    // The view's top menu carries identically-labelled items and precedes every popup in document
    // order, so both segments are looked up inside the last-opened popup only. A panel menu is a
    // SINGLE popup — hovering a group appends no second container, the level lives in the element
    // name — so the same scope serves group leaves and top-level ones alike.
    w.__menuLeaf = async (group: string | null, leaf: string) => {
      const norm = (str: string | null | undefined) => (str ?? '').trim().toLowerCase();
      const labels = () => {
        const popup = [...document.querySelectorAll('.d4-menu-popup')].pop();
        return [...(popup?.querySelectorAll('.d4-menu-item-label') ?? [])];
      };
      const visible = () => labels().map((i) => (i.textContent ?? '').trim()).join(' | ');
      const find = (text: string) =>
        labels().find((i) => norm(i.textContent) === norm(text))?.closest('.d4-menu-item') ?? null;
      const waitFor = (text: string) => w.__poll(() => find(text), (i: Element | null) => i !== null, 5000);
      if (group !== null) {
        const groupItem = await waitFor(group);
        if (!groupItem)
          throw new Error(`menu: group "${group}" not found; visible: ${visible()}`);
        const b = groupItem.getBoundingClientRect();
        for (const type of ['mouseover', 'mousemove'])
          groupItem.dispatchEvent(new MouseEvent(type, {bubbles: true, clientX: b.x + 5, clientY: b.y + 5}));
      }
      const leafItem = await waitFor(leaf);
      if (!leafItem)
        throw new Error(
          `menu: leaf "${leaf}" not found in ${group === null ? 'the popup' : `the "${group}" submenu`}; `
          + `visible: ${visible()}`);
      (leafItem as HTMLElement).click();
      return true;
    };

    w.__poll = async (read: () => any, ok: (v: any) => boolean, timeoutMs = 3000, intervalMs = 100) => {
      const deadline = Date.now() + timeoutMs;
      let value = read();
      while (!ok(value) && Date.now() < deadline) {
        await new Promise((r) => setTimeout(r, intervalMs));
        value = read();
      }
      return value;
    };
  });
}

export function finishSpec(prefix = 'Step failures'): void {
  if (stepErrors.length === 0) return;
  const summary = stepErrors.map((e: StepError) => `- ${e.step}: ${e.error}`).join('\n');
  throw new Error(`${prefix}:\n${summary}`);
}
