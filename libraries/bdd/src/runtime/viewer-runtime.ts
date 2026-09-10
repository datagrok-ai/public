/* The in-page half of the viewer runtime: `window.__bdd`, installed once per page by
   `installViewerRuntime` (viewers.ts). Self-contained — `install` is serialized into the browser.
   Nothing here sleeps: a viewer says whether a render is pending (`isRenderPending`), announces
   every render (`onViewerRendered` / `onRendered`), and renders immediately on a bdd page
   (`immediateRendering`, set on every viewer the page holds or adds). */
import type {Page} from '@playwright/test';

declare const grok: any;
declare const DG: any;

export interface Box {
  x: number;
  y: number;
  width: number;
  height: number;
}

export interface Rect {
  x: number;
  y: number;
  w: number;
  h: number;
}

/** A viewer's canvas after a change: `delta` is the number of pixels that differ from the snapshot
 * taken before the change (the histogram distance when the bitmap changed size), `ink` the painted
 * pixels now, `inkBefore` those of the snapshot, `hue` and `hueBefore` the pixels in the selection
 * hue; `selected` rows and the `view` area in device pixels size the highlight a selection must
 * paint. */
export interface CanvasChange {
  delta: number;
  ink: number;
  inkBefore: number;
  hue: number;
  hueBefore: number;
  selected: number;
  viewPx: number;
  dpr: number;
}

export interface AreaChange {
  ink: number;
  inkBefore: number;
  rect: Rect;
}

export interface AreaColor {
  hex: string;
  count: number;
}

export interface AreaColors {
  colors: AreaColor[];
  rect: Rect;
  bitmap: number[];
}

export interface ScaleRange {
  min: number;
  max: number;
}

export interface ScaleChange {
  before?: ScaleRange;
  now?: ScaleRange;
}

export interface Range {
  top: number;
  bottom: number;
  height: number;
  left: number;
  right: number;
}

export interface RangeChange {
  before?: Range;
  now?: Range;
}

export interface Balloon {
  type: string;
  message: string;
}

export interface LegendState {
  mode: string;
  slot: string;
  items: number;
  keys: string[];
  selected: string[];
  width: number;
  height: number;
}

export interface Reading {
  before?: unknown;
  now?: unknown;
  has: string[];
}

/** A test over a column's values, evaluated in the page by `rowFacts` / `setRows`. */
export type RowTest = {eq: string} | {in: string[]} | {between: [number, number]} | {contains: string} |
  {startsWith: string} | {notNull: true};

/** How the rows a test names relate to the table's selection and filter. */
export interface RowFacts {
  matching: number;
  selected: number;
  passing: number;
  selectedTotal: number;
  filteredTotal: number;
  /** Rows whose selection bit disagrees with the test. */
  wrongSelection: number;
  wrongFilter: number;
}

function install(): void {
  const w = window as any;
  if (w.__bdd)
    return;
  interface Snapshot {
    bitmap: ImageData;
    scale: number;
    renders: number;
    range?: Range;
    rects: Record<string, Box>;
    scale$?: ScaleRange;
    values: Record<string, unknown>;
    legend?: LegendState;
    /** The histogram of the bitmap, computed on the first "than before" read. */
    stats?: {colors: Map<number, number>; ink: number; hue: number};
    areaInks: Record<string, number>;
  }
  const renders = new WeakMap<Element, {count: number; last: number; sub?: any}>();
  const snapshots = new WeakMap<Element, Snapshot>();
  const rememberedValues: Record<string, unknown> = {};
  // how long a repaint the viewer says is pending may take before that is a platform failure
  const PENDING_CAP = 10000;
  const sizes = new WeakMap<Element, {width: string; height: string}>();
  const forced = new WeakMap<HTMLElement, MutationObserver>();
  const listeners = new WeakMap<Element, Record<string, {count: number; sub: any}>>();
  const armed: Record<string, Promise<unknown>> = {};
  const balloons: Balloon[] = [];
  const remembered: Record<string, Range | undefined> = {};
  let layout: any;
  let tokens = 0;
  const norm = (s: unknown) => String(s ?? '').toLowerCase().replace(/[^a-z0-9]/g, '');
  // the platform's selected-rows orange, as the pixels of a marker or a box drawn in it
  const isHue = (r: number, g: number, b: number) => r >= 150 && g >= 100 && g <= 200 && b <= 110;
  const isBlank = (r: number, g: number, b: number, a: number) => a === 0 || (r >= 250 && g >= 250 && b >= 250);

  // --- the table ------------------------------------------------------------------------------------
  const table = (): any => {
    const t = grok.shell.t;
    if (!t)
      throw new Error('no current table');
    return t;
  };
  const tableNamed = (name: string): any => {
    const t = (grok.shell.tables as any[]).find((x) => x.name === name);
    if (!t)
      throw new Error(`table "${name}" is not open; open: ${grok.shell.tables.map((x: any) => x.name).join(' | ')}`);
    return t;
  };
  const col = (name: string, df: any = table()): any => {
    const c = df.col(name);
    if (!c)
      throw new Error(`no "${name}" column in ${df.name}; it has: ${df.columns.names().join(', ')}`);
    return c;
  };
  const tester = (c: any, test: RowTest): ((i: number) => boolean) => {
    const text = (i: number) => String(c.get(i) ?? '');
    if ('eq' in test)
      return (i) => text(i) === test.eq;
    if ('in' in test)
      return (i) => test.in.includes(text(i));
    if ('contains' in test)
      return (i) => text(i).includes(test.contains);
    if ('startsWith' in test)
      return (i) => text(i).startsWith(test.startsWith);
    if ('notNull' in test)
      return (i) => !c.isNone(i);
    const [lo, hi] = test.between;
    return (i) => { const v = c.get(i); return v != null && !isNaN(v) && v >= lo && v <= hi; };
  };
  const describeTest = (name: string, test: RowTest): string =>
    'eq' in test ? `${name} = ${test.eq}` : 'in' in test ? `${name} in ${test.in.join(', ')}` :
      'contains' in test ? `${name} containing "${test.contains}"` : 'startsWith' in test ? `${name} starting with "${test.startsWith}"` :
        'notNull' in test ? `${name} not null` : `${name} between ${test.between[0]} and ${test.between[1]}`;
  /** A test that could never name a row is a mistake in the feature, not a fact about the data: a
   * category the column does not hold (a typo), a range over a text column. Checked before the
   * scan, so that an empty result stays legal — a range that keeps no row empties a chart on
   * purpose, and a claim about it is a claim about 0 rows. */
  const checkTest = (name: string, c: any, test: RowTest): void => {
    const wanted = 'eq' in test ? [test.eq] : 'in' in test ? test.in : [];
    if (wanted.length > 0) {
      const have = new Set<string>();
      for (let i = 0; i < c.length; i++)
        have.add(String(c.get(i) ?? ''));
      const missing = wanted.filter((v) => !have.has(v));
      if (missing.length > 0)
        throw new Error(`${name} has no value ${missing.map((v) => `"${v}"`).join(', ')}; it has: ${[...have].slice(0, 20).join(', ')}${have.size > 20 ? ', …' : ''}`);
    }
    if ('between' in test && !c.isNumerical)
      throw new Error(`${name} is a ${c.type} column; a range test needs a numerical one`);
  };
  /** The rows a test names against the selection and the filter. A claim about rows none of
   * which exist fails: "only rows where X starts with Z" over a table with no such row would be
   * checking that nothing is selected. */
  const rowFacts = (name: string, test: RowTest): RowFacts => {
    const df = table();
    const c = col(name, df);
    checkTest(name, c, test);
    const hit = tester(c, test);
    const facts: RowFacts = {matching: 0, selected: 0, passing: 0, selectedTotal: df.selection.trueCount, filteredTotal: df.filter.trueCount,
      wrongSelection: 0, wrongFilter: 0};
    for (let i = 0; i < df.rowCount; i++) {
      const h = hit(i);
      const s = df.selection.get(i);
      const f = df.filter.get(i);
      if (h) {
        facts.matching++;
        if (s)
          facts.selected++;
        if (f)
          facts.passing++;
      }
      if (h !== s)
        facts.wrongSelection++;
      if (h !== f)
        facts.wrongFilter++;
    }
    if (facts.matching === 0)
      throw new Error(`no row of ${df.name} has ${describeTest(name, test)}`);
    return facts;
  };
  /** The selection or the filter becomes exactly the rows the test names (or every other row);
   * an empty result is a legal one. */
  const setRows = (what: 'selection' | 'filter', name: string, test: RowTest, negate = false): number => {
    const df = table();
    const c = col(name, df);
    checkTest(name, c, test);
    const hit = tester(c, test);
    let hits = 0;
    df[what].init((i: number) => {
      const h = hit(i) !== negate;
      if (h)
        hits++;
      return h;
    });
    return hits;
  };

  // --- viewers ----------------------------------------------------------------------------------------
  const viewers = (): any[] => {
    const all: any[] = [];
    for (const view of Array.from(grok.shell.tableViews ?? []) as any[]) {
      for (const v of Array.from(view.viewers ?? []) as any[])
        all.push(v);
    }
    return all;
  };
  const findViewer = (el: Element): any => viewers().find((x) => x.root === el || x.root.contains(el) || el.contains(x.root));
  const viewerOf = (el: Element): any => {
    const v = findViewer(el);
    if (!v)
      throw new Error('the element is not a viewer of an open table view');
    return v;
  };
  const arm = (v: any): void => {
    if (!v || renders.has(v.root))
      return;
    const stamp: {count: number; last: number; sub?: any} = {count: 0, last: 0};
    renders.set(v.root, stamp);
    try {
      v.immediateRendering = true;
    } catch { /* a JS viewer without the flag */ }
    // a JS viewer announces its own renders on `onRendered` (the host's onViewerRendered never fires for it)
    try {
      stamp.sub = (v.onRendered ?? v.onViewerRendered).subscribe(() => {
        stamp.count++;
        stamp.last = Date.now();
      });
    } catch { /* a viewer without the event */ }
  };
  const stampAll = (): void => {
    for (const v of viewers())
      arm(v);
  };
  /** Whether the viewer says a change is on its way to its canvas; undefined for a viewer without
   * the signal (a JS viewer). */
  const pending = (v: any): boolean | undefined => {
    try {
      const p = v.isRenderPending;
      return typeof p === 'boolean' ? p : undefined;
    }
    catch {
      return undefined;
    }
  };
  /** Resolves with the renders since the call once the viewer has nothing on its way — through
   * every render it announces (a legend change re-anchors in up to four passes 100 ms apart). A
   * pending one that never lands is a platform failure, reported after `PENDING_CAP`. A viewer
   * without the signal falls back to a render event within `capMs`. */
  const settle = (el: Element, capMs: number): Promise<number> => {
    const v = viewerOf(el);
    arm(v);
    const stamp = renders.get(v.root)!;
    const before = stamp.count;
    const t0 = Date.now();
    return new Promise((resolve, reject) => {
      const tick = () => {
        const p = pending(v);
        if (stamp.count > before && p !== true)
          return resolve(stamp.count - before);
        if (p === false)
          return resolve(0);
        if (p === undefined) {
          if (Date.now() - t0 >= capMs)
            return resolve(0);
          return setTimeout(tick, 10);
        }
        if (Date.now() - t0 >= PENDING_CAP)
          return reject(new Error(`${v.type}: a repaint has been pending for ${PENDING_CAP} ms`));
        setTimeout(tick, 0);
      };
      setTimeout(tick, 0);
    });
  };
  /** After the frame a change would land on, and until the viewer has nothing pending. */
  const quiet = async (v: any): Promise<void> => {
    await new Promise((resolve) => requestAnimationFrame(() => setTimeout(resolve, 0)));
    const t0 = Date.now();
    while (pending(v) === true) {
      if (Date.now() - t0 >= PENDING_CAP)
        throw new Error(`${v.type}: a repaint has been pending for ${PENDING_CAP} ms`);
      await new Promise((resolve) => setTimeout(resolve, 0));
    }
  };
  const partsOf = (v: any): Record<string, Element> => v.getWidgetStatus?.()?.parts ?? {};
  const canvasOf = (v: any): HTMLCanvasElement => {
    const part = partsOf(v).canvas as HTMLCanvasElement | undefined;
    const cv = part ?? v.root.querySelector('canvas[name="canvas"]') ?? v.root.querySelector('canvas');
    if (!cv)
      throw new Error(`${v.type} has no canvas`);
    return cv;
  };
  /** What a viewer's hit areas are relative to: its canvas, else the root it reports, else its
   * root (a form, a filter panel — DOM viewers with no canvas). */
  const anchorOf = (v: any): Element => {
    const parts = partsOf(v);
    return parts.canvas ?? parts.root ?? v.root;
  };
  // a WebGL canvas (the 3D scatter plot) has no 2D context and no pixels to read: its picture is
  // the viewer's own `scene signature` reading
  const pixels = (cv: HTMLCanvasElement): ImageData => {
    const ctx = cv.getContext('2d');
    return ctx ? ctx.getImageData(0, 0, cv.width, cv.height) : new ImageData(1, 1);
  };
  /** The viewer's picture: its canvas with its `overlay` part composited on top when there is one
   * of the same size (the scatter plot draws regression lines, labels and stats on the overlay,
   * the grid its selection and current cell). */
  const pixelsOf = (v: any): ImageData => {
    const cv = canvasOf(v);
    const img = pixels(cv);
    const over = partsOf(v).overlay as HTMLCanvasElement | undefined;
    if (!(over instanceof HTMLCanvasElement) || over === cv || over.width !== cv.width || over.height !== cv.height)
      return img;
    const top = pixels(over);
    const a = img.data;
    const b = top.data;
    for (let i = 0; i < b.length; i += 4) {
      const alpha = b[i + 3];
      if (alpha === 0)
        continue;
      const w = alpha / 255;
      a[i] = Math.round(b[i] * w + a[i] * (1 - w));
      a[i + 1] = Math.round(b[i + 1] * w + a[i + 1] * (1 - w));
      a[i + 2] = Math.round(b[i + 2] * w + a[i + 2] * (1 - w));
      a[i + 3] = Math.max(a[i + 3], alpha);
    }
    return img;
  };
  const histogram = (img: ImageData) => {
    const data = img.data;
    const colors = new Map<number, number>();
    let ink = 0;
    let hue = 0;
    for (let i = 0; i < data.length; i += 4) {
      const key = (data[i] << 16) | (data[i + 1] << 8) | data[i + 2];
      colors.set(key, (colors.get(key) ?? 0) + 1);
      if (isBlank(data[i], data[i + 1], data[i + 2], data[i + 3]))
        continue;
      ink++;
      if (isHue(data[i], data[i + 1], data[i + 2]))
        hue++;
    }
    return {colors, ink, hue};
  };
  const areasOf = (v: any): Record<string, Box> => v.getWidgetStatus()?.hitAreas ?? {};
  const areaKey = (v: any, name: string): string => {
    const areas = areasOf(v);
    const key = Object.keys(areas).find((k) => norm(k) === norm(name));
    if (!key)
      throw new Error(`${v.type} has no "${name}" area right now; it has: ${Object.keys(areas).join(', ') || 'none'}`);
    return key;
  };
  const scaleOf = (cv: HTMLCanvasElement): number => cv.width / cv.getBoundingClientRect().width;
  /** A hit area (canvas coordinates) as device pixels of the canvas bitmap. */
  const deviceRect = (r: Box, scale: number): Rect =>
    ({x: Math.floor(r.x * scale), y: Math.floor(r.y * scale), w: Math.max(1, Math.floor(r.width * scale)), h: Math.max(1, Math.floor(r.height * scale))});
  /** The colors drawn inside a rectangle of the bitmap, by pixel count, blanks aside. */
  const colorsIn = (img: ImageData, r: Rect): Map<number, number> => {
    const colors = new Map<number, number>();
    const data = img.data;
    for (let y = r.y; y < Math.min(r.y + r.h, img.height); y++) {
      for (let x = r.x; x < Math.min(r.x + r.w, img.width); x++) {
        const i = (y * img.width + x) * 4;
        if (isBlank(data[i], data[i + 1], data[i + 2], data[i + 3]))
          continue;
        const key = (data[i] << 16) | (data[i + 1] << 8) | data[i + 2];
        colors.set(key, (colors.get(key) ?? 0) + 1);
      }
    }
    return colors;
  };
  const inkIn = (img: ImageData, r: Rect): number => {
    let n = 0;
    const data = img.data;
    for (let y = r.y; y < Math.min(r.y + r.h, img.height); y++) {
      for (let x = r.x; x < Math.min(r.x + r.w, img.width); x++) {
        const i = (y * img.width + x) * 4;
        if (!isBlank(data[i], data[i + 1], data[i + 2], data[i + 3]))
          n++;
      }
    }
    return n;
  };
  /** The colors drawn in at least `minPx` pixels, blanks and near-whites aside. */
  const palette = (el: Element, minPx: number): number => {
    const data = histogram(pixelsOf(viewerOf(el)));
    let n = 0;
    for (const [c, count] of data.colors) {
      const r = (c >> 16) & 255;
      const g = (c >> 8) & 255;
      const b = c & 255;
      if (count >= minPx && !(r >= 250 && g >= 250 && b >= 250) && c !== 0)
        n++;
    }
    return n;
  };
  const hex = (c: number): string => '#' + c.toString(16).padStart(6, '0').toUpperCase();
  /** The colors drawn inside a hit area, most pixels first, with the bitmap rectangle read. */
  const areaColors = (el: Element, name: string): AreaColors => {
    const v = viewerOf(el);
    const cv = canvasOf(v);
    const rect = deviceRect(areasOf(v)[areaKey(v, name)], scaleOf(cv));
    const colors = colorsIn(pixelsOf(v), rect);
    return {colors: [...colors.entries()].sort((a, b) => b[1] - a[1]).map(([c, count]) => ({hex: hex(c), count})), rect, bitmap: [cv.width, cv.height]};
  };
  /** The legend a viewer hosts, as the attributes it publishes on every commit. A legend shown in
   * the tooltip is re-parented out of the viewer, so it is looked up there. */
  const legendState = (v: any): LegendState | undefined => {
    const l: HTMLElement | null = v.root.querySelector('[name="legend"]') ?? document.querySelector('.d4-tooltip [name="legend"]');
    if (!l)
      return undefined;
    const items = Array.from(l.querySelectorAll('[name="legend-item"]'));
    const key = (i: Element) => i.getAttribute('data-item-key') ?? i.getAttribute('aria-label') ?? '';
    return {mode: l.dataset.legendMode ?? '', slot: l.dataset.legendSlot ?? '', items: Number(l.dataset.legendItems ?? items.length),
      keys: items.map(key), selected: items.filter((i) => i.getAttribute('aria-selected') === 'true' || i.getAttribute('data-item-selected') === 'true').map(key),
      width: l.offsetWidth, height: l.offsetHeight};
  };
  const legendChange = (el: Element): {before?: LegendState; now?: LegendState} => {
    const v = viewerOf(el);
    return {before: snapshots.get(v.root)?.legend, now: legendState(v)};
  };
  const valuesOf = (v: any): Record<string, unknown> => v.getWidgetStatus?.()?.values ?? {};
  /** A named reading now and at the snapshot before the last change. */
  const valueChange = (el: Element, name: string): Reading => {
    const v = viewerOf(el);
    const now = valuesOf(v);
    const key = Object.keys(now).find((k) => norm(k) === norm(name));
    const before = snapshots.get(v.root)?.values ?? {};
    const keyBefore = Object.keys(before).find((k) => norm(k) === norm(name));
    return {before: keyBefore === undefined ? undefined : before[keyBefore], now: key === undefined ? undefined : now[key], has: Object.keys(now)};
  };
  /** A reading kept by viewer type and name, for "as remembered" across a close and a reopen. */
  const rememberValue = (el: Element, name: string): void => {
    const v = viewerOf(el);
    const r = valueChange(el, name);
    if (r.now === undefined || r.now === null)
      throw new Error(`${v.type} has no "${name}" reading; it reports: ${r.has.join(', ') || 'no readings'}`);
    rememberedValues[`${v.type}|${norm(name)}`] = r.now;
  };
  const rememberedValue = (el: Element, name: string): Reading => {
    const v = viewerOf(el);
    const r = valueChange(el, name);
    return {before: rememberedValues[`${v.type}|${norm(name)}`], now: r.now, has: r.has};
  };
  /** A hit area's rectangle now and at the snapshot before the last change (anchor coordinates). */
  const areaRectChange = (el: Element, name: string): {before?: Box; now?: Box; has: string[]} => {
    const v = viewerOf(el);
    const areas = areasOf(v);
    const key = Object.keys(areas).find((k) => norm(k) === norm(name));
    const before = snapshots.get(v.root)?.rects ?? {};
    return {before: before[norm(name)], now: key === undefined ? undefined : areas[key], has: Object.keys(areas)};
  };
  const scaleRange = (v: any): ScaleRange | undefined => {
    const values = valuesOf(v);
    const min = values['color scale min'];
    const max = values['color scale max'];
    return typeof min === 'number' && typeof max === 'number' ? {min, max} : undefined;
  };
  const tableOf = (el: Element): string => String(viewerOf(el).dataFrame?.name ?? '');
  const rangeOf = (v: any): Range | undefined => {
    try {
      const vp = v.viewport;
      return vp ? {top: vp.top, bottom: vp.bottom, height: vp.height, left: vp.left, right: vp.right} : undefined;
    }
    catch {
      return undefined;
    }
  };
  const property = (v: any, caption: string): any => {
    const want = norm(caption);
    const short = norm(caption.replace(/\s+columns?$/i, ''));
    const friendly = (name: string) => name.replace(/ColumnNames?$/, '');
    const props: any[] = v.getProperties();
    // "Marker Size" is the markerSize property; "Marker Size Column" is markerSizeColumnName
    const p = props.find((x) => norm(x.caption) === want) ?? props.find((x) => norm(x.name) === want) ??
      props.find((x) => norm(friendly(x.name)) === want) ??
      props.find((x) => /ColumnNames?$/.test(x.name) && norm(friendly(x.name)) === short);
    if (!p) {
      // a property's caption defaults to its name; the grid shows the name split into words
      const title = (x: any): string => x.caption && x.caption !== x.name ? x.caption :
        x.name.replace(/ColumnNames?$/, ' Column').replace(/[a-z\d](?=[A-Z])|[a-z](?=\d)|[A-Z](?=[A-Z][a-z])/g, '$& ')
          .replace(/^./, (c: string) => c.toUpperCase());
      const distance = (a: string, b: string): number => {
        const row = [...Array(b.length + 1).keys()];
        for (let i = 1; i <= a.length; i++) {
          let prev = row[0]++;
          for (let j = 1; j <= b.length; j++) {
            const cur = row[j];
            row[j] = Math.min(row[j] + 1, row[j - 1] + 1, prev + (a[i - 1] === b[j - 1] ? 0 : 1));
            prev = cur;
          }
        }
        return row[b.length];
      };
      const titles = props.map(title);
      const nearest = titles.map((t) => ({t, d: distance(norm(t), want)})).sort((x, y) => x.d - y.d).slice(0, 3).map((x) => x.t);
      throw new Error(`${v.type} has no "${caption}" property; nearest: ${nearest.join(', ')}. All: ${titles.join(', ')}`);
    }
    return p;
  };
  const HEX = /^#?([0-9a-f]{6})$/i;
  const convert = (p: any, text: string): unknown => {
    const s = text.replace(/\\n/g, '\n');
    const type = String(p.propertyType ?? '');
    if (type === 'bool')
      return /^(true|yes|on|checked|1)$/i.test(s);
    if (type === 'int' || type === 'double' || type === 'num' || type === 'bigint') {
      if (s === '')
        return null;
      const hex = HEX.exec(s);
      if (hex && type === 'int')
        return (0xFF000000 | parseInt(hex[1], 16)) >>> 0;
      return Number(s);
    }
    if (type === 'string_list' || type === 'list' || type === 'column_list')
      return s === '' ? [] : s.split(/\s*,\s*/);
    return s;
  };
  /** The property as text: "" for none, lists comma-joined, an int read against a `#rrggbb`
   * expectation as `#rrggbb` (colors are ints in the property bag). */
  const readProperty = (el: Element, caption: string, expected = ''): string => {
    const v = viewerOf(el);
    const p = property(v, caption);
    const value = v.props[p.name];
    if (value == null)
      return '';
    if (Array.isArray(value))
      return value.join(', ');
    if (typeof value === 'number' && String(p.propertyType) === 'int' && HEX.test(expected))
      return '#' + (value & 0xFFFFFF).toString(16).padStart(6, '0').toUpperCase();
    return String(value);
  };
  /** The state the "than before" checks compare with: the bitmap, the hit areas, the value range,
   * the color scale, the readings and the legend. The bitmap's histogram and the ink of every area
   * are computed on the first read that needs them — most changes are never asked "than before",
   * and a histogram of a million pixels per property set was the cost of that. */
  const snapshot = (el: Element): void => {
    const v = viewerOf(el);
    arm(v);
    const hit = areasOf(v);
    const rects: Record<string, Box> = {};
    for (const key of Object.keys(hit))
      rects[norm(key)] = {...hit[key]};
    let cv: HTMLCanvasElement | undefined;
    try {
      cv = canvasOf(v);
    } catch { /* a DOM viewer: its areas, readings and legend are the snapshot */ }
    snapshots.set(v.root, {bitmap: cv ? pixelsOf(v) : new ImageData(1, 1), scale: cv ? scaleOf(cv) : 1, renders: renders.get(v.root)!.count,
      range: rangeOf(v), rects, scale$: scaleRange(v), values: {...valuesOf(v)}, legend: legendState(v), areaInks: {}});
  };
  const snapshotOf = (v: any): Snapshot => {
    const s = snapshots.get(v.root);
    if (!s)
      throw new Error(`${v.type}: no snapshot to compare with`);
    return s;
  };
  const statsOf = (s: Snapshot) => s.stats ??= histogram(s.bitmap);
  const areaInkBefore = (s: Snapshot, key: string): number | undefined => {
    const rect = s.rects[key];
    if (!rect)
      return undefined;
    return s.areaInks[key] ??= inkIn(s.bitmap, deviceRect(rect, s.scale));
  };
  /** Pixels that differ between two bitmaps of the same size — a reorder of equal bars keeps the
   * color histogram and is still a repaint; bitmaps of different sizes differ by their histogram
   * and their size. */
  const pixelDelta = (before: Snapshot, now: ImageData, colorsNow: Map<number, number>): number => {
    if (before.bitmap.width !== now.width || before.bitmap.height !== now.height) {
      const colorsBefore = statsOf(before).colors;
      let delta = Math.abs(before.bitmap.width * before.bitmap.height - now.width * now.height);
      for (const [c, n] of colorsNow)
        delta += Math.abs(n - (colorsBefore.get(c) ?? 0));
      for (const [c, n] of colorsBefore) {
        if (!colorsNow.has(c))
          delta += n;
      }
      return delta;
    }
    const a = before.bitmap.data;
    const b = now.data;
    let delta = 0;
    for (let i = 0; i < a.length; i += 4) {
      if (a[i] !== b[i] || a[i + 1] !== b[i + 1] || a[i + 2] !== b[i + 2] || a[i + 3] !== b[i + 3])
        delta++;
    }
    return delta;
  };
  /** Where the picture differs from the snapshot, for the failure to say: the changed pixels'
   * box in CSS px and the hit areas it touches, or the two canvas sizes when they differ. */
  const changedWhere = (v: any, cv: HTMLCanvasElement, before: ImageData, now: ImageData): string => {
    if (before.width !== now.width || before.height !== now.height)
      return `canvas ${before.width}x${before.height} → ${now.width}x${now.height}`;
    const a = before.data;
    const b = now.data;
    let x0 = now.width, y0 = now.height, x1 = -1, y1 = -1;
    for (let i = 0, p = 0; i < a.length; i += 4, p++) {
      if (a[i] !== b[i] || a[i + 1] !== b[i + 1] || a[i + 2] !== b[i + 2] || a[i + 3] !== b[i + 3]) {
        const x = p % now.width;
        const y = (p - x) / now.width;
        if (x < x0) x0 = x;
        if (x > x1) x1 = x;
        if (y < y0) y0 = y;
        if (y > y1) y1 = y;
      }
    }
    if (x1 < 0)
      return '';
    const s = scaleOf(cv);
    const box = {x: x0 / s, y: y0 / s, width: (x1 - x0 + 1) / s, height: (y1 - y0 + 1) / s};
    const touched = Object.entries(areasOf(v)).filter(([n, r]) => n !== 'view' && r.x < box.x + box.width && box.x < r.x + r.width && r.y < box.y + box.height && box.y < r.y + r.height).map(([n]) => n);
    return `within ${Math.round(box.x)},${Math.round(box.y)} ${Math.round(box.width)}x${Math.round(box.height)} of the canvas`
      + (touched.length ? `, over: ${touched.join(', ')}` : '');
  };
  /** The painted pixels now, without moving the snapshot. */
  const ink = (el: Element): number => histogram(pixelsOf(viewerOf(el))).ink;
  const baseline = (el: Element): void => {
    try {
      snapshot(el);
    } catch { /* not a viewer of an open view */ }
  };
  /** Every viewer's baseline at once — before a change that reaches them all (a filter, a
   * selection, a column's colors). */
  const baselineAll = (): void => {
    for (const v of viewers())
      baseline(v.root);
  };
  const writeProperties = async (el: Element, entries: [string, string][], capMs: number): Promise<number> => {
    const v = viewerOf(el);
    arm(v);
    // the state a "than before" claim compares with is the one the viewer had FINISHED, not
    // whatever it happened to be showing mid-render: a word cloud between two layouts reports no
    // words at all, and the claim after the write then has nothing to compare with
    await quiet(v).catch(() => undefined);
    baseline(el);
    const settled = settle(el, capMs);
    for (const [caption, text] of entries) {
      const p = property(v, caption);
      v.props[p.name] = convert(p, text);
    }
    return settled;
  };
  const canvasBox = (v: any): Box => {
    const r = anchorOf(v).getBoundingClientRect();
    return {x: r.x, y: r.y, width: r.width, height: r.height};
  };
  /** The named hit area in client coordinates, or the names the viewer reports instead. */
  const findArea = (el: Element, name: string, beforeChange = false): {box?: Box; has: string[]} => {
    const v = viewerOf(el);
    if (beforeChange)
      baseline(el);
    const areas: Record<string, Box> = v.getWidgetStatus()?.hitAreas ?? {};
    const has = Object.keys(areas);
    const key = has.find((k) => norm(k) === norm(name));
    if (!key)
      return {has};
    const r = areas[key];
    const cv = canvasBox(v);
    return {box: {x: cv.x + r.x, y: cv.y + r.y, width: r.width, height: r.height}, has};
  };
  const hitArea = (el: Element, name: string, beforeChange = false): Box => {
    const found = findArea(el, name, beforeChange);
    if (!found.box)
      throw new Error(`${viewerOf(el).type} has no "${name}" area right now; it has: ${found.has.join(', ') || 'none'}`);
    return found.box;
  };
  /** Every hit area the viewer reports, in client coordinates. */
  const areas = (el: Element): Record<string, Box> => {
    const v = viewerOf(el);
    const cv = canvasBox(v);
    const result: Record<string, Box> = {};
    for (const [key, r] of Object.entries(areasOf(v)))
      result[key] = {x: cv.x + r.x, y: cv.y + r.y, width: r.width, height: r.height};
    return result;
  };
  /** Painted pixels inside a hit area (the canvas may be scaled to the device). */
  const areaInk = (el: Element, name: string): number => {
    const v = viewerOf(el);
    const cv = canvasOf(v);
    return inkIn(pixelsOf(v), deviceRect(areasOf(v)[areaKey(v, name)], scaleOf(cv)));
  };
  /** A hit area's ink now against the snapshot's (the area must have been reported then too). */
  const areaChange = (el: Element, name: string): AreaChange => {
    const v = viewerOf(el);
    const before = snapshotOf(v);
    const inkBefore = areaInkBefore(before, norm(name));
    if (inkBefore === undefined)
      throw new Error(`${v.type} reported no "${name}" area at the snapshot; it had: ${Object.keys(before.rects).join(', ') || 'none'}`);
    const cv = canvasOf(v);
    return {ink: areaInk(el, name), inkBefore, rect: deviceRect(areasOf(v)[areaKey(v, name)], scaleOf(cv))};
  };
  /** The pixels that changed inside one hit area since the snapshot — the area's own repaint,
   * as opposed to the whole canvas'. */
  const areaDelta = (el: Element, name: string): number => {
    const v = viewerOf(el);
    const before = snapshotOf(v);
    const cv = canvasOf(v);
    const now = pixelsOf(v);
    if (before.bitmap.width !== now.width || before.bitmap.height !== now.height)
      throw new Error(`${v.type}: the canvas resized since the snapshot`);
    const r = deviceRect(areasOf(v)[areaKey(v, name)], scaleOf(cv));
    const a = before.bitmap.data;
    const b = now.data;
    let delta = 0;
    for (let y = r.y; y < Math.min(r.y + r.h, now.height); y++) {
      for (let x = r.x; x < Math.min(r.x + r.w, now.width); x++) {
        const i = (y * now.width + x) * 4;
        if (a[i] !== b[i] || a[i + 1] !== b[i + 1] || a[i + 2] !== b[i + 2] || a[i + 3] !== b[i + 3])
          delta++;
      }
    }
    return delta;
  };
  const change = (el: Element): CanvasChange => {
    const v = viewerOf(el);
    const cv = canvasOf(v);
    const img = pixelsOf(v);
    const now = histogram(img);
    const before = snapshotOf(v);
    const stats = statsOf(before);
    const delta = pixelDelta(before, img, now.colors);
    const view = areasOf(v)['view'];
    const scale = scaleOf(cv);
    const viewPx = view ? deviceRect(view, scale).w * deviceRect(view, scale).h : cv.width * cv.height;
    let selected = 0;
    try {
      selected = v.dataFrame.selection.trueCount;
    } catch { /* a viewer without a table */ }
    return {delta, ink: now.ink, inkBefore: stats.ink, hue: now.hue, hueBefore: stats.hue, selected, viewPx, dpr: window.devicePixelRatio};
  };
  const rangeChange = (el: Element): RangeChange => {
    const v = viewerOf(el);
    return {before: snapshots.get(v.root)?.range, now: rangeOf(v)};
  };
  const scaleChange = (el: Element): ScaleChange => {
    const v = viewerOf(el);
    return {before: snapshots.get(v.root)?.scale$, now: scaleRange(v)};
  };
  /** A range kept by the viewer's type, so it survives the viewer being closed and reopened (a
   * project round-trip). */
  const rememberRange = (el: Element): void => {
    const v = viewerOf(el);
    remembered[String(v.type)] = rangeOf(v);
  };
  const rememberedRange = (el: Element): RangeChange => {
    const v = viewerOf(el);
    return {before: remembered[String(v.type)], now: rangeOf(v)};
  };
  /** Whether the viewer painted since its snapshot, read once it is quiet: after the frame a
   * change would land on, and after whatever it says is still pending. */
  const stillness = async (el: Element): Promise<{renders: number; delta: number; where: string}> => {
    const v = viewerOf(el);
    const before = snapshotOf(v);
    await quiet(v);
    const cv = canvasOf(v);
    const delta = change(el).delta;
    return {renders: renders.get(v.root)!.count - before.renders, delta, where: delta ? changedWhere(v, cv, before.bitmap, pixelsOf(v)) : ''};
  };
  /** The value range against the snapshot's, read once the viewer is quiet (a reset that lands a
   * tick after the change is read, not missed). */
  const quietRangeChange = async (el: Element): Promise<RangeChange> => {
    await quiet(viewerOf(el));
    return rangeChange(el);
  };
  const quietValueChange = async (el: Element, name: string): Promise<Reading> => {
    await quiet(viewerOf(el));
    return valueChange(el, name);
  };
  /** One subscription per viewer and event, alive from "listens for" until `unlisten` (the
   * "should have fired" read, or the viewer closing). */
  const unlisten = (v: any, event?: string): void => {
    const all = listeners.get(v.root);
    if (!all)
      return;
    for (const e of event === undefined ? Object.keys(all) : [event]) {
      all[e]?.sub.unsubscribe();
      delete all[e];
    }
  };
  const listen = (el: Element, event: string): void => {
    const v = viewerOf(el);
    unlisten(v, event);
    const all = listeners.get(v.root) ?? {};
    listeners.set(v.root, all);
    const entry = {count: 0, sub: undefined as any};
    entry.sub = v.onEvent(event).subscribe(() => { entry.count++; });
    all[event] = entry;
  };
  const firedCount = (el: Element, event: string): number => listeners.get(viewerOf(el).root)?.[event]?.count ?? -1;
  const forget = (v: any): void => {
    unlisten(v);
    renders.get(v.root)?.sub?.unsubscribe();
    renders.delete(v.root);
  };
  /** Resolves once [read] says the same thing on [frames] consecutive frames after the first — what
   * it describes has stopped moving — or after [capMs]. */
  const stable = (read: () => string, capMs: number, frames: number): Promise<void> => new Promise((done) => {
    const t0 = Date.now();
    let last = '';
    let same = 0;
    const tick = (): void => {
      const now = read();
      same = now === last ? same + 1 : 0;
      last = now;
      if (same >= frames || Date.now() - t0 >= capMs)
        return done();
      requestAnimationFrame(tick);
    };
    requestAnimationFrame(tick);
  });
  /** The dock manager sizes the element it hosts, so a size written while it is still laying a
   * freshly docked viewer out is overwritten by the pass that follows. */
  const stableBox = (root: HTMLElement, capMs: number): Promise<void> => stable(() => {
    const r = root.getBoundingClientRect();
    return `${Math.round(r.width)}x${Math.round(r.height)}`;
  }, capMs, 2);
  /** A gesture aims at where the viewer has finished putting the thing, and a finished render is not
   * the end of that: a layout pass on the next frame can still move the whole viewer (a title just
   * set moved the pivot's grid down a row, and the right-click meant for a header landed on a cell).
   * The areas are relative to the anchor, so it is the anchor's box that must agree on two
   * consecutive frames. */
  const stableArea = (el: Element, capMs: number): Promise<void> => stable(() => {
    const r = anchorOf(viewerOf(el)).getBoundingClientRect();
    return [r.x, r.y, r.width, r.height].map(Math.round).join(',');
  }, capMs, 1);
  // the repaint a resize causes lands on the next task, so the settle is armed before the event
  const resize = async (el: Element, width: number | null, height: number | null, capMs: number): Promise<number> => {
    const root = viewerOf(el).root as HTMLElement;
    if (!sizes.has(root))
      sizes.set(root, {width: root.style.width, height: root.style.height});
    await stableBox(root, 1000);
    await quiet(viewerOf(el)).catch(() => undefined);
    baseline(el);
    const settled = settle(el, capMs);
    // and the dock keeps sizing it afterwards, whenever anything else in the view is laid out: the
    // size a step asked for is the fixture every claim after it is made against, so it is held
    // until the step that gives it back
    forced.get(root)?.disconnect();
    const want = {
      width: width === null ? null : `${width}px`,
      height: height === null ? null : `${height}px`,
    };
    const apply = (): void => {
      if (want.width !== null && root.style.width !== want.width)
        root.style.width = want.width;
      if (want.height !== null && root.style.height !== want.height)
        root.style.height = want.height;
    };
    const observer = new MutationObserver(apply);
    observer.observe(root, {attributes: true, attributeFilter: ['style']});
    forced.set(root, observer);
    apply();
    window.dispatchEvent(new Event('resize'));
    return settled;
  };
  const restoreSize = (el: Element, capMs: number): Promise<number> => {
    const root = viewerOf(el).root as HTMLElement;
    const size = sizes.get(root);
    if (!size)
      return Promise.resolve(0);
    forced.get(root)?.disconnect();
    forced.delete(root);
    baseline(el);
    const settled = settle(el, capMs);
    root.style.width = size.width;
    root.style.height = size.height;
    sizes.delete(root);
    window.dispatchEvent(new Event('resize'));
    return settled;
  };
  /** Subscribes to a `grok.events` stream before a gesture; `waitArmed` collects the outcome after. */
  const armEvent = (name: string, capMs: number): string => {
    const token = `t${++tokens}`;
    armed[token] = new Promise((resolve) => {
      let sub: any;
      try {
        sub = grok.events[name].subscribe((args: unknown) => {
          sub.unsubscribe();
          resolve(args ?? true);
        });
      } catch {
        resolve(undefined);
        return;
      }
      setTimeout(() => {
        try {
          sub?.unsubscribe();
        } catch { /* already gone */ }
        resolve(undefined);
      }, capMs);
    });
    return token;
  };
  const waitArmed = async (token: string): Promise<boolean> => {
    const result = await armed[token];
    delete armed[token];
    return result !== undefined;
  };
  /** Closes the open popup menu the platform's way (a click outside) and waits for it to be gone. */
  const closeMenu = async (): Promise<void> => {
    if (!document.querySelector('.d4-menu-popup'))
      return;
    const closed = armEvent('onContextMenuClosed', 1500);
    document.body.click();
    await waitArmed(closed);
  };
  /** Ready for a right-click: no popup open, `onContextMenuShown` armed (`onContextMenu` fires
   * before the menu exists). */
  const openMenu = async (capMs: number): Promise<string> => {
    await closeMenu();
    return armEvent('onContextMenuShown', capMs);
  };
  /** Where to right-click an element for its context menu — a named hit area, else the viewer's
   * `view` area, else the element's centre — with the canvas baseline taken and the menu armed. A
   * tree node, a card or a list row has no render to wait through and no areas: its own centre
   * (`settle` and `stableArea` throw for a non-viewer). */
  const menuPoint = async (el: Element, area: string | null, capMs: number): Promise<{x: number; y: number; token: string}> => {
    const v = findViewer(el);
    if (!v && area !== null)
      throw new Error(`"${area}" names a hit area, and the element is not a viewer of an open table view`);
    let box: Box | undefined;
    if (v) {
      await settle(el, 300).catch(() => undefined);
      await stableArea(el, 1000);
      try {
        box = hitArea(el, area ?? 'view', true);
      } catch (e) {
        if (area !== null)
          throw e;
      }
    }
    if (!box) {
      const r = el.getBoundingClientRect();
      box = {x: r.x, y: r.y, width: r.width, height: r.height};
    }
    return {x: box.x + box.width / 2, y: box.y + box.height / 2, token: await openMenu(capMs)};
  };
  const addViewer = (type: string): void => {
    const types: string[] = DG.Viewer.getViewerTypes();
    const exact = types.find((t) => norm(t) === norm(type));
    if (!exact)
      throw new Error(`no viewer type "${type}"; the platform has: ${types.join(', ')}`);
    arm(grok.shell.tv.addViewer(exact));
  };
  /** The balloons shown since the last read, and clears them. */
  const takeBalloons = (): Balloon[] => balloons.splice(0, balloons.length);
  const saveLayout = (): void => { layout = grok.shell.tv.saveLayout(); };
  /** Saves through the server and keeps only the id: "loads the saved layout" then fetches what
   * the server stored, so the round-trip covers the serialization too. */
  const saveLayoutToServer = async (): Promise<string> => {
    const l = grok.shell.tv.saveLayout();
    await grok.dapi.layouts.save(l);
    layout = {serverId: l.id};
    return l.id;
  };
  const loadLayout = async (): Promise<void> => {
    if (!layout)
      throw new Error('no layout saved in this feature');
    const saved = layout.serverId ? await grok.dapi.layouts.find(layout.serverId) : layout;
    if (!saved)
      throw new Error(`the server has no layout ${layout.serverId}`);
    grok.shell.tv.loadLayout(saved);
  };
  const deleteLayout = async (id: string): Promise<void> => {
    const saved = await grok.dapi.layouts.find(id);
    if (saved)
      await grok.dapi.layouts.delete(saved);
  };

  /** The function call a menu command starts, and the current table's columns before it: the
   * platform announces every call (`onBeforeRunAction` / `onAfterRunAction`), and the one whose
   * function is registered under the picked menu path is the command — its end is when the
   * command is done, whatever dialog it showed in between. */
  let command: {name: string; done: Promise<void>; started: number} | undefined;
  let columnsBefore: {table: any; names: string[]} | undefined;
  let commandArm: any;
  const armCommand = (path: string): void => {
    const t = grok.shell.t;
    columnsBefore = t ? {table: t.dart, names: t.columns.names()} : undefined;
    command = undefined;
    commandArm?.unsubscribe();
    const want = norm(path);
    const sub = grok.functions.onBeforeRunAction.subscribe((fc: any) => {
      const menu = fc?.func?.topMenu;
      if (!menu || norm(menu) !== want)
        return;
      sub.unsubscribe();
      if (commandArm === sub)
        commandArm = undefined;
      let resolve!: () => void;
      const done = new Promise<void>((r) => { resolve = r; });
      const after = grok.functions.onAfterRunAction.subscribe((ended: any) => {
        if (ended?.dart === fc.dart || ended?.id === fc.id) {
          after.unsubscribe();
          resolve();
        }
      });
      command = {name: String(fc.func?.nqName ?? fc.func?.name ?? path), done, started: Date.now()};
    });
    commandArm = sub;
    // a pick that started no call within 5 s is disarmed — this arm only, never a later pick's
    setTimeout(() => {
      sub.unsubscribe();
      if (commandArm === sub)
        commandArm = undefined;
    }, 5000);
  };
  const waitCommand = async (capMs: number): Promise<string> => {
    // the click's handler may start the call a task or two later
    const t0 = Date.now();
    while (!command && commandArm && Date.now() - t0 < 5000)
      await new Promise((r) => setTimeout(r, 10));
    if (!command)
      throw new Error('no menu command has started a function call in this scenario (a Dart command, or the menu item ran nothing)');
    const c = command;
    const timeout = new Promise<'timeout'>((r) => setTimeout(() => r('timeout'), capMs));
    if (await Promise.race([c.done.then(() => 'done'), timeout]) === 'timeout')
      throw new Error(`${c.name} has been running for ${Math.round((Date.now() - c.started) / 1000)} s`);
    return c.name;
  };
  /** The current table's columns now and before the last menu command; `same` says whether the
   * table is still the one the command started on. */
  const columnsSince = (): {before: string[] | null; now: string[]; same: boolean} => {
    const t = grok.shell.t;
    return {before: columnsBefore?.names ?? null, now: t ? t.columns.names() : [], same: !!t && columnsBefore?.table === t.dart};
  };

  /** Custom platform events (`grok.events.fireCustomEvent`) by id, counted from "listens for"
   * until the page resets; a read takes the count and the last arguments and zeroes them. */
  const customEvents = new Map<string, {count: number; last: unknown; sub: any}>();
  const listenCustom = (id: string): void => {
    if (customEvents.has(id))
      return;
    const entry = {count: 0, last: undefined as unknown, sub: undefined as any};
    entry.sub = grok.events.onCustomEvent(id).subscribe((args: unknown) => { entry.count++; entry.last = args; });
    customEvents.set(id, entry);
  };
  const customFired = (id: string, take: boolean): {count: number; last: unknown} => {
    const entry = customEvents.get(id);
    if (!entry)
      return {count: -1, last: undefined};
    let last: unknown;
    try {
      last = JSON.parse(JSON.stringify(entry.last ?? null));
    }
    catch {
      last = String(entry.last);
    }
    const read = {count: entry.count, last};
    if (take) {
      entry.count = 0;
      entry.last = undefined;
    }
    return read;
  };

  w.__bdd = {table, tableNamed, col, rowFacts, setRows,
    viewerOf, arm, stampAll, settle, quiet, readProperty, writeProperties, findArea, hitArea, areas, areaInk, areaChange, areaDelta, areaColors,
    areaRectChange, legendState: (el: Element) => legendState(viewerOf(el)), legendChange, rememberValue, rememberedValue,
    snapshot, baselineAll, change, rangeChange, quietRangeChange, scaleChange, valueChange, quietValueChange, rememberRange, rememberedRange, stillness,
    palette, tableOf, listen, unlisten, firedCount, resize, restoreSize, armEvent, waitArmed, closeMenu, openMenu, menuPoint, stableArea, addViewer,
    takeBalloons, saveLayout, saveLayoutToServer, loadLayout, deleteLayout, ink, armCommand, waitCommand, columnsSince, listenCustom, customFired};
  stampAll();
  grok.events.onViewerAdded.subscribe((a: any) => arm(a?.args?.viewer));
  grok.events.onViewerClosed.subscribe((a: any) => a?.args?.viewer && forget(a.args.viewer));
  grok.events.onEvent('d4-balloon-shown').subscribe((a: any) => balloons.push({type: String(a?.args?.type ?? ''), message: String(a?.args?.message ?? '')}));
}

const installed = new WeakSet<Page>();
const watched = new WeakSet<Page>();

/** Installs the in-page side once the shell is up. Cheap on every step: the page is remembered,
 * and forgotten when its main frame navigates (a reload drops the in-page side). */
export async function installViewerRuntime(page: Page): Promise<void> {
  if (installed.has(page))
    return;
  await page.evaluate(install);
  installed.add(page);
  if (!watched.has(page)) {
    watched.add(page);
    page.on('framenavigated', (frame) => {
      if (frame === page.mainFrame())
        installed.delete(page);
    });
  }
}

/** `page.evaluate` with the in-page runtime installed — for a callback that uses `__bdd`. */
export async function evaluate<R>(page: Page, fn: (arg: any) => R | Promise<R>, arg?: unknown): Promise<R> {
  await installViewerRuntime(page);
  return page.evaluate(fn as (arg: unknown) => R, arg);
}
