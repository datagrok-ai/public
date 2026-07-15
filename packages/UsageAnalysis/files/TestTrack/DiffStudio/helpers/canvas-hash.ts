import { Page, Locator } from '@playwright/test';
import { createHash } from 'crypto';

/**
 * Hash the contents of a chart so two snapshots can be compared without pixel-by-pixel diffing.
 *
 * D4 viewers typically render to **two** canvases: a `canvas` for static data (the chart line)
 * and an `overlay` for interactive highlights. Selecting `.d4-viewer canvas` non-deterministically
 * picks one — often the overlay, which is blank and identical across redraws. To avoid that
 * pitfall, this helper concatenates the dataURLs of **every** canvas under the selector and
 * hashes the join — any one of them changing produces a different hash.
 *
 * Returns a short hex digest, or `'<missing>'` if no canvas was found.
 *
 * Usage:
 *   const h1 = await canvasHash(page, '.d4-viewer');
 *   await page.locator(inputEditor('switch-at')).fill('150');
 *   await page.waitForTimeout(2000);
 *   const h2 = await canvasHash(page, '.d4-viewer');
 *   expect(h2).not.toBe(h1);  // chart redrew
 */
export async function canvasHash(page: Page, scopeSelector: string): Promise<string> {
  const joined = await page.evaluate((sel) => {
    const scope = document.querySelector(sel);
    if (!scope) return null;
    const canvases = scope.tagName === 'CANVAS'
      ? [scope as HTMLCanvasElement]
      : Array.from(scope.querySelectorAll('canvas')) as HTMLCanvasElement[];
    if (canvases.length === 0) return null;
    const parts: string[] = [];
    for (const c of canvases) {
      try { parts.push(c.toDataURL('image/png')); }
      catch { /* tainted / cross-origin — skip */ }
    }
    return parts.join('|');
  }, scopeSelector);
  if (!joined) return '<missing>';
  return createHash('sha1').update(joined).digest('hex').slice(0, 16);
}

/** Hash a single canvas located via Playwright Locator. */
export async function canvasHashLocator(loc: Locator): Promise<string> {
  const dataUrl = await loc.evaluate((c: Element) => {
    const canvas = c as HTMLCanvasElement;
    try { return canvas.toDataURL('image/png'); }
    catch { return null; }
  });
  if (!dataUrl) return '<missing>';
  return createHash('sha1').update(dataUrl).digest('hex').slice(0, 16);
}

/**
 * Hash the PNG pixel data of the bounding box of an element via Playwright's `element.screenshot()`.
 * This works for canvas-based charts inside Vue components where the canvas itself is not
 * directly addressable — Playwright captures whatever the browser is actually painting on
 * screen, which is exactly what the user sees.
 *
 * Selector should resolve to a single visible element. Falls back to `<missing>` if not found
 * or not visible.
 */
export async function elementScreenshotHash(page: Page, selector: string): Promise<string> {
  const loc = page.locator(selector).first();
  if (await loc.count() === 0) return '<missing>';
  try {
    const buf = await loc.screenshot({ timeout: 5_000 });
    return createHash('sha1').update(buf).digest('hex').slice(0, 16);
  } catch {
    return '<missing>';
  }
}

/**
 * Fallback signal for "the chart redrew" when the canvas hash is unreliable — e.g. inside the
 * Compute2 RichFunctionView used by saved models in Model Hub. Sums every finite numeric value
 * across all dataframes the platform currently knows about. If the solver re-ran on input
 * change, the sum changes.
 *
 * Searches `grok.shell.t`, `grok.shell.v.dataFrame`, and every entry in `grok.shell.tables` —
 * Compute2 RFV usually exposes its result table via one of the latter two.
 *
 * Uses page.evaluate to reach grok.shell because the dataframe is part of platform state, not
 * the DOM — this is the same kind of read-only API access the user permitted as a fallback.
 */
export async function dataframeFingerprint(page: Page): Promise<number | null> {
  return await page.evaluate(() => {
    const win = window as any;
    const candidates: any[] = [];
    const t = win.grok?.shell?.t;
    if (t) candidates.push(t);
    const v = win.grok?.shell?.v;
    if (v?.dataFrame) candidates.push(v.dataFrame);
    const all = win.grok?.shell?.tables;
    if (all && typeof all[Symbol.iterator] === 'function') {
      for (const tbl of all) if (tbl) candidates.push(tbl);
    }
    if (candidates.length === 0) return null;

    let s = 0;
    let counted = 0;
    const seenKeys = new Set<string>();
    for (const df of candidates) {
      // Dedupe by name+rowCount+columnCount signature
      const key = `${df.name ?? '?'}|${df.rowCount}|${df.columns?.length ?? 0}`;
      if (seenKeys.has(key)) continue;
      seenKeys.add(key);
      const cols = df.columns?.toList?.() ?? df.columns;
      if (!cols) continue;
      for (const col of cols) {
        const type: string = col.type;
        if (type !== 'double' && type !== 'float' && type !== 'int' && type !== 'qnum') continue;
        const len = col.length ?? 0;
        for (let i = 0; i < len; i++) {
          const x = col.get(i);
          if (typeof x === 'number' && isFinite(x)) { s += x; counted++; }
        }
      }
    }
    return counted > 0 ? s : null;
  });
}
