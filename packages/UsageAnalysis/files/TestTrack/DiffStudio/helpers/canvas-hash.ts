import { Page, Locator } from '@playwright/test';
import { createHash } from 'crypto';

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
      catch {  }
    }
    return parts.join('|');
  }, scopeSelector);
  if (!joined) return '<missing>';
  return createHash('sha1').update(joined).digest('hex').slice(0, 16);
}

export async function canvasHashLocator(loc: Locator): Promise<string> {
  const dataUrl = await loc.evaluate((c: Element) => {
    const canvas = c as HTMLCanvasElement;
    try { return canvas.toDataURL('image/png'); }
    catch { return null; }
  });
  if (!dataUrl) return '<missing>';
  return createHash('sha1').update(dataUrl).digest('hex').slice(0, 16);
}

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
