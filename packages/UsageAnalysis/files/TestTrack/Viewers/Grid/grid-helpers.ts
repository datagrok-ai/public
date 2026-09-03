import {Page} from '@playwright/test';
import * as v from '../../helpers/viewers';

declare const grok: any;

export const DEMOG = 'System:DemoFiles/demog.csv';
export const OVERLAY = '[name="viewer-Grid"] canvas[name="overlay"]';

export interface Point { x: number; y: number; }

export function cellCenter(page: Page, col: string, row: number): Promise<Point> {
  return page.evaluate(({c, r}) => {
    const db = grok.shell.tv.grid.cell(c, r).documentBounds;
    return {x: db.x + db.width / 2, y: db.y + db.height / 2};
  }, {c: col, r: row});
}

export function headerCenter(page: Page, col: string): Promise<Point> {
  return page.evaluate((c) => {
    const grid = grok.shell.tv.grid;
    const db = grid.cell(c, 0).documentBounds;
    return {x: db.x + db.width / 2, y: db.y - grid.colHeaderHeight / 2};
  }, col);
}

export function focusGrid(page: Page): Promise<void> {
  return page.evaluate((sel) => (document.querySelector(sel) as HTMLElement).focus(), OVERLAY);
}

export async function openGridMenu(page: Page, at: Point): Promise<void> {
  await page.evaluate(({x, y, sel}) => {
    const overlay = document.querySelector(sel) as HTMLElement;
    const cm = {bubbles: true, cancelable: true, clientX: x, clientY: y, button: 2, buttons: 2} as any;
    overlay.dispatchEvent(new MouseEvent('mousedown', cm));
    overlay.dispatchEvent(new MouseEvent('mouseup', cm));
    overlay.dispatchEvent(new MouseEvent('contextmenu', cm));
  }, {x: at.x, y: at.y, sel: OVERLAY});
  await page.locator('.d4-menu-popup .d4-menu-item').first().waitFor({timeout: 5000});
}

export async function closeGridMenu(page: Page): Promise<void> {
  await page.evaluate(() => {
    document.body.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, clientX: 5, clientY: 5}));
  });
  await v.pollValue(() => page.evaluate(() => Array.from(document.querySelectorAll('.d4-menu-popup'))
    .every((p) => (p as HTMLElement).offsetParent === null)), (gone) => gone, 1000, 50);
  await page.evaluate(() => { for (const p of Array.from(document.querySelectorAll('.d4-menu-popup'))) p.remove(); });
}

/**
 * Opens the grid context menu at `at`, walks the `groups` chain and clicks `leaf`.
 * Submenu containers are shown directly instead of hover-revealed, so the walk has no timing;
 * the caller waits for whatever the leaf sets in motion.
 */
export async function clickMenuLeaf(page: Page, at: Point, groups: string[], leaf: string): Promise<boolean> {
  await openGridMenu(page, at);
  const clicked = await page.evaluate(async ({groups, leaf}) => {
    const w = window as any;
    const popup = () => Array.from(document.querySelectorAll('.d4-menu-popup')).pop();
    for (const g of groups) {
      const group = await w.__poll(() => popup()?.querySelector(`[name="${g}"]`) ?? null, (e: any) => !!e, 3000, 50);
      if (!group) return false;
      const b = group.getBoundingClientRect();
      for (const t of ['mouseover', 'mousemove'])
        group.dispatchEvent(new MouseEvent(t, {bubbles: true, clientX: b.x + 5, clientY: b.y + 5}));
      const container = group.querySelector('.d4-menu-item-container.d4-vert-menu') as HTMLElement | null;
      if (container) container.style.display = 'flex';
    }
    const item = await w.__poll(() => {
      const el = popup()?.querySelector(`[name="${leaf}"]`) as HTMLElement | null;
      return el && el.getBoundingClientRect().width > 0 ? el : null;
    }, (e: any) => !!e, 3000, 50);
    if (!item) {
      console.error(`[clickMenuLeaf] no leaf ${leaf}; items: ` + Array.from(popup()?.querySelectorAll('[name^="div-"]') ?? [])
        .map((e) => e.getAttribute('name')).filter((n) => n && n.startsWith(groups[groups.length - 1] ?? 'div-')).join(', '));
      return false;
    }
    const r = item.getBoundingClientRect();
    const o = {bubbles: true, cancelable: true, clientX: r.x + r.width / 2, clientY: r.y + r.height / 2, button: 0} as any;
    for (const t of ['mousedown', 'mouseup', 'click']) item.dispatchEvent(new MouseEvent(t, o));
    return true;
  }, {groups, leaf});
  await closeGridMenu(page);
  return clicked;
}

export function pinViaMenu(page: Page, at: Point, leaf: string): Promise<boolean> {
  return clickMenuLeaf(page, at, ['div-Pin'], leaf);
}

/** Opens the grid's own settings (the gear at the grid corner) in the property panel. */
export async function openGridSettings(page: Page, probe = 'prop-row-height'): Promise<boolean> {
  const rows = page.locator(`[name="${probe}"]`);
  if (await rows.count() > 0) return true;
  const gearAt = () => page.evaluate(() => {
    const gear = document.querySelector('.d4-grid-settings-icon') as HTMLElement | null;
    if (!gear) return null;
    const r = gear.getBoundingClientRect();
    return r.width > 0 ? {x: r.x + r.width / 2, y: r.y + r.height / 2} : null;
  });
  const built = () => rows.first().waitFor({state: 'attached', timeout: 4000}).then(() => true).catch(() => false);
  const box = await gearAt();
  if (box) {
    await page.mouse.move(box.x - 40, box.y + 20);
    await page.mouse.move(box.x, box.y, {steps: 6});
    await page.mouse.click(box.x, box.y);
    if (await built()) return true;
  }
  await page.evaluate(() => {
    const gear = document.querySelector('.d4-grid-settings-icon') as HTMLElement | null;
    if (!gear) return;
    const r = gear.getBoundingClientRect();
    const o = {bubbles: true, cancelable: true, clientX: r.x + r.width / 2, clientY: r.y + r.height / 2, button: 0} as any;
    for (const t of ['mouseover', 'mousedown', 'mouseup', 'click']) gear.dispatchEvent(new MouseEvent(t, o));
  });
  return built();
}

/**
 * Picks `name` in the column combobox found by `comboSelector` (the element holding the
 * .d4-column-selector-column label): a trusted click opens the popup, the typed text is
 * confirmed in the popup's search input before Enter, and the label is read back.
 */
export async function pickColumnInCombo(page: Page, comboSelector: string, name: string): Promise<boolean> {
  const labelOf = (sel: string) => {
    const combo = document.querySelector(sel) as HTMLElement | null;
    const label = combo?.querySelector('.d4-column-selector-column') as HTMLElement | null;
    return (label?.textContent ?? '').trim();
  };
  for (let attempt = 0; attempt < 2; attempt++) {
    const at = await page.evaluate((sel) => {
      const host = document.querySelector(sel) as HTMLElement | null;
      const combo = (host?.matches('[name^="div-column-combobox"]') ? host : host?.querySelector('[name^="div-column-combobox"]') ?? host) as HTMLElement | null;
      if (!combo) return null;
      const label = combo.querySelector('.d4-column-selector-column') as HTMLElement | null;
      const lr = label?.getBoundingClientRect();
      const r = lr && lr.width > 0 && lr.height > 0 ? lr : combo.getBoundingClientRect();
      return {x: r.x + r.width / 2, y: r.y + r.height / 2};
    }, comboSelector);
    if (!at) return false;
    await page.mouse.click(at.x, at.y);
    const opened = await page.waitForFunction(() => !!document.querySelector('.d4-column-selector-backdrop'),
      null, {timeout: 2000}).then(() => true).catch(() => false);
    if (!opened) continue;
    await page.keyboard.type(name.toLowerCase());
    const typed = await v.pollValue(() => page.evaluate(() => {
      const el = document.activeElement as HTMLInputElement | null;
      return el?.classList.contains('d4-column-selector-search-input') ? el.value : null;
    }), (t) => t === name.toLowerCase(), 1000, 50);
    if (typed !== name.toLowerCase()) {
      await page.keyboard.press('Escape').catch(() => {});
      continue;
    }
    await page.keyboard.press('Enter');
    const picked = await v.pollValue(() => page.evaluate(labelOf, comboSelector), (t) => t === name, 1500, 50);
    if (picked === name) return true;
    await page.keyboard.press('Escape').catch(() => {});
  }
  return false;
}

export interface ErrorTracker { list: string[]; count(): number; stop(): void; }

/** Collects console errors (and page errors) until stop(); every spec on the shared page must stop it. */
export function trackErrors(page: Page, isNoise: (t: string) => boolean = () => false): ErrorTracker {
  const list: string[] = [];
  const onConsole = (m: any) => { if (m.type() === 'error' && !isNoise(m.text())) list.push(m.text()); };
  const onPageError = (e: any) => { if (!isNoise(String(e))) list.push(String(e)); };
  page.on('console', onConsole);
  page.on('pageerror', onPageError);
  return {
    list,
    count: () => list.length,
    stop: () => { page.off('console', onConsole); page.off('pageerror', onPageError); },
  };
}

export const BENIGN_NOISE = (t: string): boolean =>
  /Unable to find element in cloned iframe/i.test(t) ||
  /NullError: method not found: '[a-zA-Z]+' on null/i.test(t) ||
  /Stack trace [A-Za-z0-9]+/.test(t);

export interface ShellFlags { showContextPanel: boolean; showToolbox: boolean; showProperties: boolean; showHelp: boolean; }

export function readShellFlags(page: Page): Promise<ShellFlags> {
  return page.evaluate(() => {
    const w = grok.shell.windows;
    return {showContextPanel: w.showContextPanel, showToolbox: w.showToolbox, showProperties: w.showProperties, showHelp: w.showHelp};
  });
}

/** Closes every view and puts the shell window flags back: what a spec on the shared page owes the next one. */
export async function leaveShellClean(page: Page, flags: ShellFlags): Promise<void> {
  await v.closeAllAndWait(page).catch(() => {});
  await page.evaluate((f) => {
    const w = grok.shell.windows;
    for (const k of Object.keys(f)) { try { w[k] = (f as any)[k]; } catch (_) {} }
    grok.shell.o = null;
  }, flags).catch(() => {});
}
