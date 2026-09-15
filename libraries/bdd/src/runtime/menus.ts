/* The Dart top menu (the ribbon's menu bar: File, Edit, ..., a package's group) driven by path
   the way the context menu is by `pickMenuPath`. The items carry the platform's names —
   `div-Bio---Analyze---Sequence-Space...` for "Bio > Analyze > Sequence Space..." — so every
   segment is found by name, never by scanning labels. A top-level group opens on the pointer
   entering it (a click toggles it — a second click on an open one closes it); a vertical group
   opens on a pointer move inside it, and a pointer that already rests on it moves nowhere, so
   the entry is two moves within the item. A menu bar too narrow for its groups folds the rest
   into a "more" group, where the folded items are vertical. */
import {Locator, Page} from '@playwright/test';
import {expect} from './patience.js';
import {installViewerRuntime, baselineAll} from './viewers.js';

const MORE = '[role="menubar"] .d4-menu-item-more';

/** `Bio > Analyze > Sequence Space...` → the segments and the platform names of the path so far. */
export function menuNames(path: string): {segments: string[]; names: string[]} {
  const segments = path.split(/\s*[>|]\s*/).filter((s) => s.length > 0);
  const names = segments.map((_, i) => 'div-' + segments.slice(0, i + 1).map((s) => s.trim().replace(/\s+/g, '-')).join('---'));
  return {segments, names};
}

/** A pointer move within the item opens a vertical group. */
async function enterGroup(item: Locator): Promise<void> {
  const box = await item.boundingBox();
  if (!box)
    throw new Error('the menu item has no box');
  const cx = box.x + box.width / 2;
  const cy = box.y + box.height / 2;
  await item.page().mouse.move(cx + 1, cy);
  await item.page().mouse.move(cx, cy);
}

export async function visibleLabels(page: Page, groupName: string): Promise<string> {
  const labels = await page.locator(`[name="${groupName}"] > .d4-menu-item-container .d4-menu-item-label`).filter({visible: true}).allTextContents();
  return labels.map((s) => s.trim()).filter(Boolean).join(' | ');
}

/** Walks the path: the top-level group by hover (through the "more" group when the bar folded
 * it), every group below by a move inside it, and the leaf by a click when `pick` is set,
 * else left open with its group. */
export async function openTopMenu(page: Page, path: string, pick: boolean): Promise<void> {
  const {segments, names} = menuNames(path);
  if (segments.length === 0)
    throw new Error('an empty menu path');
  // a package's group shows once the table has a column it applies to, a moment after the
  // detection the open step waited for: the bar rebuilds on its own schedule
  const top = page.locator(`.d4-menu-item-horz[name="${names[0]}"]`).filter({visible: true});
  await top.first().waitFor({state: 'visible', timeout: 5000}).catch(() => undefined);
  if (await top.count() > 0) {
    await top.first().hover();
  }
  else {
    const more = page.locator(MORE).filter({visible: true});
    if (await more.count() === 0) {
      const bar = await page.locator('[role="menubar"] .d4-menu-item-horz > .d4-menu-item-label').filter({visible: true}).allTextContents();
      throw new Error(`no "${segments[0]}" in the top menu; it shows: ${bar.map((s) => s.trim()).filter(Boolean).join(' | ') || 'nothing'}`);
    }
    await more.first().hover();
    const folded = page.locator(`[name="${names[0]}"]`).filter({visible: true}).first();
    await folded.waitFor({state: 'visible', timeout: 5000}).catch(() => {
      throw new Error(`no "${segments[0]}" in the top menu's overflow group`);
    });
    await enterGroup(folded);
  }
  for (let i = 1; i < segments.length; i++) {
    const item = page.locator(`[name="${names[i]}"]`).first();
    await item.waitFor({state: 'visible', timeout: 5000}).catch(async () => {
      throw new Error(`no "${segments[i]}" in the ${segments.slice(0, i).join(' > ')} menu; it shows: ${await visibleLabels(page, names[i - 1]) || 'nothing'}`);
    });
    if (i < segments.length - 1 || !pick)
      await enterGroup(item);
    else
      await item.click();
  }
}

/** Picks a command: the viewers' baseline and the table's columns are taken first, and the
 * function call the command starts is watched (`waitCommand`). */
export async function pickTopMenu(page: Page, path: string): Promise<void> {
  await installViewerRuntime(page);
  await baselineAll(page);
  const {segments, names} = menuNames(path);
  if (segments.length < 2)
    throw new Error(`"${path}" names a group, not a command: a command is "Group > Item"`);
  await openTopMenu(page, segments.slice(0, -1).join(' > '), false);
  const leaf = page.locator(`[name="${names[names.length - 1]}"]`).first();
  await leaf.waitFor({state: 'visible', timeout: 5000}).catch(async () => {
    throw new Error(`no "${segments[segments.length - 1]}" in the ${segments.slice(0, -1).join(' > ')} menu; it shows: ${await visibleLabels(page, names[names.length - 2]) || 'nothing'}`);
  });
  await page.evaluate((p) => (window as any).__bdd.armCommand(p), segments.join(' | '));
  await leaf.click();
}

/** The bar's group closes when the pointer leaves it (Escape is not a key it listens to). */
export async function closeTopMenu(page: Page): Promise<void> {
  const size = page.viewportSize() ?? {width: 1280, height: 800};
  await page.mouse.move(size.width / 2, size.height / 2);
  await expect(page.locator('[role="menubar"] .d4-vert-menu').filter({visible: true})).toHaveCount(0);
}

/** Resolves once the last menu command's function call has ended (the platform's
 * `onAfterRunAction` for it); the name of the function. */
export async function waitCommand(page: Page, capMs = 120000): Promise<string> {
  await installViewerRuntime(page);
  return page.evaluate((cap) => (window as any).__bdd.waitCommand(cap), capMs);
}

export interface ColumnsSince {
  before: string[] | null;
  now: string[];
  same: boolean;
}

export async function columnsSince(page: Page): Promise<ColumnsSince> {
  await installViewerRuntime(page);
  return page.evaluate(() => (window as any).__bdd.columnsSince());
}
