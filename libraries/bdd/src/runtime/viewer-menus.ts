/* The Dart context menu of a viewer, by path: a real right-click at a point the viewer reports,
   the popup awaited through `onContextMenuShown` (it fires once the menu is in the DOM), groups
   opened by pointer moves, the leaf clicked. Items match their own label (a group item contains
   its children's labels too), and a group's children live in the same popup, hidden until it
   opens — so "is it shown" is whether any match is visible. */
import {Locator, Page} from '@playwright/test';
import {expect, pollMs} from './patience.js';
import type {ElementRef} from './args.js';
import {exactText} from './locate.js';
import {evaluate} from './viewer-runtime.js';
import {onViewer} from './viewers.js';

const POPUP = '.d4-menu-popup';
const SEP = /\s*[>|]\s*/;
const MENU_ITEM = 'xpath=ancestor::*[contains(concat(" ", normalize-space(@class), " "), " d4-menu-item ")][1]';

/** Arms a `grok.events.<name>` subscription before a gesture; the returned function waits for it. */
export async function armEvent(page: Page, name: string, capMs = 3000): Promise<() => Promise<boolean>> {
  const token: string = await evaluate(page, ([n, cap]) => (window as any).__bdd.armEvent(n, cap), [name, capMs] as [string, number]);
  return () => page.evaluate((t) => (window as any).__bdd.waitArmed(t), token);
}

export async function closeContextMenu(page: Page): Promise<void> {
  await evaluate(page, () => (window as any).__bdd.closeMenu(), undefined);
  await expect(page.locator(POPUP)).toHaveCount(0);
}

async function rightClickArmed(page: Page, x: number, y: number, token: string): Promise<Locator> {
  await page.mouse.click(x, y, {button: 'right'});
  if (!await page.evaluate((t) => (window as any).__bdd.waitArmed(t), token))
    throw new Error(`no context menu opened at (${Math.round(x)}, ${Math.round(y)})`);
  return page.locator(POPUP).last();
}

export async function openContextMenuAt(page: Page, x: number, y: number): Promise<Locator> {
  const token: string = await evaluate(page, (cap) => (window as any).__bdd.openMenu(cap), 3000);
  return rightClickArmed(page, x, y, token);
}

/** The context menu of an element: at a named hit area, else its `view` area when it reports
 * one, else its centre. Three roundtrips: the point and the arming in one, the click, the wait. */
export async function openContextMenuOf(page: Page, target: ElementRef, area?: string): Promise<Locator> {
  const {x, y, token} = await onViewer(page, target, (el, [a, cap]) => (window as any).__bdd.menuPoint(el, a, cap), [area ?? null, 3000] as [string | null, number]);
  return rightClickArmed(page, x, y, token);
}

/** A group item holds its children under its own label: the labels inside it are the children's. */
const CHILDREN = '.d4-menu-item-container';

/** Every label of the open menus reading exactly that — a viewer menu can hold two groups of the
 * same name (the axis "Annotations" and the viewer's own) — or, inside a group, of its children. */
export function menuLabelsMatching(page: Page, label: string, within?: Locator): Locator {
  return (within ? within.locator(CHILDREN) : page.locator(POPUP)).locator('.d4-menu-item-label', {hasText: exactText(label)});
}

export function menuItems(page: Page, label: string, within?: Locator): Locator {
  return menuLabelsMatching(page, label, within).locator(MENU_ITEM);
}

async function menuShows(page: Page, label: string, within?: Locator): Promise<boolean> {
  return (await menuLabelsMatching(page, label, within).filter({visible: true}).count()) > 0;
}

export function menuItem(page: Page, label: string): Locator {
  return menuItems(page, label).first();
}

/** What the open menus show, by item: its label, or the caption of an element item (a property
 * editor the menu embeds — the PC plot's column picker under "Columns" — carries it as `data-source`). */
async function visibleMenuLabels(page: Page, within?: Locator): Promise<string> {
  const items = (within ? within.locator(CHILDREN) : page.locator(POPUP)).locator('.d4-menu-item').filter({visible: true});
  const labels = await items.evaluateAll((els) => els.map((e) => e.querySelector(':scope > .d4-menu-item-label')?.textContent ?? e.getAttribute('data-source') ?? ''));
  return labels.map((t) => t.trim()).filter(Boolean).join(' | ');
}

/** A group item opens on a pointer move over it: the move enters from the left, since a pointer
 * already resting on the item (a hover before the menu was reopened) would move nowhere. When
 * several items share the label, each is tried until [wanted] shows up among its children — so
 * "HEIGHT > Chart type" is the series' group, not the chart-wide "Chart Type" beside it. Returns
 * the group item, the scope of the next segment. */
async function openGroup(page: Page, label: string, wanted: string | undefined, within?: Locator): Promise<Locator> {
  const candidates = menuItems(page, label, within);
  await candidates.filter({visible: true}).first().waitFor({state: 'visible', timeout: 5000});
  const count = await candidates.count();
  // a group already open, or one the menu renders inline, needs no hover
  for (let i = 0; wanted !== undefined && i < count; i++) {
    if (await menuShows(page, wanted, candidates.nth(i)))
      return candidates.nth(i);
  }
  for (let i = 0; i < count; i++) {
    const item = candidates.nth(i);
    if (!(await item.isVisible()))
      continue;
    const box = await item.boundingBox();
    if (!box)
      continue;
    const cy = box.y + box.height / 2;
    await page.mouse.move(Math.max(0, box.x - 8), cy);
    await page.mouse.move(box.x + box.width * 0.75, cy);
    // a move back to the left is what the menu takes as leaving a sibling's open submenu (a move
    // to the right stays inside the triangle that keeps it open), and the same move opens this one
    await page.mouse.move(box.x + box.width / 2, cy);
    if (wanted === undefined)
      return item;
    const opened = await expect.poll(() => menuShows(page, wanted, item), {timeout: pollMs(2000)}).toBe(true).then(() => true, () => false);
    if (opened) {
      // step into the flyout along the group's own row: the menu hides a submenu when the pointer
      // leaves the group item at a steep angle, which would take the item away mid-click
      const target = await menuLabelsMatching(page, wanted, item).filter({visible: true}).first().boundingBox();
      if (target != null) {
        await page.mouse.move(target.x + 4, cy);
        await page.mouse.move(target.x + 4, target.y + target.height / 2);
      }
      return item;
    }
  }
  throw new Error(`the "${label}" group did not show "${wanted}"; the menu shows: ${await visibleMenuLabels(page, within)}`);
}

/** `Misc > Show Inside Values`: hovers the groups, each found inside the one before, clicks the
 * leaf. The item ancestor is the usual click target, but an inline group lays the box out on the
 * label itself — decided before the click, and clicked once: a click whose handler rebuilds the
 * viewer synchronously can outlast a short cap with its work done, and a second click would undo
 * a toggle. */
export async function pickMenuPath(page: Page, path: string): Promise<void> {
  const segments = path.split(SEP).filter((s) => s.length > 0);
  let group: Locator | undefined;
  for (let i = 0; i < segments.length; i++) {
    const last = i === segments.length - 1;
    if (!last) {
      group = await openGroup(page, segments[i], segments[i + 1], group).catch(async () => {
        throw new Error(`no "${segments[i]}" in the menu; it shows: ${await visibleMenuLabels(page, group)}`);
      });
      continue;
    }
    const label = menuLabelsMatching(page, segments[i], group).filter({visible: true}).first();
    await label.waitFor({state: 'visible', timeout: 3000}).catch(async () => {
      throw new Error(`no "${segments[i]}" in the menu; it shows: ${await visibleMenuLabels(page, group)}`);
    });
    const item = label.locator(MENU_ITEM);
    await (await item.count() > 0 ? item : label).click();
  }
}

/** The labels under [path], opening every group of it first ("" reads the top level). */
export async function menuLabels(page: Page, path: string): Promise<string[]> {
  const segments = path.split(SEP).filter((s) => s.length > 0);
  let group: Locator | undefined;
  for (let i = 0; i < segments.length; i++)
    group = await openGroup(page, segments[i], segments[i + 1], group);
  return (await visibleMenuLabels(page, group)).split(' | ').filter(Boolean);
}
