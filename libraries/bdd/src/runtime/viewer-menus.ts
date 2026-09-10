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

/** Every label of the open menus reading exactly that — a viewer menu can hold two groups of the
 * same name (the axis "Annotations" and the viewer's own). */
export function menuLabelsMatching(page: Page, label: string): Locator {
  return page.locator(POPUP).locator('.d4-menu-item-label', {hasText: exactText(label)});
}

export function menuItems(page: Page, label: string): Locator {
  return menuLabelsMatching(page, label).locator(MENU_ITEM);
}

async function menuShows(page: Page, label: string): Promise<boolean> {
  return (await menuLabelsMatching(page, label).filter({visible: true}).count()) > 0;
}

export function menuItem(page: Page, label: string): Locator {
  return menuItems(page, label).first();
}

async function visibleMenuLabels(page: Page): Promise<string> {
  const labels = await page.locator(POPUP).locator('.d4-menu-item-label').filter({visible: true}).allTextContents();
  return labels.map((t) => t.trim()).filter(Boolean).join(' | ');
}

/** A group item opens on a pointer move over it: the move enters from the left, since a pointer
 * already resting on the item (a hover before the menu was reopened) would move nowhere. When
 * several items share the label, each is tried until [wanted] shows up. */
async function openGroup(page: Page, label: string, wanted?: string): Promise<void> {
  // a group already open, or one the menu renders inline, needs no hover
  if (wanted !== undefined && await menuShows(page, wanted))
    return;
  const candidates = menuItems(page, label);
  await candidates.filter({visible: true}).first().waitFor({state: 'visible', timeout: 5000});
  const count = await candidates.count();
  for (let i = 0; i < count; i++) {
    const item = candidates.nth(i);
    if (!(await item.isVisible()))
      continue;
    const box = await item.boundingBox();
    if (!box)
      continue;
    const cy = box.y + box.height / 2;
    await page.mouse.move(Math.max(0, box.x - 8), cy);
    await page.mouse.move(box.x + box.width / 2, cy);
    if (wanted === undefined)
      return;
    const opened = await expect.poll(() => menuShows(page, wanted), {timeout: pollMs(2000)}).toBe(true).then(() => true, () => false);
    if (opened) {
      // step into the flyout along the group's own row: the menu hides a submenu when the pointer
      // leaves the group item at a steep angle, which would take the item away mid-click
      const target = await menuLabelsMatching(page, wanted).filter({visible: true}).first().boundingBox();
      if (target != null) {
        await page.mouse.move(target.x + 4, cy);
        await page.mouse.move(target.x + 4, target.y + target.height / 2);
      }
      return;
    }
  }
  throw new Error(`the "${label}" group did not show "${wanted}"; the menu shows: ${await visibleMenuLabels(page)}`);
}

/** `Misc > Show Inside Values`: hovers the groups, clicks the leaf. The item ancestor is the usual
 * click target, but an inline group lays the box out on the label itself. */
export async function pickMenuPath(page: Page, path: string): Promise<void> {
  const segments = path.split(SEP).filter((s) => s.length > 0);
  for (let i = 0; i < segments.length; i++) {
    const last = i === segments.length - 1;
    const act = last
      ? menuItems(page, segments[i]).filter({visible: true}).first().click({timeout: 3000})
        .catch(() => menuLabelsMatching(page, segments[i]).filter({visible: true}).first().click({timeout: 3000}))
      : openGroup(page, segments[i], segments[i + 1]);
    await act.catch(async () => {
      throw new Error(`no "${segments[i]}" in the menu; it shows: ${await visibleMenuLabels(page)}`);
    });
  }
}

/** The labels of the open menu, opening every group of [path] first ("" reads the top level). */
export async function menuLabels(page: Page, path: string): Promise<string[]> {
  const segments = path.split(SEP).filter((s) => s.length > 0);
  for (let i = 0; i < segments.length; i++)
    await openGroup(page, segments[i], segments[i + 1]);
  return (await visibleMenuLabels(page)).split(' | ').filter(Boolean);
}
