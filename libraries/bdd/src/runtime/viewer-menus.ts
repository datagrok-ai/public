import {Locator, Page} from '@playwright/test';
import {expect, pollMs} from './patience.js';
import type {ElementRef} from './args.js';
import {exactText} from './locate.js';
import {evaluate} from './viewer-runtime.js';
import {onReachable, onViewer} from './viewers.js';
import * as guide from './guide.js';

/** A Dart popup menu, or a u2 one (a u2 submenu is a `.u2-menu` inside its group item, as a Dart group's children
 * are inside its item). */
const POPUP = '.d4-menu-popup, .u2-menu';
const LABEL = '.d4-menu-item-label, .u2-menu-label';
const OWN_LABEL = ':scope > .d4-menu-item-label, :scope > .u2-menu-label';
const ITEM = '.d4-menu-item, .u2-menu-item';
const SEP = /\s*[>|]\s*/;
/** A menu is shown once its items are built, and a tree node's may wait for the server: under load that
 * took past 3 s; the wait ends as the menu shows, so the cap costs only a real failure. */
const MENU_SHOWN_MS = (): number => pollMs(10000);
const MENU_ITEM = 'xpath=ancestor::*[contains(concat(" ", normalize-space(@class), " "), " d4-menu-item ") or ' +
  'contains(concat(" ", normalize-space(@class), " "), " u2-menu-item ")][1]';

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
  if (!await page.evaluate((t) => (window as any).__bdd.waitArmed(t), token)) {
    const under = await page.evaluate(([px, py]) => {
      const e = document.elementFromPoint(px, py);
      return e ? `${e.tagName.toLowerCase()}${e.getAttribute('name') ? `[name="${e.getAttribute('name')}"]` : ''}.${[...e.classList].join('.')} "${(e.textContent ?? '').trim().slice(0, 40)}"` : 'nothing';
    }, [x, y] as [number, number]).catch(() => 'unknown');
    throw new Error(`no context menu opened at (${Math.round(x)}, ${Math.round(y)}); under the point now: ${under}`);
  }
  return page.locator(POPUP).last();
}

export async function openContextMenuAt(page: Page, x: number, y: number): Promise<Locator> {
  const token: string = await evaluate(page, (cap) => (window as any).__bdd.openMenu(cap), MENU_SHOWN_MS());
  return rightClickArmed(page, x, y, token);
}

/** The context menu of an element: at a named hit area, else its `view` area when it reports
 * one, else its centre. Three roundtrips: the point and the arming in one, the click, the wait.
 * A list that reflows between the point and the click (a gallery relaid out, a tree refilled) puts
 * another element under it, whose menu opens instead: the point also arms a check of what the
 * right-click reached, and a click that missed the element is aimed again. */
export async function openContextMenuOf(page: Page, target: ElementRef, area?: string): Promise<Locator> {
  for (let attempt = 1; ; attempt++) {
    const {x, y, token} = await onReachable(page, target, (el, [a, cap]) => {
      const w = window as any;
      w.__bddMenuHit = undefined;
      w.__bddMenuReached = 'no contextmenu event';
      document.addEventListener('contextmenu', (e) => {
        const t = e.target as Element;
        w.__bddMenuHit = el.contains(t);
        w.__bddMenuReached = `${t.tagName?.toLowerCase()}${t.getAttribute?.('name') ? `[name="${t.getAttribute('name')}"]` : ''} ` +
          `"${(t.textContent ?? '').trim().slice(0, 40)}"${w.__bddMenuHit ? '' : ', outside the target'}`;
      }, {capture: true, once: true});
      return w.__bdd.menuPoint(el, a, cap);
    }, [area ?? null, MENU_SHOWN_MS()] as [string | null, number]);
    const menu = await rightClickArmed(page, x, y, token);
    if (attempt === 3 || await page.evaluate(() => (window as any).__bddMenuHit !== false)) {
      await page.evaluate((n) => (window as any).__bddMenuReached += ` (try ${n})`, attempt).catch(() => undefined);
      return menu;
    }
    await closeContextMenu(page);
  }
}

/** A group item holds its children under its own label: the labels inside it are the children's. */
const CHILDREN = '.d4-menu-item-container, .u2-menu';

/** Every label of the open menus reading exactly that — a viewer menu can hold two groups of the
 * same name (the axis "Annotations" and the viewer's own) — or, inside a group, of its children. */
export function menuLabelsMatching(page: Page, label: string, within?: Locator): Locator {
  return (within ? within.locator(CHILDREN) : page.locator(POPUP)).locator(LABEL, {hasText: exactText(label)});
}

export function menuItems(page: Page, label: string, within?: Locator): Locator {
  return menuLabelsMatching(page, label, within).locator(MENU_ITEM);
}

// a group's children live in the same popup, hidden until it opens: shown is any visible match
async function menuShows(page: Page, label: string, within?: Locator): Promise<boolean> {
  return (await menuLabelsMatching(page, label, within).filter({visible: true}).count()) > 0;
}

export function menuItem(page: Page, label: string): Locator {
  return menuItems(page, label).first();
}

/** What the open menus show, by item: its label, or the caption of an element item (a property
 * editor the menu embeds — the PC plot's column picker under "Columns" — carries it as `data-source`). */
async function visibleMenuLabels(page: Page, within?: Locator): Promise<string> {
  const items = (within ? within.locator(CHILDREN) : page.locator(POPUP)).locator(ITEM).filter({visible: true});
  const labels = await items.evaluateAll((els, own) => els.map((e) => e.querySelector(own)?.textContent ?? e.getAttribute('data-source') ?? ''), OWN_LABEL);
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
    if (await menuShows(page, wanted, candidates.nth(i))) {
      await guide.hop(page, candidates.nth(i));
      return candidates.nth(i);
    }
  }
  for (let i = 0; i < count; i++) {
    const item = candidates.nth(i);
    if (!(await item.isVisible()))
      continue;
    const box = await item.boundingBox();
    if (!box)
      continue;
    await guide.hop(page, item);
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
  // a context menu's own: what its right-click reached, when the item is not in it
  const reached = async () => {
    const r: string | undefined = await page.evaluate(() => {
      const w = window as any;
      const value = w.__bddMenuReached;
      delete w.__bddMenuReached;
      return value;
    }).catch(() => undefined);
    return r ? `; the right-click reached ${r}` : '';
  };
  let group: Locator | undefined;
  for (let i = 0; i < segments.length; i++) {
    const last = i === segments.length - 1;
    if (!last) {
      group = await openGroup(page, segments[i], segments[i + 1], group).catch(async () => {
        throw new Error(`no "${segments[i]}" in the menu; it shows: ${await visibleMenuLabels(page, group)}${await reached()}`);
      });
      continue;
    }
    const label = menuLabelsMatching(page, segments[i], group).filter({visible: true}).first();
    await label.waitFor({state: 'visible', timeout: 3000}).catch(async () => {
      throw new Error(`no "${segments[i]}" in the menu; it shows: ${await visibleMenuLabels(page, group)}${await reached()}`);
    });
    const item = label.locator(MENU_ITEM);
    const leaf = await item.count() > 0 ? item : label;
    await guide.hop(page, leaf);
    await leaf.click();
  }
  await page.evaluate(() => delete (window as any).__bddMenuReached).catch(() => undefined);
}

/** The labels under [path], opening every group of it first ("" reads the top level). */
/** Whether the last item of [path] is `enabled` or `disabled` (its `aria-disabled`), the groups before it
 * opened; `not listed` while the menu does not show it. */
export async function menuItemState(page: Page, path: string): Promise<string> {
  const segments = path.split(SEP).filter((s) => s.length > 0);
  const item = segments.pop()!;
  let group: Locator | undefined;
  for (let i = 0; i < segments.length; i++)
    group = await openGroup(page, segments[i], segments[i + 1] ?? item, group);
  const items = (group ? group.locator(CHILDREN) : page.locator(POPUP)).locator(ITEM).filter({visible: true});
  const states = await items.evaluateAll((els, [label, own]) => els
    .filter((e) => (e.querySelector(own)?.textContent ?? '').trim() === label)
    .map((e) => e.getAttribute('aria-disabled') === 'true' ? 'disabled' : 'enabled'), [item, OWN_LABEL] as [string, string]);
  return states[0] ?? 'not listed';
}

export async function menuLabels(page: Page, path: string): Promise<string[]> {
  const segments = path.split(SEP).filter((s) => s.length > 0);
  let group: Locator | undefined;
  for (let i = 0; i < segments.length; i++)
    group = await openGroup(page, segments[i], segments[i + 1], group);
  return (await visibleMenuLabels(page, group)).split(' | ').filter(Boolean);
}
