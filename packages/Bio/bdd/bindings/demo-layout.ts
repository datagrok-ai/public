// Steps of the search and analysis demo clips (features/demos/15-18b): the view arranged off camera. A Bio command
// docks what it adds where it chooses (a search viewer under the grid, the cliffs table along the whole bottom of
// the view); a clip moves it beside the grid, or close to what the pointer works on next, before the next picture.
import type {Page} from '@playwright/test';
import {When} from '@datagrok-libraries/bdd';
import type {ElementRef} from '@datagrok-libraries/bdd/runtime';
import {silent, viewers} from '@datagrok-libraries/bdd/runtime';

type Side = 'right' | 'down';

/** Docks the viewer the phrase names (its element found among the current view's viewers) against the view's grid,
 * or against another viewer, taking that share of the space it splits. A viewer already docked is moved: the dock
 * manager takes it out of its old place (what Bio's Search commands do right after adding theirs). */
async function dockViewer(page: Page, target: ElementRef, side: Side, other: ElementRef | null, share: number): Promise<void> {
  silent(page);
  const mark = 'data-bdd-dock';
  const tag = async (ref: ElementRef, value: string): Promise<void> =>
    viewers.onViewer(page, ref, (el, v) => { el.setAttribute('data-bdd-dock', v as string); }, value);
  await tag(target, 'moved');
  if (other)
    await tag(other, 'anchor');
  try {
    const why = await page.evaluate(([s, ratio, m]) => {
      const tv = (window as any).grok.shell.tv;
      if (!tv)
        return 'no table view is current';
      const owner = (el: Element | null) => el == null ? null :
        Array.from(tv.viewers as Iterable<any>).find((v: any) => v.root === el || v.root.contains(el) || el.contains(v.root));
      const moved = owner(document.querySelector(`[${m}="moved"]`));
      if (!moved)
        return 'the viewer is not a viewer of the current view';
      const anchorEl = document.querySelector(`[${m}="anchor"]`);
      const anchor = anchorEl ? owner(anchorEl) : null;
      if (anchorEl && !anchor)
        return 'the other viewer is not a viewer of the current view';
      const node = tv.dockManager.findNode((anchor ?? tv.grid).root);
      if (!node)
        return `no dock node holds ${anchor ? 'the other viewer' : 'the grid'}`;
      tv.dockManager.dock(moved, s, node, undefined, ratio);
      return '';
    }, [side, share, mark] as [Side, number, string]);
    if (why)
      throw new Error(`${target.phrase}: ${why}`);
  } finally {
    await page.evaluate((m) => document.querySelectorAll(`[${m}]`).forEach((e) => e.removeAttribute(m)), mark).catch(() => undefined);
  }
  await viewers.settleAll(page);
}

export const dockedRightOfGrid = When('{widget} is docked to the right of the grid, {int}% of the view wide',
  (page: Page, target: ElementRef, percent: number) => dockViewer(page, target, 'right', null, percent / 100),
  {tier: 'api', description: 'set-up a filmed demo leaves out: the viewer moved beside the view\'s grid, on its right, taking that share of the width'});

export const dockedUnderGrid = When('{widget} is docked under the grid, {int}% of the view tall',
  (page: Page, target: ElementRef, percent: number) => dockViewer(page, target, 'down', null, percent / 100),
  {tier: 'api', description: 'set-up a filmed demo leaves out: the viewer moved under the view\'s grid (not under the whole view), taking that share of the height'});

export const dockedRightOfViewer = When('{widget} is docked to the right of {widget}',
  (page: Page, target: ElementRef, other: ElementRef) => dockViewer(page, target, 'right', other, 0.5),
  {tier: 'api', description: 'set-up a filmed demo leaves out: the viewer moved beside the other one, on its right, the two sharing its width'});

export const dockedUnderViewer = When('{widget} is docked under {widget}',
  (page: Page, target: ElementRef, other: ElementRef) => dockViewer(page, target, 'down', other, 0.5),
  {tier: 'api', description: 'set-up a filmed demo leaves out: the viewer moved under the other one, the two sharing its height'});
