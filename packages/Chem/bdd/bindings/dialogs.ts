/* A dialog button the way an impatient user meets it: the pointer lands on it as soon as it is on
   screen, whether or not it has become enabled (the library's click waits for an enabled target,
   which is what hides a button that is clickable for its first seconds). The button's state is
   recorded when it first appears and when the click lands, and what follows is watched over a
   window rather than read once. */
import {Page} from '@playwright/test';
import {Given, Then, When} from '@datagrok-libraries/bdd';
import {type ElementRef, el, expect, locate} from '@datagrok-libraries/bdd/runtime';
import {check, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {userSettingsPutBack} from '@datagrok-libraries/bdd/bindings/platform/steps';

declare global {
  interface Window {__bddOkStates?: {appeared?: boolean; clicked?: boolean; hit?: boolean; under?: string}}
}

export const watchNextOk = Given('user watches the OK button of the next dialog', async (page: Page) => {
  await page.evaluate(() => {
    const disabled = (e: Element) => (e as HTMLButtonElement).disabled ||
      e.getAttribute('aria-disabled') === 'true' || e.classList.contains('disabled');
    window.__bddOkStates = {};
    const seen = () => {
      const ok = document.querySelector('.d4-dialog [name="button-OK"]');
      if (!ok)
        return false;
      window.__bddOkStates!.appeared = disabled(ok);
      return true;
    };
    if (seen())
      throw new Error('a dialog with an OK button is already open');
    const observer = new MutationObserver(() => {
      if (seen())
        observer.disconnect();
    });
    observer.observe(document.body, {childList: true, subtree: true});
  });
}, {tier: 'api', description: 'records whether the OK button of the next dialog is disabled in the same task the button enters the page'});

export const clickAtOnce = When('user clicks on {element} at once', async (page: Page, target: ElementRef) => {
  const loc = await locate(page, target);
  await loc.waitFor({state: 'visible'});
  // a dialog still laying itself out moves its buttons: the point is taken from the element as it is
  // now, and what is under it is recorded with the state. A disabled button takes no pointer events
  // (ui.css), so the hit test falls through it to the bar it sits in: that is still a click on it,
  // anything else covering it is not.
  const state = await loc.evaluate((e) => {
    const r = e.getBoundingClientRect();
    const [x, y] = [r.x + r.width / 2, r.y + r.height / 2];
    const under = document.elementFromPoint(x, y);
    return {x, y, clicked: (e as HTMLButtonElement).disabled || e.getAttribute('aria-disabled') === 'true' || e.classList.contains('disabled'),
      hit: under !== null && (e.contains(under) || under.contains(e)),
      under: under ? `${under.tagName.toLowerCase()}.${under.className} "${(under.textContent ?? '').trim().slice(0, 30)}"` : 'nothing'};
  });
  await page.mouse.click(state.x, state.y);
  await page.evaluate((s) => { window.__bddOkStates = {...window.__bddOkStates, clicked: s.clicked, hit: s.hit, under: s.under}; }, state);
}, {tier: 'ui', description: 'a pointer click on the element the moment it is visible, without waiting for it to be enabled; its state at the click, and what is under the click point, are recorded'});

export const okDisabledThroughout = Then('the OK button should have been disabled when it appeared and when it was clicked', async (page: Page) => {
  const {under, ...states} = await page.evaluate(() => window.__bddOkStates ?? {});
  expect(states, `the OK button's state when it entered the page and when the click landed on it (under the click: ${under})`)
    .toEqual({appeared: true, clicked: true, hit: true});
});

// ---------------------------------------------------------------- the sketcher dialog's "Filter as you draw"

/* The substructure filter's "Filter as you draw" is the account's (the `chem-filter` user settings, on the server), and a
   feature leaves the stand as it found it (the library's hard rule): the first toggle in a feature remembers the whole map
   and puts it back at feature end, the server read back, as the library's `the "…" user settings are put back at feature
   end` does. These literal lines are the library's `user checks {element}` and `user unchecks {element}` for that one
   input, which they take over (more literal text wins). */
const FILTER_AS_YOU_DRAW = '"Filter as you draw" input in sketcher dialog';

export const checkFilterAsYouDraw = When('user checks "Filter as you draw" input in sketcher dialog', async (page: Page) => {
  await userSettingsPutBack(page, 'chem-filter');
  await check(page, el(FILTER_AS_YOU_DRAW));
}, {tier: 'ui', description: 'checks the option, the chem-filter settings of the account put back at feature end'});

export const uncheckFilterAsYouDraw = When('user unchecks "Filter as you draw" input in sketcher dialog', async (page: Page) => {
  await userSettingsPutBack(page, 'chem-filter');
  await uncheck(page, el(FILTER_AS_YOU_DRAW));
}, {tier: 'ui', description: 'unchecks the option, the chem-filter settings of the account put back at feature end'});
