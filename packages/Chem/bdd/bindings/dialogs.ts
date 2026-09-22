/* A dialog button the way an impatient user meets it: the pointer lands on it as soon as it is on
   screen, whether or not it has become enabled (the library's click waits for an enabled target,
   which is what hides a button that is clickable for its first seconds). The button's state is
   recorded when it first appears and when the click lands, and what follows is watched over a
   window rather than read once. */
import {Page} from '@playwright/test';
import {Given, Then, When} from '@datagrok-libraries/bdd';
import {type ElementRef, expect, locate, takeErrors, viewers} from '@datagrok-libraries/bdd/runtime';

declare global {
  interface Window {__bddOkStates?: {appeared?: boolean; clicked?: boolean}}
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
  const box = await loc.boundingBox();
  if (!box)
    throw new Error(`${target.phrase}: no bounding box to click`);
  const state = await loc.evaluate((e) => (e as HTMLButtonElement).disabled ||
    e.getAttribute('aria-disabled') === 'true' || e.classList.contains('disabled'));
  await page.mouse.click(box.x + box.width / 2, box.y + box.height / 2);
  await page.evaluate((s) => { window.__bddOkStates = {...window.__bddOkStates, clicked: s}; }, state);
}, {tier: 'ui', description: 'a pointer click on the element the moment it is visible, without waiting for it to be enabled; its state at the click is recorded'});

export const okDisabledThroughout = Then('the OK button should have been disabled when it appeared and when it was clicked', async (page: Page) => {
  const states = await page.evaluate(() => window.__bddOkStates ?? {});
  expect(states, 'the OK button\'s state when it entered the page and when the click landed').toEqual({appeared: true, clicked: true});
});

export const quietWindow = Then('no error or warning balloon and no error should appear for {int} seconds', async (page: Page, seconds: number) => {
  const seen: string[] = [];
  const end = Date.now() + seconds * 1000;
  while (Date.now() < end) {
    seen.push(...(await viewers.takeBalloons(page)).filter((b) => b.type === 'error' || b.type === 'warning').map((b) => `${b.type}: ${b.message}`));
    seen.push(...takeErrors(page));
    if (seen.length > 0)
      break;
    await page.waitForTimeout(200);
  }
  expect(seen, `error and warning balloons and errors over ${seconds} s`).toEqual([]);
}, {description: 'a zero held over the window, not read once: a balloon or error that lands a few tasks after the gesture is caught'});
