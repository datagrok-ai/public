/* Tour geometry in a real browser — the one thing the headless shim cannot answer: the SVG mask
   hole really covers the target, the popup really sits clear of it, and both really follow a
   scroll. Runs against the static gallery (`#tour`), not the designer, so the tour is exercised
   platform-free; the last check brings the client back for the suites that follow. */
import {BOOT_TIMEOUT, LOCAL_URL, ok, shot} from '../local.mjs';
import {reopenApp} from '../lib.mjs';
import {startGalleryServer} from '../gallery-server.mjs';

let server = null;

const state = (page) => page.evaluate(() => {
  const rect = (el) => {
    const b = el.getBoundingClientRect();
    return {x: b.x, y: b.y, width: b.width, height: b.height};
  };
  const hole = document.querySelector('.u2-tour-hole');
  const popup = document.querySelector('.u2-tour-popup');
  const target = document.querySelector('[data-u2-name="nameInput"]');
  return {
    overlay: !!document.querySelector('.u2-tour-overlay'),
    hole: hole && {x: Number(hole.getAttribute('x')), y: Number(hole.getAttribute('y')),
      width: Number(hole.getAttribute('width')), height: Number(hole.getAttribute('height'))},
    popup: popup && rect(popup),
    target: target && rect(target),
    counter: document.querySelector('.u2-tour-counter')?.textContent ?? '',
    focus: document.activeElement?.className ?? '',
  };
});

const frames = (page) => page.evaluate(() =>
  new Promise((resolve) => requestAnimationFrame(() => requestAnimationFrame(resolve))));

const near = (a, b, tolerance = 6) => Math.abs(a - b) <= tolerance;
const overlap = (a, b) => Math.max(0, Math.min(a.x + a.width, b.x + b.width) - Math.max(a.x, b.x)) *
  Math.max(0, Math.min(a.y + a.height, b.y + b.height) - Math.max(a.y, b.y));

async function startTour(page) {
  await page.getByRole('button', {name: 'Start tour'}).click();
  await page.waitForSelector('.u2-tour-popup', {timeout: 5000});
  await frames(page);
}

export async function fixture(page) {
  server = await startGalleryServer();
  await page.goto(`${server.url}/gallery/#tour`, {waitUntil: 'load'});
  await page.waitForSelector('[data-u2-name="nameInput"]', {timeout: 15000});
}

async function checkSpotlight(page) {
  await startTour(page);
  const {hole, popup, target, focus} = await state(page);
  await shot(page, 'tour-1-spotlight');
  const fits = near(hole.x, target.x - 4) && near(hole.y, target.y - 4) &&
    near(hole.width, target.width + 8) && near(hole.height, target.height + 8);
  ok('tour/1a/the hole covers the target with 4px of padding', fits,
    `hole=${JSON.stringify(hole)} target=${JSON.stringify(target)}`);
  ok('tour/1b/the popup stays clear of the spotlight', overlap(popup, hole) === 0,
    `popup=${JSON.stringify(popup)} hole=${JSON.stringify(hole)}`);
  ok('tour/1c/focus lands in the popup', focus.includes('u2-tour-next'), `activeElement="${focus}"`);
}

/** The gallery page fits the e2e viewport, so it is shrunk first — which also proves the resize
 * path — and only then scrolled. */
async function checkScroll(page) {
  const size = page.viewportSize();
  await page.setViewportSize({width: size.width, height: 420});
  await frames(page);
  const before = await state(page);
  ok('tour/2a/a resize re-anchors the hole', near(before.hole.y, before.target.y - 4),
    `hole=${JSON.stringify(before.hole)} target=${JSON.stringify(before.target)}`);

  await page.evaluate(() => document.getElementById('main').scrollBy(0, 200));
  await frames(page);
  const after = await state(page);
  const moved = before.hole.y - after.hole.y;
  ok('tour/2b/the hole follows the target through a scroll',
    moved > 0 && near(moved, before.target.y - after.target.y) && near(after.hole.y, after.target.y - 4),
    `moved=${moved} hole=${JSON.stringify(after.hole)} target=${JSON.stringify(after.target)}`);
  ok('tour/2c/the popup is re-anchored too and still clear',
    near(before.popup.y - after.popup.y, moved, 12) && overlap(after.popup, after.hole) === 0,
    `popup ${before.popup.y} -> ${after.popup.y}, hole moved ${moved}`);
  await page.setViewportSize(size);
  await frames(page);
}

async function checkSkipForward(page) {
  await page.locator('.u2-tour-next').click();
  await frames(page);
  const second = await state(page);
  await page.locator('.u2-tour-next').click();
  await frames(page);
  const fourth = await state(page);
  await shot(page, 'tour-3-skip-forward');
  ok('tour/3a/NEXT walks the steps', second.counter === '2 / 4', `counter="${second.counter}"`);
  ok('tour/3b/the step whose target does not exist is skipped forward',
    fourth.counter === '4 / 4', `counter="${fourth.counter}"`);
}

/** The last step spotlights Run and advances on the signal a click on it sets — which only works
 * because the dim layer never takes the pointer. */
async function checkClickThrough(page) {
  await page.locator('[data-u2-name="runButton"]').click();
  await frames(page);
  const after = await state(page);
  const result = await page.evaluate(() => document.getElementById('tour-result')?.textContent ?? '');
  ok('tour/4a/the dim layer does not block the spotlit control', result.includes('done'),
    `readout="${result}"`);
  ok('tour/4b/advanceOn finished the tour and cleared the layer', !after.overlay && after.popup === null,
    `overlay=${after.overlay} popup=${after.popup === null ? 'gone' : 'present'}`);
}

async function checkEscape(page) {
  await startTour(page);
  await page.keyboard.press('Escape');
  await frames(page);
  const after = await state(page);
  const result = await page.evaluate(() => document.getElementById('tour-result')?.textContent ?? '');
  ok('tour/5a/Esc ends the tour and removes the overlay', !after.overlay && after.popup === null,
    `overlay=${after.overlay} popup=${after.popup === null ? 'gone' : 'present'}`);
  ok('tour/5b/onFinish reports the skip', result.includes('skipped'), `readout="${result}"`);
}

/** The gallery lives on its own origin, so the client has to be re-booted, not just reopened. */
async function restoreClient(page) {
  await server?.close();
  server = null;
  await page.goto(`${LOCAL_URL}/login.html?mode=local`, {waitUntil: 'load', timeout: BOOT_TIMEOUT});
  await page.waitForFunction(() => {
    try {
      return !!(window.DG && DG.Func && grok.shell.user);
    } catch (e) {
      return false;
    }
  }, null, {timeout: BOOT_TIMEOUT});
  await reopenApp(page);
  ok('tour/6/the designer is back for the suites that follow', true);
}

export const checks = [
  {id: 'tour/1 spotlight geometry and focus', run: checkSpotlight},
  {id: 'tour/2 resize and scroll re-anchor both layers', run: checkScroll},
  {id: 'tour/3 NEXT walks the steps and skips the missing target', run: checkSkipForward},
  {id: 'tour/4 advanceOn through the click-through layer', run: checkClickThrough},
  {id: 'tour/5 Esc ends the tour', run: checkEscape},
  {id: 'tour/6 back to the client', run: restoreClient},
];
