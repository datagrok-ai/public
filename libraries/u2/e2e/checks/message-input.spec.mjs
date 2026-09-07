/* Real-browser caret and chip behaviour of the compose box, against the gallery's static provider:
   contenteditable, selection and Range are exactly what the headless suite cannot stand in for.
   Gallery-hosted (`e2e/gallery-server.mjs`), so it ends by reopening the app for the suites after it. */
import {BOOT_TIMEOUT, LOCAL_URL, ok, shot} from '../local.mjs';
import {reopenApp} from '../lib.mjs';
import {startGalleryServer} from '../gallery-server.mjs';

const EDITOR = '.u2-msg .u2-msg-editor';
const POPUP = '.u2-msg-popup';

let server;

export async function fixture(page) {
  server = await startGalleryServer();
  await page.goto(`${server.url}/gallery/#message-input`);
  await page.waitForSelector(EDITOR, {timeout: 30000});
}

const editor = (page) => page.locator(EDITOR).first();
const value = (page) => page.locator('.u2-gallery-status code').first().textContent();

async function typeInEditor(page, text) {
  await editor(page).click();
  await page.keyboard.type(text, {delay: 20});
}

async function checkMentionFlow(page) {
  await typeInEditor(page, 'hello @al');
  await page.waitForSelector(POPUP, {timeout: 5000});
  const query = await page.locator('.u2-msg-query').first().boundingBox();
  const popup = await page.locator(POPUP).first().boundingBox();
  const box = await editor(page).boundingBox();
  ok('message-input/1a/popup-is-anchored-to-the-caret-not-the-box',
    Math.abs(popup.x - query.x) < 8 && popup.x > box.x + 8,
    `query.x=${query?.x} popup.x=${popup?.x} editor.x=${box?.x}`);

  await page.keyboard.press('ArrowDown');
  await page.keyboard.press('Enter');
  await page.waitForTimeout(200);
  const chips = await page.locator('.u2-msg-chip').count();
  const editable = await page.locator('.u2-msg-chip').first().getAttribute('contenteditable');
  await shot(page, 'message-input-1-chip-inserted');
  ok('message-input/1b/pick-inserts-one-atomic-chip', chips === 1 && editable === 'false',
    `chips=${chips} contenteditable=${editable}`);

  await page.keyboard.type(', ping');
  await page.waitForTimeout(200);
  const text = await editor(page).textContent();
  const serialized = await value(page);
  ok('message-input/1c/typing-continues-after-the-chip-and-the-value-carries-the-token',
    text.endsWith(', ping') && /<span>#\{x\..+?\."[^"]+"\}<\/span>/.test(serialized),
    `text="${text}" value="${serialized}"`);
}

async function checkBackspaceEatsTheChip(page) {
  const before = await page.locator('.u2-msg-chip').count();
  for (let i = 0; i < 6; i++)
    await page.keyboard.press('Backspace');
  await page.waitForTimeout(100);
  const midway = await page.locator('.u2-msg-chip').count();
  await page.keyboard.press('Backspace');
  await page.keyboard.press('Backspace');
  await page.waitForTimeout(200);
  const after = await page.locator('.u2-msg-chip').count();
  await shot(page, 'message-input-2-chip-deleted');
  ok('message-input/2a/backspace-deletes-the-chip-whole',
    before === 1 && midway === 1 && after === 0,
    `chips before=${before} after-text=${midway} after-chip=${after}`);
}

async function checkEnterModes(page) {
  const quick = page.locator('.u2-msg .u2-msg-editor').nth(1);
  await quick.click();
  await page.keyboard.type('first');
  await page.keyboard.press('Shift+Enter');
  await page.keyboard.type('second');
  const lines = (await quick.textContent()).split('\n').length;
  await page.keyboard.press('Enter');
  await page.waitForTimeout(200);
  const cleared = await quick.textContent();
  ok('message-input/3a/shift-enter-breaks-the-line-and-enter-sends',
    lines === 2 && cleared.trim() === '', `lines=${lines} left="${cleared}"`);
}

/** The first gallery box carries a `draftKey`, so what is left unsent in it is what a reload has to
 * bring back — text and chip, the chip's token byte-for-byte. */
async function checkDraftSurvivesReload(page) {
  await page.locator('#main').getByRole('button', {name: 'Clear', exact: true}).first().click();
  await typeInEditor(page, 'ship it @bru');
  await page.waitForSelector(POPUP, {timeout: 5000});
  await page.keyboard.press('ArrowDown');
  await page.keyboard.press('Enter');
  await page.waitForTimeout(200);
  const before = await page.locator('.u2-msg-chip').first().getAttribute('data-token');

  await page.reload();
  await page.waitForSelector(EDITOR, {timeout: 30000});
  const text = await editor(page).textContent();
  const chips = await page.locator('.u2-msg-chip').count();
  const after = chips > 0 ? await page.locator('.u2-msg-chip').first().getAttribute('data-token') : null;
  await shot(page, 'message-input-4-draft-restored');
  ok('message-input/4a/the-unsent-text-and-its-chip-come-back-after-a-reload',
    text.startsWith('ship it ') && chips === 1 && after === before && before !== null,
    `text="${text}" chips=${chips} token before=${before} after=${after}`);
  await page.evaluate(() => localStorage.removeItem('u2-msg-draft-e2e'));

  await restoreClient(page);
}

/** The gallery lives on its own origin, so the client has to be re-booted, not just reopened. */
async function restoreClient(page) {
  await server.close();
  await page.goto(`${LOCAL_URL}/login.html?mode=local`, {waitUntil: 'load', timeout: BOOT_TIMEOUT});
  await page.waitForFunction(() => {
    try {
      return !!(window.DG && DG.Func && grok.shell.user);
    } catch (e) {
      return false;
    }
  }, null, {timeout: BOOT_TIMEOUT});
  await reopenApp(page);
}

export const checks = [
  {id: 'message-input/1 @ opens a caret-anchored popup and picks a chip', run: checkMentionFlow},
  {id: 'message-input/2 backspace deletes a chip whole', run: checkBackspaceEatsTheChip},
  {id: 'message-input/3 enter modes', run: checkEnterModes},
  {id: 'message-input/4 draft persistence', run: checkDraftSurvivesReload},
];
