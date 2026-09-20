/* The domain seam in a real browser, platform-free: the gallery's `#domain` page over the memory
   backend — rows render, selecting one puts it in the form, an edit through the form arms Save,
   Save writes the batch and the list shows it, a row the caller may not edit is text, Delete is
   offered only where the row says so, a draft is pristine until typed into. Gallery-hosted
   (`e2e/gallery-server.mjs`), so it ends by rebooting the client for the suites after it. */
import {BOOT_TIMEOUT, LOCAL_URL, consoleErrors, ok, pageErrors, shot} from '../local.mjs';
import {reopenApp} from '../lib.mjs';
import {startGalleryServer} from '../gallery-server.mjs';

let server = null;

const LIST = '[data-u2="list"][data-u2-name="tasks"]';
const rows = (page) => page.locator(`${LIST} .u2-list-row`);
const summary = (page) => page.locator('[data-u2-name="summary"]').textContent();
const input = (page, name) => page.locator(`.u2-domain-form [data-u2-name="${name}"] input`).first();
const saveButton = (page) => page.locator('[data-u2="save-button"]');
const errors = () => consoleErrors.length + pageErrors.length;

export async function fixture(page) {
  server = await startGalleryServer();
  await page.goto(`${server.url}/gallery/#domain`, {waitUntil: 'load'});
  await page.waitForSelector(`${LIST} .u2-list-row`, {timeout: 30000});
}

async function checkLoad(page) {
  const count = await rows(page).count();
  const text = await summary(page);
  ok('domain/1a/rows render from the memory backend', count > 0 && count <= 15, `rows=${count}`);
  ok('domain/1b/the summary counts the page against the total', text.includes('15 of 40'), text);
  ok('domain/1c/Save is disabled while clean', await saveButton(page).isDisabled());
  await shot(page, 'domain-1-loaded');
}

async function checkEditAndSave(page) {
  await rows(page).nth(0).click();
  await page.waitForSelector('.u2-domain-form [data-u2-name="title"] input');
  const title = input(page, 'title');
  ok('domain/2a/the selected row is the form\'s row', (await title.inputValue()).startsWith('Task #1:'));
  await title.fill('Task #1: renamed through the form');
  await page.waitForTimeout(100);
  ok('domain/2b/an edit arms Save', !(await saveButton(page).isDisabled()));
  ok('domain/2c/the summary counts the change', (await summary(page)).includes('1 change'));
  await saveButton(page).click();
  await page.waitForTimeout(200);
  ok('domain/2d/Save wrote the batch and the list shows it',
    (await rows(page).nth(0).textContent()).includes('renamed through the form') && await saveButton(page).isDisabled());
  await shot(page, 'domain-2-saved');
}

async function checkAccess(page) {
  await rows(page).nth(3).click();
  await page.waitForTimeout(100);
  const readonly = await page.locator('.u2-domain-form [data-u2="readonly-field"][data-u2-name="title"]').count();
  ok('domain/3a/a row carrying ~can_edit = false is text', readonly === 1, `readonly rows=${readonly}`);
  const deletable = await rows(page).nth(0).locator('.u2-row-actions button').count();
  const kept = await rows(page).nth(1).locator('.u2-row-actions button').count();
  ok('domain/3b/Delete follows the row\'s ~can_delete', deletable === 1 && kept === 0, `even=${deletable} odd=${kept}`);
}

async function checkDraft(page) {
  await page.locator('[data-u2-name="newTask"]').click();
  await page.waitForTimeout(100);
  ok('domain/4a/a draft is pristine', await saveButton(page).isDisabled() && (await summary(page)).includes('0 change'));
  await input(page, 'title').fill('A brand new task');
  await page.waitForTimeout(100);
  ok('domain/4b/typing into the draft arms Save', !(await saveButton(page).isDisabled()));
  await saveButton(page).click();
  await page.waitForTimeout(200);
  ok('domain/4c/the draft was inserted', (await summary(page)).includes('of 41'), await summary(page));
  ok('domain/4d/no console or page errors', errors() === 0,
    [...consoleErrors, ...pageErrors].join(' | ').slice(0, 300));
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
  ok('domain/5/the designer is back for the suites that follow', true);
}

export const checks = [
  {id: 'domain/1 load', run: checkLoad},
  {id: 'domain/2 edit and save', run: checkEditAndSave},
  {id: 'domain/3 access per row', run: checkAccess},
  {id: 'domain/4 draft', run: checkDraft},
  {id: 'domain/5 restore', run: restoreClient},
];
