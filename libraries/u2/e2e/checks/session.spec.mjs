/* The shared session in a real browser, platform-free: the gallery's `#session` page over the
   memory backend — two sources under one `SharedSession`, an edit in each counted together and
   saved as one batch, a draft parent referenced by a child draft landing in one Save, the
   unsaved-changes gate in front of a project switch, the search box narrowing the child rows.
   Gallery-hosted (`e2e/gallery-server.mjs`), so it ends by rebooting the client for the suites
   after it. */
import {BOOT_TIMEOUT, LOCAL_URL, consoleErrors, ok, pageErrors, shot} from '../local.mjs';
import {reopenApp} from '../lib.mjs';
import {startGalleryServer} from '../gallery-server.mjs';

let server = null;

const LIST = '[data-u2="list"][data-u2-name="projects"]';
const projects = (page) => page.locator(`${LIST} .u2-list-row`);
const issues = (page) => page.locator('[data-u2="table"][data-u2-name="issues"] tbody tr');
const summary = (page) => page.locator('[data-u2-name="summary"]').textContent();
const input = (page, name) => page.locator(`.u2-domain-form [data-u2-name="${name}"] input`);
const saveButton = (page) => page.locator('[data-u2="save-button"]');
const dialogButton = (page, text) => page.locator('.u2-dialog button', {hasText: text});
const errors = () => consoleErrors.length + pageErrors.length;

export async function fixture(page) {
  server = await startGalleryServer();
  await page.goto(`${server.url}/gallery/#session`, {waitUntil: 'load'});
  await page.waitForSelector(`${LIST} .u2-list-row`, {timeout: 30000});
  await page.waitForSelector('[data-u2="table"][data-u2-name="issues"] tbody tr', {timeout: 30000});
}

async function checkLoad(page) {
  ok('session/1a/the projects render from the memory backend', await projects(page).count() === 3);
  ok('session/1b/the first project is current and its issues are shown', await issues(page).count() === 4,
    `issues=${await issues(page).count()}`);
  ok('session/1c/Save is disabled while clean', await saveButton(page).isDisabled());
  await shot(page, 'session-1-loaded');
}

async function checkOneSaveTwoTables(page) {
  await input(page, 'name').first().fill('Grit renamed');
  await issues(page).nth(0).click();
  await page.waitForSelector('.u2-domain-form [data-u2-name="title"] input');
  await input(page, 'title').first().fill('Issue #1: renamed through the form');
  await page.waitForTimeout(100);
  ok('session/2a/the summary counts both tables', (await summary(page)).includes('2 unsaved changes in 2 tables'),
    await summary(page));
  await saveButton(page).click();
  await page.waitForTimeout(300);
  ok('session/2b/one Save landed both', await saveButton(page).isDisabled() &&
    (await projects(page).nth(0).textContent()).includes('Grit renamed') &&
    (await issues(page).nth(0).textContent()).includes('renamed through the form'));
  await shot(page, 'session-2-saved');
}

async function checkDraftParentAndChild(page) {
  await page.locator('[data-u2-name="newProject"]').click();
  await page.waitForTimeout(100);
  ok('session/3a/a draft project is current and pristine', await projects(page).count() === 4 &&
    await saveButton(page).isDisabled());
  await page.locator('[data-u2-name="newIssue"]').click();
  await page.waitForTimeout(100);
  await input(page, 'title').first().fill('First issue of the new project');
  await input(page, 'key').first().fill('NEW');
  await input(page, 'name').first().fill('New project');
  await page.waitForTimeout(100);
  ok('session/3b/the child refers to the draft parent by its name', (await summary(page)).includes('in 2 tables') &&
    (await page.locator('.u2-domain-form [data-u2-name="project_id"] input').first().inputValue()).includes('New project'));
  await saveButton(page).click();
  await page.waitForTimeout(300);
  ok('session/3c/parent and child landed as one transaction', await saveButton(page).isDisabled() &&
    await issues(page).count() === 1 && (await issues(page).nth(0).textContent()).includes('First issue'),
    `issues=${await issues(page).count()} · ${await summary(page)}`);
  await shot(page, 'session-3-draft-saved');
}

async function checkGate(page) {
  await issues(page).nth(0).click();
  await input(page, 'title').first().fill('An edit to abandon');
  await page.waitForTimeout(100);
  await projects(page).nth(0).click();
  await page.waitForSelector('.u2-dialog', {timeout: 5000});
  ok('session/4a/switching projects while dirty asks first', await dialogButton(page, 'DISCARD').count() === 1);
  await dialogButton(page, 'DISCARD').click();
  await page.waitForTimeout(300);
  ok('session/4b/Discard drops the edit and the switch goes through', await saveButton(page).isDisabled() &&
    await issues(page).count() === 4, `issues=${await issues(page).count()}`);
}

async function checkSearch(page) {
  await page.locator('[data-u2="domain-search"] input').fill('#4');
  await page.keyboard.press('Enter');
  await page.waitForTimeout(300);
  ok('session/5a/the search narrows the issues within the project', await issues(page).count() === 1,
    `issues=${await issues(page).count()}`);
  ok('session/5b/no console or page errors', errors() === 0, [...consoleErrors, ...pageErrors].join(' | ').slice(0, 300));
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
  ok('session/6/the designer is back for the suites that follow', true);
}

export const checks = [
  {id: 'session/1 load', run: checkLoad},
  {id: 'session/2 one save, two tables', run: checkOneSaveTwoTables},
  {id: 'session/3 draft parent and child', run: checkDraftParentAndChild},
  {id: 'session/4 gate', run: checkGate},
  {id: 'session/5 search', run: checkSearch},
  {id: 'session/6 restore', run: restoreClient},
];
