/* The steps only Tutorials defines: starting a tutorial from its card, the running tutorial's
   progress and completion, the completion record it writes, and closing it. Everything a tutorial
   asks the learner to do is ordinary platform UI and uses the library's vocabulary; step entries,
   cards and the panel are named in elements.ts. */
import {type Page, test} from '@playwright/test';
import {Given, Then, When} from '@datagrok-libraries/bdd';
import {type ElementRef, atFeatureEnd, expect, locate, pollMs, reportedServices, serviceGap} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;
declare const DG: any;

const DATA_STORAGE_KEY = 'tutorials';

/** Closes the running tutorial through its own Close button (which rejects the step it waits on and
 * closes its table) and the Tutorials panel with it: the panel is not a view, so the shell reset
 * leaves it — and a tutorial left running keeps listening to the platform's events and would take
 * the next feature's gestures as its steps. */
const TUTORIAL_CLOSE = '.tutorials-root-header button:has(.fa-times-circle)';

async function closeTutorials(page: Page): Promise<void> {
  await page.evaluate((selector) => {
    const close = document.querySelector(selector) as HTMLElement | null;
    close?.click();
    const root = document.querySelector('.tutorials-root') as HTMLElement | null;
    const node = root ? grok.shell.dockManager.findNode(root) : null;
    if (node)
      grok.shell.dockManager.close(node);
  }, TUTORIAL_CLOSE);
  await expect.poll(() => page.evaluate(() => document.querySelector('.tutorials-root') == null &&
    document.querySelectorAll('.grok-tutorial-entry').length === 0), {message: 'the Tutorials panel closed'}).toBe(true);
}

const closesAtEnd = new WeakSet<Page>();

function closeAtFeatureEnd(page: Page): void {
  if (closesAtEnd.has(page))
    return;
  closesAtEnd.add(page);
  atFeatureEnd(page, async () => {
    closesAtEnd.delete(page);
    await closeTutorials(page);
  });
}

export const tutorialsOpen = Given('the Tutorials app is open', async (page: Page) => {
  await closeTutorials(page);
  closeAtFeatureEnd(page);
  await page.evaluate(() => grok.functions.call('Tutorials:trackOverview'));
  await expect(page.locator('.tutorials-root .tutorials-card').first(), 'the Tutorials panel with its cards').toBeVisible({timeout: pollMs(30000)});
}, {tier: 'api', description: 'the app\'s function docks the panel (no running tutorial left from before); the panel closes at feature end'});

export const tutorialsClosed = Given('the Tutorials app is closed', async (page: Page) => {
  await closeTutorials(page);
  closeAtFeatureEnd(page);
  await expect(page.locator('.tutorials-root'), 'the Tutorials panel').toHaveCount(0);
}, {tier: 'api', description: 'no Tutorials panel on the page, so the feature opens it the way a user does; closed again at feature end'});

export const startTutorial = When('user starts the {string} tutorial', async (page: Page, name: string) => {
  const card = page.locator(`.tutorials-root .tutorials-card[aria-label="${name.replace(/"/g, '\\"')}"]`);
  if (!await card.isVisible()) {
    await closeTutorials(page);
    closeAtFeatureEnd(page);
    await page.evaluate(() => grok.functions.call('Tutorials:trackOverview'));
  }
  await card.click({timeout: pollMs(30000)});
  await expect(page.locator('.tutorials-root-header h1'), 'the running tutorial\'s title').toHaveText(name, {timeout: pollMs(60000)});
  await expect(page.locator('.grok-tutorial-entry[aria-current="step"]').first(), 'the first step of the tutorial').toBeVisible({timeout: pollMs(60000)});
}, {tier: 'ui', description: 'clicks the tutorial\'s card in the Tutorials panel (opening the panel first if needed) and waits for its first step'});

/* A tutorial with a service prerequisite reads the service health the stand reports before it starts
   (Tutorial.checkService) and refuses to start unless the service is reported running, so a dev stack,
   which reports no health at all, refuses too. This gate is as strict as that check; the library's
   service gate lets such a stand go on. */
export const tutorialServiceReported = Given('the stand reports the {string} service the tutorial requires', async (page: Page, service: string) => {
  const gap = serviceGap(await reportedServices(page), service, true);
  test.skip(gap !== '', `the tutorial will not start: the stand does not report the ${service} service running (${gap})`);
}, {tier: 'api', description: 'a capability gate as strict as the tutorial\'s own prerequisite check: a stand that reports no service health skips too'});

export const tutorialNotCompleted = Given('the {string} tutorial is not completed yet', async (page: Page, name: string) => {
  await page.evaluate(([key, n]) => {
    if (grok.userSettings.getValue(key, n) != null)
      grok.userSettings.delete(key, n);
  }, [DATA_STORAGE_KEY, name]);
}, {tier: 'api', description: 'drops the tutorial\'s completion record in the page; pair with `the "tutorials" user settings are put back at feature end` first'});

export const tutorialProgress = Then('the tutorial progress should be {int} of {int}', async (page: Page, value: number, max: number) => {
  const read = () => page.evaluate(() => {
    const bar = document.querySelector('.tutorials-root-progress [role="progressbar"]');
    const caption = document.querySelector('.tutorials-root-progress > div')?.textContent ?? '';
    return bar ? `${bar.getAttribute('aria-valuenow')} of ${bar.getAttribute('aria-valuemax')} | ${caption}` : 'no progress bar';
  });
  await expect.poll(read, {message: 'the tutorial progress (bar | caption)'}).toBe(`${value} of ${max} | Step: ${value} of ${max}`);
}, {description: 'the progress bar\'s value and maximum and the "Step: N of M" caption agree on N of M'});

export const tutorialCompleted = Then('the {string} tutorial should be completed', async (page: Page, name: string) => {
  await expect(page.locator('.tutorials-root h3', {hasText: 'Congratulations!'}), 'the congratulations').toBeVisible();
  await expect.poll(() => page.evaluate(() => Array.from(document.querySelectorAll('.grok-tutorial-entry[role="checkbox"]'))
    .filter((e) => e.getAttribute('aria-checked') !== 'true').map((e) => e.getAttribute('aria-label'))), {message: 'steps not done'}).toEqual([]);
  await expect.poll(() => page.evaluate(([key, n]) => {
    const record = grok.userSettings.getValue(key, n);
    return record ? JSON.parse(record).isCompleted === true : false;
  }, [DATA_STORAGE_KEY, name]), {message: `the completion record of "${name}"`}).toBe(true);
}, {description: 'the congratulations are shown, every step listed is done, and the tutorial\'s record says completed'});

export const tutorialStepsListed = Then('the tutorial should have listed {int} steps', async (page: Page, count: number) => {
  await expect.poll(() => page.locator('.grok-tutorial-entry[role="checkbox"]').count(), {message: 'step entries listed'}).toBe(count);
}, {description: 'the number of step entries the running tutorial has shown so far'});

export const cardDone = Then('the {string} tutorial card should show it is done', async (page: Page, name: string) => {
  await expect(page.locator(`.tutorials-card[aria-label="${name.replace(/"/g, '\\"')}"]`), `the "${name}" card`)
    .toHaveAttribute('data-status', 'done', {timeout: pollMs(30000)});
}, {description: 'the card\'s own status, as the runner read the completion record'});

export const closeTutorial = When('user closes the tutorial', async (page: Page) => {
  await page.locator(TUTORIAL_CLOSE).click();
  await expect(page.locator('.grok-tutorial-entry'), 'the step entries of the closed tutorial').toHaveCount(0);
}, {tier: 'ui', description: 'the Close button in the running tutorial\'s header'});


/** A step entry by its instruction exactly as shown — a {string}, since instructions quote what the
 * learner types ('Name a column "BMI"'), which an element phrase cannot hold. A tutorial names a
 * modifier as the learner's keyboard labels it (`platformKeyMap`, src/tracks/shortcuts.ts), so on a Mac
 * the Ctrl and Alt a feature writes read Command and Option. */
async function stepEntry(page: Page, instruction: string) {
  const mac = await page.evaluate(() => navigator.platform.toLowerCase().includes('mac'));
  const shown = mac ? instruction.replace(/\bCtrl\b/g, 'Command').replace(/\bAlt\b/g, 'Option') : instruction;
  return page.locator(`.grok-tutorial-entry[role="checkbox"][aria-label="${shown.replace(/"/g, '\\"')}"]`);
}

export const stepDone = Then('the tutorial step {string} should be done', async (page: Page, instruction: string) => {
  const entries = await stepEntry(page, instruction);
  await expect.poll(async () => {
    const states = await entries.evaluateAll((els) => els.map((e) => `${e.getAttribute('aria-checked')}${e.getAttribute('aria-invalid') === 'true' ? ' invalid' : ''}`));
    return states.length === 0 ? 'not listed' : states.includes('true') ? 'done' : states.join(', ');
  }, {message: `the tutorial step "${instruction}"`}).toBe('done');
}, {description: 'the entry with exactly this instruction is listed and checked (aria-checked) — not shown as could-not-complete'});

export const stepDoneTimes = Then('the tutorial step {string} should be done {int} times', async (page: Page, instruction: string, times: number) => {
  const entries = await stepEntry(page, instruction);
  await expect.poll(() => entries.evaluateAll((els) => els.filter((e) => e.getAttribute('aria-checked') === 'true').length),
    {message: `checked entries "${instruction}"`}).toBe(times);
}, {description: 'for an instruction a tutorial repeats ("Open scatter plot"): that many of its entries are checked'});

export const stepNotDone = Then('the tutorial step {string} should not be done yet', async (page: Page, instruction: string) => {
  const entries = await stepEntry(page, instruction);
  await expect.poll(() => entries.evaluateAll((els) => els.length === 0 ? 'not listed' : els.every((e) => e.getAttribute('aria-checked') === 'false') ? 'pending' : 'done'),
    {message: `the tutorial step "${instruction}"`}).toBe('pending');
}, {description: 'the entry is listed and still unchecked — the claim that pairs with a gesture which must not tick it'});


/* The compute tutorials explain a view with a guided tour (ui-describer.ts): one popup per page, with
   "next" until the last, which has "done" (labelled "OK", "ok" or "clear" on a one-page tour). How many
   pages there are follows the view (one per viewer), so the step reads the tour, not a count. */
export const walkTour = When('user goes through the tour to its end', async (page: Page) => {
  const next = page.locator('[name="button-tour-next"]').filter({visible: true});
  const done = page.locator('[name="button-tour-done"]').filter({visible: true});
  await expect(next.or(done).first(), 'a tour popup').toBeVisible();
  for (let pages = 1; await next.count() > 0; pages++) {
    if (pages > 20)
      throw new Error('the tour still offers "next" after 20 pages');
    await next.first().click();
    await expect(next.or(done).first(), `tour page ${pages + 1}`).toBeVisible();
  }
  await done.first().click();
  await expect(next.or(done), 'the tour after "done"').toHaveCount(0);
}, {tier: 'ui', description: '"next" on every page of the guided tour, then "done" — the tour closes'});

export const stepListedTimes = Then('the tutorial step {string} should be listed {int} times', async (page: Page, instruction: string, times: number) => {
  const entries = await stepEntry(page, instruction);
  await expect.poll(() => entries.count(), {message: `entries "${instruction}"`}).toBe(times);
}, {description: 'for an instruction a tutorial repeats: its Nth entry is on the list — the tutorial has prepared that step'});

export const entityTypeExists = Then('the entity type {string} should exist', async (page: Page, type: string) => {
  await expect.poll(() => page.evaluate(async (type) => (await grok.dapi.entityTypes.list()).some((t: any) => t.name === type), type),
    {message: `the "${type}" entity type on the server`}).toBe(true);
}, {tier: 'api', description: 'the entity type is saved on the server'});

export const stickySchemaExists = Then('the Sticky Meta schema {string} should exist', async (page: Page, schema: string) => {
  await expect.poll(() => page.evaluate(async (schema) => (await grok.dapi.stickyMeta.getSchemas()).some((s: any) => s.name === schema), schema),
    {message: `the "${schema}" schema on the server`}).toBe(true);
}, {tier: 'api', description: 'the schema is saved on the server'});

export const cliffMoleculeIsCurrent = Then('{element} should show the current row', async (page: Page, target: ElementRef) => {
  const label = await (await locate(page, target)).getAttribute('aria-label') ?? '';
  const row = Number(/ of row (\d+)$/.exec(label)?.[1] ?? 0);
  if (row === 0)
    throw new Error(`${target.phrase} names no row: its label is "${label}"`);
  await expect.poll(() => page.evaluate(() => grok.shell.t?.currentRowIdx ?? -1), {message: `the current row (${target.phrase} is row ${row})`})
    .toBe(row - 1);
}, {description: 'an element named by its row ("molecule of row 12", numbered as the grid numbers rows) shows the table\'s current row'});
