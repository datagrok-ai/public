/* The steps only Tutorials defines: starting a tutorial from its card, the running tutorial's
   progress and completion, the completion record it writes, and closing it. Everything a tutorial
   asks the learner to do is ordinary platform UI and uses the library's vocabulary; step entries,
   cards and the panel are named in elements.ts. */
import type {Page} from '@playwright/test';
import {Given, Then, When} from '@datagrok-libraries/bdd';
import {atFeatureEnd, expect, pollMs} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

const DATA_STORAGE_KEY = 'tutorials';

/** Closes the running tutorial through its own Close button (which rejects the step it waits on and
 * closes its table) and the Tutorials panel with it: the panel is not a view, so the shell reset
 * leaves it — and a tutorial left running keeps listening to the platform's events and would take
 * the next feature's gestures as its steps. */
async function closeTutorials(page: Page): Promise<void> {
  await page.evaluate(() => {
    const close = document.querySelector('.tutorials-root-header button[aria-label="Close"], .tutorials-root-header .ui-btn') as HTMLElement | null;
    close?.click();
    const root = document.querySelector('.tutorials-root') as HTMLElement | null;
    const node = root ? grok.shell.dockManager.findNode(root) : null;
    if (node)
      grok.shell.dockManager.close(node);
  });
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
  await expect.poll(() => page.evaluate(() => [...document.querySelectorAll('.grok-tutorial-entry[role="checkbox"]')]
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
  await page.locator('.tutorials-root-header button').last().click();
  await expect(page.locator('.grok-tutorial-entry'), 'the step entries of the closed tutorial').toHaveCount(0);
}, {tier: 'ui', description: 'the Close button in the running tutorial\'s header'});


/** A step entry by its instruction exactly as shown — a {string}, since instructions quote what the
 * learner types ('Name a column "BMI"'), which an element phrase cannot hold. */
function stepEntry(page: Page, instruction: string) {
  return page.locator(`.grok-tutorial-entry[role="checkbox"][aria-label="${instruction.replace(/"/g, '\\"')}"]`);
}

export const stepDone = Then('the tutorial step {string} should be done', async (page: Page, instruction: string) => {
  const entries = stepEntry(page, instruction);
  await expect.poll(async () => {
    const states = await entries.evaluateAll((els) => els.map((e) => `${e.getAttribute('aria-checked')}${e.getAttribute('aria-invalid') === 'true' ? ' invalid' : ''}`));
    return states.length === 0 ? 'not listed' : states.includes('true') ? 'done' : states.join(', ');
  }, {message: `the tutorial step "${instruction}"`}).toBe('done');
}, {description: 'the entry with exactly this instruction is listed and checked (aria-checked) — not shown as could-not-complete'});

export const stepDoneTimes = Then('the tutorial step {string} should be done {int} times', async (page: Page, instruction: string, times: number) => {
  const entries = stepEntry(page, instruction);
  await expect.poll(() => entries.evaluateAll((els) => els.filter((e) => e.getAttribute('aria-checked') === 'true').length),
    {message: `checked entries "${instruction}"`}).toBe(times);
}, {description: 'for an instruction a tutorial repeats ("Open scatter plot"): that many of its entries are checked'});

export const stepNotDone = Then('the tutorial step {string} should not be done yet', async (page: Page, instruction: string) => {
  const entries = stepEntry(page, instruction);
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
  const entries = stepEntry(page, instruction);
  await expect.poll(() => entries.count(), {message: `entries "${instruction}"`}).toBe(times);
}, {description: 'for an instruction a tutorial repeats: its Nth entry is on the list — the tutorial has prepared that step'});
