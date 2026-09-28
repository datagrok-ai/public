/* The steps only Tutorials defines: starting a tutorial from its card, the running tutorial's
   progress and completion, the completion record it writes, and closing it. Everything a tutorial
   asks the learner to do is ordinary platform UI and uses the library's vocabulary; step entries,
   cards and the panel are named in elements.ts. */
import type {Page} from '@playwright/test';
import {Given, Then, When} from '@datagrok-libraries/bdd';
import {atFeatureEnd, expect, pollMs} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;
declare const DG: any;

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

/* The Sticky Meta tutorial saves an entity type and a schema under fixed names. They are removed at
   feature end and swept when the feature starts, since a killed run never reached its end; the schema
   goes first (a type a schema still uses is not deleted), and the server is read back after. The
   annotated value belongs to the schema and goes with it. */
async function removeStickyMetaFixtures(page: Page, schema: string, type: string): Promise<void> {
  await page.evaluate(async ([schema, type]) => {
    // the JS Schema has no id getter although deleteSchema takes the id: read it off the Dart entity
    const idOf = (s: any) => (window as any).grok_Entity_Get_Id(s.dart);
    for (const s of await grok.dapi.stickyMeta.getSchemas())
      if (s.name === schema)
        await grok.dapi.stickyMeta.deleteSchema(idOf(s));
    for (const t of await grok.dapi.entityTypes.list())
      if (t.name === type)
        await grok.dapi.entityTypes.delete(t);
  }, [schema, type]);
  await expect.poll(() => page.evaluate(async ([schema, type]) =>
    (await grok.dapi.stickyMeta.getSchemas()).filter((s: any) => s.name === schema).length +
    (await grok.dapi.entityTypes.list()).filter((t: any) => t.name === type).length, [schema, type]),
  {message: `the "${schema}" schema and the "${type}" entity type on the server`}).toBe(0);
}

export const stickyMetaFixturesGone = Given('the Sticky Meta schema {string} and entity type {string} are removed now and at feature end',
  async (page: Page, schema: string, type: string) => {
    await removeStickyMetaFixtures(page, schema, type);
    atFeatureEnd(page, () => removeStickyMetaFixtures(page, schema, type));
  }, {tier: 'api', description: 'swept before the walk and removed after it, each time read back from the server'});

export const entityTypeExists = Then('the entity type {string} should exist', async (page: Page, type: string) => {
  await expect.poll(() => page.evaluate(async (type) => (await grok.dapi.entityTypes.list()).some((t: any) => t.name === type), type),
    {message: `the "${type}" entity type on the server`}).toBe(true);
}, {tier: 'api', description: 'the entity type is saved on the server'});

export const stickySchemaExists = Then('the Sticky Meta schema {string} should exist', async (page: Page, schema: string) => {
  await expect.poll(() => page.evaluate(async (schema) => (await grok.dapi.stickyMeta.getSchemas()).some((s: any) => s.name === schema), schema),
    {message: `the "${schema}" schema on the server`}).toBe(true);
}, {tier: 'api', description: 'the schema is saved on the server'});

/* The Data Connectors tutorial saves a connection and a query under fixed names that learners share: on
   a stand where people took the tutorial by hand there are "Starbucks" connections of other users. Only
   the running user's own are removed — the query first, then its connection — now and at feature end. */
async function removeOwnConnection(page: Page, connection: string, query: string): Promise<void> {
  await page.evaluate(async ([connection, query]) => {
    const me = (await grok.dapi.users.current()).id;
    const mine = (e: any) => e.author?.id === me;
    for (const q of await grok.dapi.queries.list())
      if ((q.friendlyName === query || q.name === query) && mine(q))
        await grok.dapi.queries.delete(q);
    for (const c of await grok.dapi.connections.list())
      if ((c.friendlyName === connection || c.name === connection) && mine(c))
        await grok.dapi.connections.delete(c);
  }, [connection, query]);
  await expect.poll(() => page.evaluate(async ([connection, query]) => {
    const me = (await grok.dapi.users.current()).id;
    const mine = (e: any) => e.author?.id === me;
    return (await grok.dapi.queries.list()).filter((q: any) => (q.friendlyName === query || q.name === query) && mine(q)).length +
      (await grok.dapi.connections.list()).filter((c: any) => (c.friendlyName === connection || c.name === connection) && mine(c)).length;
  }, [connection, query]), {message: `the user's own "${connection}" connection and "${query}" query`, timeout: pollMs(30000)}).toBe(0);
}

export const ownConnectionGone = Given('the user\'s own connection {string} and query {string} are removed now and at feature end',
  async (page: Page, connection: string, query: string) => {
    await removeOwnConnection(page, connection, query);
    atFeatureEnd(page, () => removeOwnConnection(page, connection, query));
  }, {tier: 'api', description: 'only what the running user authored: other users\' entities of the same name stay'});

async function removeOwnProject(page: Page, project: string): Promise<void> {
  await page.evaluate(async (project) => {
    const me = (await grok.dapi.users.current()).id;
    for (const p of await grok.dapi.projects.list())
      if ((p.friendlyName === project || p.name === project) && p.author?.id === me)
        await grok.dapi.projects.delete(p);
  }, project);
  await expect.poll(() => page.evaluate(async (project) => {
    const me = (await grok.dapi.users.current()).id;
    return (await grok.dapi.projects.list()).filter((p: any) => (p.friendlyName === project || p.name === project) && p.author?.id === me).length;
  }, project), {message: `the user's own "${project}" project`, timeout: pollMs(30000)}).toBe(0);
}

export const ownProjectGone = Given('the user\'s own project {string} is removed now and at feature end', async (page: Page, project: string) => {
  await removeOwnProject(page, project);
  atFeatureEnd(page, () => removeOwnProject(page, project));
}, {tier: 'api', description: 'only what the running user authored: other users\' dashboards of the same name stay'});
