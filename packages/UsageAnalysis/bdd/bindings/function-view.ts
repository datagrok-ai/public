/* The function view's own parts (core `function_view.dart`, `func_runs_panel.dart`): the result header with
   the run's number and actions, the results under the form, the docked Run history panel and its rows,
   and the switches of the analysis views (Compute2's fitting and sensitivity analysis), whose switch sits
   beside an input instead of carrying its label. The runs a feature makes are deleted with it. */
import type {Page} from '@playwright/test';
import {element, kind, Given, Then, When} from '@datagrok-libraries/bdd';
import {atFeatureEnd, expect, pollMs, serverRequests} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

element('run header', {selector: '.grok-function-view-run-header',
  description: 'the line above the results: star, the run\'s title or number, when and who, and the run\'s actions'});
element('run header star', {selector: '.grok-function-view-run-header .grok-favorite-star .d4-star',
  description: 'the star of the run shown; starring keeps the run with its tables'});
element('function results', {selector: '.grok-function-view-outputs',
  description: 'the results under the form, with the "Inputs changed" note while they are stale'});
element('run history', {selector: '.grok-func-runs',
  description: 'the Run history panel docked in the right column'});
element('run history search', {selector: '.grok-func-runs .d4-search-input',
  description: 'the search box at the top of the Run history panel'});

kind('run', {
  selector: '.grok-func-runs-row',
  match: ['label'],
  labelSelector: '.grok-func-runs-number, .grok-func-runs-title',
  parts: {
    star: '.d4-star',
    changes: '.grok-func-runs-tags',
    delete: '[name="icon-trash-alt"]',
    edit: '[name="icon-pencil"]',
  },
  description: 'a row of the Run history panel, by its title or "Run #N"; "delete of X run" is its trash icon',
});

kind('run filter', {
  selector: '.grok-func-runs-chip',
  match: ['text'],
  description: 'a filter chip of the Run history panel: Mine, Starred, Failed',
});

async function deleteRuns(page: Page, func: string): Promise<void> {
  await page.evaluate(async (name) => {
    const calls = grok.dapi.functions.calls.allPackageVersions();
    for (const call of await calls.filter(`func.name = "${name}"`).list({pageSize: 500}))
      await calls.delete(call);
  }, func);
  await expect.poll(() => page.evaluate(async (name) =>
    (await grok.dapi.functions.calls.allPackageVersions().filter(`func.name = "${name}"`).list()).length, func),
  {message: `runs of ${func} left on the server`, timeout: pollMs(30000)}).toBe(0);
}

/* Registered first, so the runs go before the script the fixture step deletes. */
export const runsDeleted = Given('the runs of the {string} function are deleted at the end', async (page: Page, func: string) => {
  atFeatureEnd(page, () => deleteRuns(page, func), true);
  await deleteRuns(page, func);
}, {tier: 'api', description: 'every run of the function, with the tables starring kept, at feature start and end'});

async function ownRunIds(page: Page, func: string): Promise<string[]> {
  return page.evaluate(async (name) => {
    const filter = `func.name = "${name}" and session.user.id = "${grok.shell.user.id}"`;
    return (await grok.dapi.functions.calls.allPackageVersions().filter(filter).list({pageSize: 5000})).map((c: any) => c.id);
  }, func);
}

/* For a shared function (a package's model): only the runs this user made after the step, by their ids. */
export const ownRunsDeleted = Given('the runs of the {string} function made from now on are deleted at the end',
  async (page: Page, func: string) => {
    const before = new Set(await ownRunIds(page, func));
    atFeatureEnd(page, async () => {
      const made = (await ownRunIds(page, func)).filter((id) => !before.has(id));
      await page.evaluate(async (ids) => {
        const calls = grok.dapi.functions.calls.allPackageVersions();
        for (const id of ids)
          await calls.delete(await calls.find(id));
      }, made);
    }, true);
  }, {tier: 'api', description: 'the runs of the function this user makes after the step, deleted by id at feature end'});

/* An analysis form puts a switch with each input instead of in its label; the switch turns the input into a
   range (sensitivity analysis), a fitted input or a target (fitting). */
async function setAnalysisSwitch(page: Page, caption: string, on: boolean): Promise<void> {
  const label = page.locator('label', {hasText: new RegExp(`^${caption.replace(/[.*+?^${}()|[\]\\]/g, '\\$&')}$`)})
    .filter({visible: true}).first();
  const root = label.locator('xpath=ancestor::*[contains(@class, "ui-input-root")][1]');
  // sensitivity analysis puts the switch inside the input, fitting right before it
  const inside = root.locator('.ui-input-switch');
  const toggle = (await inside.count()) > 0 ? inside.first() : root.locator('xpath=preceding-sibling::*[1]').locator('.ui-input-switch').first();
  await toggle.scrollIntoViewIfNeeded();
  if (await toggle.evaluate((e) => e.classList.contains('ui-input-switch-on')) === on)
    return;
  const clicked = await toggle.elementHandle();
  await toggle.click();
  // switching on rebuilds the input (a range shows "<caption> min" and "max"), so the switch may be replaced
  await expect.poll(() => clicked!.evaluate((e, want) => !e.isConnected || e.classList.contains('ui-input-switch-on') === want, on),
    {message: `the switch of "${caption}" turned ${on ? 'on' : 'off'}`}).toBe(true);
}

export const analysisSwitch = When('user switches on {string} in the analysis form',
  (page: Page, caption: string) => setAnalysisSwitch(page, caption, true),
  {description: 'the switch of the input captioned so in a fitting or sensitivity analysis view; one already on stays on'});

export const analysisSwitchOff = When('user switches off {string} in the analysis form',
  (page: Page, caption: string) => setAnalysisSwitch(page, caption, false),
  {description: 'the switch of the input captioned so in a fitting or sensitivity analysis view; one already off stays off'});

/* A space holds a run as a link (a project relation); the JS API lists a space's children without them. */
export const spaceLinksRun = Then('the space {string} on the server should link a run of the {string} function',
  async (page: Page, space: string, func: string) => {
    const api = await serverRequests(page);
    const id = await page.evaluate(async (n) =>
      (await grok.dapi.spaces.list()).find((p: any) => p.name === n || p.friendlyName === n)?.id ?? null, space);
    expect(id, `the space "${space}" on the server`).not.toBeNull();
    await expect.poll(async () => (await api.get<any[]>(`/projects/relations?projectId=${id}`))
      .filter((r) => r.entity?.['#type'] === 'FuncCall' && String(r.entity?.name ?? '').startsWith(func)).length,
    {message: `runs of ${func} the space "${space}" links`, timeout: pollMs(15000)}).toBeGreaterThan(0);
  }, {tier: 'api', description: 'a FuncCall relation of the space whose run belongs to the function'});
