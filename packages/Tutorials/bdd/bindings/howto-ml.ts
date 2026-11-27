/* What the "Machine learning and statistics" how-to clips (features/howto/10-ml-and-statistics) need besides the
   library's vocabulary: demo tables with many numeric columns and a category to color by, opened off camera like the
   library's datasets, and the sweep of what a Forge training leaves on the server (Forge's own bdd steps are not
   loaded in this project). */
import type {Page} from '@playwright/test';
import {dataset, Given, kind, Then, When} from '@datagrok-libraries/bdd';
import {atFeatureEnd, expect, pollMs} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

/* A panel of the Forge Train view's results by its title ("Performance" results panel), as Forge's bdd project names
   it (Forge/bdd/bindings/elements.ts), which this project does not load. */
kind('results panel', {
  selector: '.forge-results-panel, .forge-runs',
  match: ['label'],
  labelSelector: ':scope > .forge-results-title',
  parts: {title: ':scope > .forge-results-title'},
  description: 'a panel of the Forge Train view\'s results by its title',
});

/* A Forge collapsible group (the Apply dialog's Columns, the Train view's Data): aria-expanded on its header, as Forge's
   bdd project names it with this alias. */
kind('collapsible group', {
  selector: '.forge-group',
  match: ['label'],
  labelSelector: '.forge-group-header label',
  parts: {header: '.forge-group-header', summary: '.forge-group-summary', body: '.forge-group-body'},
  description: 'a Forge collapsible group by its caption',
});

/** A categorically colored column paints whole cells in the grid; for a clip that colors a plot by it, the grid
 * paints only the letters, so the column does not shout from the table. */
export const gridColorsText = Given('the grid paints the colors of {string} column on its text', async (page: Page,
  column: string) => {
  await page.evaluate((c) => {
    const gc = grok.shell.tv?.grid?.col(c);
    if (!gc)
      throw new Error(`no "${c}" column in the grid`);
    gc.isTextColorCoded = true;
  }, column);
}, {tier: 'api', description: 'the grid column\'s isTextColorCoded, as "Apply to: Text" of the Color Coding menu'});

/** The Train view publishes its state on its results (Forge ARCHITECTURE.md, "Automation surface"): ready is
 * aria-busy false without aria-invalid, and a shown training carries the id of its recorded run. */
export const trainingShown = Then('the Train view should show a finished training', async (page: Page) => {
  const results = page.locator('.forge-train-results').filter({visible: true}).first();
  await expect.poll(async () => results.evaluate((e) =>
    `${e.getAttribute('aria-busy')}|${e.getAttribute('aria-invalid') ?? 'false'}|${e.hasAttribute('data-forge-run-id')}`)
    .catch(() => 'no Train view'),
  {message: 'the Train view\'s results: busy|invalid|a shown run', timeout: pollMs(60000)}).toBe('false|false|true');
}, {tier: 'ui', description: 'the results area of the Train view ready (aria-busy false, not aria-invalid) and showing a ' +
  'recorded run (data-forge-run-id)'});

const shownRuns = new WeakMap<Page, string>();

async function shownRun(page: Page): Promise<string> {
  return await page.locator('.forge-train-results').filter({visible: true}).first()
    .getAttribute('data-forge-run-id') ?? '';
}

export const rememberShownTraining = When('user remembers the training the Train view shows', async (page: Page) => {
  shownRuns.set(page, await shownRun(page));
}, {tier: 'api', description: 'the run id on the Train view\'s results (data-forge-run-id)'});

export const newTrainingShown = Then('the Train view should show a new training', async (page: Page) => {
  const before = shownRuns.get(page) ?? '';
  const results = page.locator('.forge-train-results').filter({visible: true}).first();
  await expect.poll(async () => results.evaluate((e, old) => {
    const run = e.getAttribute('data-forge-run-id');
    return `${e.getAttribute('aria-busy')}|${e.getAttribute('aria-invalid') ?? 'false'}|${run !== null && run !== old}`;
  }, before).catch(() => 'no Train view'),
  {message: 'the Train view\'s results: busy|invalid|a run other than the remembered one', timeout: pollMs(60000)})
    .toBe('false|false|true');
}, {tier: 'ui', description: 'ready, and the shown run id differs from the remembered one'});

dataset('winequality', {path: 'System:DemoFiles/winequality.csv', aliases: ['wine quality'],
  description: '6497 wines: eleven lab measurements, a quality score and the type, red or white (the view is "winequality")'});

/** Deletes the Forge models of exactly these names the Forge way (`Forge:deleteModel`: the row, its model file, its
 * training copy, its applications), then the training runs recorded for these tables (a run outlives its model), and
 * reads both back gone. Runs are matched by the table's name, so a clip trains on a table of its own name. */
async function sweepForge(page: Page, models: string[], tables: string[]): Promise<void> {
  const left = (): Promise<string[]> => page.evaluate(async ({models, tables}) => {
    const named = async (table: string, column: string, value: string): Promise<number> =>
      await grok.dapi.domains.table(table).count({property: column, operator: '=', value});
    const found: string[] = [];
    for (const name of models)
      if (await named('forge.model', 'name', name) > 0)
        found.push(`model "${name}"`);
    for (const table of tables)
      if (await named('forge.training_run', 'dataset_name', table) > 0)
        found.push(`runs of "${table}"`);
    return found;
  }, {models, tables});
  await page.evaluate(async ({models, tables}) => {
    for (const name of models) {
      const ids: string[] = (await grok.dapi.domains.table('forge.model').query(
        {filter: {property: 'name', operator: '=', value: name}})).map((r: any) => String(r.id));
      for (const id of ids)
        await grok.functions.call('Forge:deleteModel', {model: id});
    }
    for (const table of tables) {
      const runs = grok.dapi.domains.table('forge.training_run');
      for (let more = true, pages = 0; more && pages < 20; pages++)
        more = (await runs.deleteWhere({property: 'dataset_name', operator: '=', value: table}, {limit: 1000})).hasMore;
    }
  }, {models, tables});
  await expect.poll(left, {message: 'Forge models and training runs left on the server', timeout: pollMs(30000)})
    .toEqual([]);
}

/** A model to apply, trained and saved through Forge's API on a demo file loaded in the page but not opened as a view,
 * so the clip starts with the model already in the catalog. Its training run is recorded under [table]: sweep that
 * table's runs with the step below. */
export const forgeModelFixture = Given('a Forge model {string} predicting {string} from the demo file {string} as the ' +
  'table {string} is on the server', async (page: Page, model: string, target: string, file: string, table: string) => {
  const id = await page.evaluate(async ({model, target, file, table}) => {
    const df = await grok.data.files.openTable(`System:DemoFiles/${file}`);
    df.name = table;
    const result = await grok.functions.call('Forge:trainModel', {table: df, target, name: model, showProgress: false});
    return String(result?.modelId ?? '');
  }, {model, target, file, table});
  expect(id, `the id of the saved Forge model "${model}"`).not.toBe('');
}, {tier: 'api', description: 'Forge:trainModel with the method Forge suggests and no data kept, on the demo file ' +
  'renamed to the table name; the saved model\'s id read back from the result'});

export const forgeSweep = Given('the Forge model {string} and the training runs of the table {string} are removed now ' +
  'and when the feature ends', async (page: Page, model: string, table: string) => {
  await sweepForge(page, [model], [table]);
  atFeatureEnd(page, () => sweepForge(page, [model], [table]));
}, {tier: 'api', description: 'Forge:deleteModel for every model of exactly that name, then every forge.training_run ' +
  'row of that table name deleted through grok.dapi.domains; both read back gone, now and at feature end'});
