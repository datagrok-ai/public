import {readFileSync} from 'node:fs';
import type {Page} from '@playwright/test';
import {Given, Then, When} from '@datagrok-libraries/bdd';
import {type ElementRef, atFeatureEnd, el, expect, locate, silent} from '@datagrok-libraries/bdd/runtime';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';

declare const grok: any;

const SCHEMA = 'moltrack_ems';
const remembered = new WeakMap<Page, string>();

const BULK_EXAMPLE = new URL('../fixtures/bulk-registration.csv', import.meta.url);
const BULK_NAMES = ['Hentriacontanoic acid', 'Hexacosane-1,26-diamine', 'Triacontane-1,30-diol'];

async function openView(page: Page, view: string): Promise<void> {
  if (!await page.evaluate(() => typeof (window as any).grok?.shell?.closeAll === 'function').catch(() => false))
    await loggedIn(page);
  await page.evaluate(async (path) => {
    grok.shell.addView(await grok.functions.call('Moltrack:molTrackEmsApp', {path}));
  }, view);
  // loaded: the first row, or the empty state a table with no rows shows
  await page.locator('[data-u2="domain-list"] :is([data-u2-row], .u2-domain-list-empty)').filter({visible: true})
    .first().waitFor({timeout: 60000});
}

export const openEms = Given('user opens the MolTrack EMS app', (page: Page) => openView(page, ''),
  {tier: 'api', enters: 'MolTrack EMS',
    description: 'runs the app function and waits for the compound list to load (its first row, or its empty state); ' +
      'needs this package published'});

export const openEmsView = Given('user opens the {string} view of MolTrack EMS',
  (page: Page, view: string) => openView(page, view.replace(/ /g, '')),
  {tier: 'api', enters: 'MolTrack EMS',
    description: 'runs the app function on that view (Compounds, Batches, Additions, Assays, Assay Runs, Assay Results) ' +
      'and waits for its list to load'});

export const demoData = Given('the MolTrack EMS demo data is on the server', async (page: Page) => {
  silent(page);
  const missing = await page.evaluate(async (schema) => {
    const count = () => Promise.all(['compound', 'assay_result']
      .map((table) => grok.dapi.domains.table(`${schema}.${table}`).count() as Promise<number>));
    if ((await count()).every((n) => n > 0))
      return '';
    const loaded: string = await grok.functions.call('Moltrack:loadEmsDemoData');
    return (await count()).every((n) => n > 0) ? '' : loaded;
  }, SCHEMA);
  if (missing !== '')
    throw new Error(`the MolTrack EMS demo data did not load: ${missing}`);
}, {tier: 'api', description: 'loads the package\'s demo set (Moltrack:loadEmsDemoData: the compounds and batches of ' +
  'files/ems-demo, the additions, the hepatocyte stability assay, its runs and results) where the schema has no ' +
  'compounds or no results; a fixture loaded once and kept, as a dataset is'});

/** Soft-deletes the live compounds of that structure and their batches; throws if any is left. */
async function sweep(page: Page, structure: string): Promise<void> {
  const left = await page.evaluate(async ([schema, s]) => {
    const compounds = grok.dapi.domains.table(`${schema}.compound`);
    const filter = {property: 'structure', operator: '=', value: s};
    const ids: string[] = (await compounds.query({filter, columns: ['structure']})).map((r: any) => r.id);
    if (ids.length > 0) {
      const batches = grok.dapi.domains.table(`${schema}.compound_batch`);
      const owned = {property: 'compound_id', operator: '=', value: ids};
      if (await batches.count({filter: owned}) > 0)
        await batches.deleteWhere(owned);
      await compounds.deleteWhere(filter);
    }
    return compounds.count({filter});
  }, [SCHEMA, structure] as const);
  if (left !== 0)
    throw new Error(`${left} compound(s) with the structure ${structure} are still registered after the delete`);
}

export const noCompound = Given('no compound with the structure {string} is registered, now and at feature end',
  async (page: Page, structure: string) => {
    silent(page);
    await sweep(page, structure);
    atFeatureEnd(page, () => sweep(page, structure));
  }, {tier: 'api', description: 'deletes the live compounds of that standardized structure and their batches (EMS deletes are soft: ' +
    'the rows leave every live query) and proves none is left, now and again when the feature ends'});

export const rememberCompound = When('user remembers the compound the breadcrumbs name', async (page: Page) => {
  const current = (await locate(page, el('breadcrumbs'))).locator('[aria-current="page"]').first();
  await expect(current).toHaveText(/^DG-\d+$/);
  if (!remembered.has(page))
    atFeatureEnd(page, async () => {
      remembered.delete(page);
    });
  remembered.set(page, (await current.textContent())!.trim());
}, {tier: 'ui', description: 'the corporate ID (DG-…) the breadcrumbs show for the record the form holds'});

export const namesRemembered = Then('{element} should name the remembered compound', async (page: Page, target: ElementRef) => {
  const id = remembered.get(page);
  if (id === undefined)
    throw new Error('no compound was remembered: "user remembers the compound the breadcrumbs name" comes first');
  await expect((await locate(page, target)).filter({visible: true}).first()).toContainText(id);
}, {description: 'the element\'s text holds the corporate ID remembered from the breadcrumbs'});

export const compoundCount = Then('{int} compound(s) with the structure {string} should be registered on the server', async (page: Page, n: number, structure: string) => {
  await expect.poll(() => page.evaluate(([table, s]) =>
    grok.dapi.domains.table(table).count({filter: {property: 'structure', operator: '=', value: s}}),
    [`${SCHEMA}.compound`, structure] as const), {message: `live compounds with the structure ${structure}`}).toBe(n);
}, {tier: 'api', description: 'reads the compound table on the server: the live rows of that standardized structure'});

const rememberedTexts = new WeakMap<Page, Map<string, string>>();

export const rememberText = When('user remembers the text of {element}', async (page: Page, target: ElementRef) => {
  const shown = (await locate(page, target)).filter({visible: true}).first();
  await expect(shown).not.toHaveText('');
  let texts = rememberedTexts.get(page);
  if (texts === undefined) {
    rememberedTexts.set(page, texts = new Map<string, string>());
    atFeatureEnd(page, async () => {
      rememberedTexts.delete(page);
    });
  }
  texts.set(target.phrase, (await shown.textContent())!.trim());
}, {tier: 'ui', description: 'the visible text of a plain element (a status slot), for a later claim that it came back'});

function rememberedText(page: Page, target: ElementRef): string {
  const text = rememberedTexts.get(page)?.get(target.phrase);
  if (text === undefined)
    throw new Error(`the text of ${target.phrase} was not remembered: "user remembers the text of ${target.phrase}" comes first`);
  return text;
}

export const hasRememberedText = Then('{element} should have the remembered text', async (page: Page, target: ElementRef) => {
  await expect((await locate(page, target)).filter({visible: true}).first()).toHaveText(rememberedText(page, target));
}, {description: 'the visible text is the one remembered for the same phrase'});

export const lostRememberedText = Then('{element} should not have the remembered text',
  async (page: Page, target: ElementRef) => {
    await expect((await locate(page, target)).filter({visible: true}).first()).not.toHaveText(rememberedText(page, target));
  }, {description: 'the visible text has moved off the one remembered for the same phrase'});

export const batchCount = Then('the compound with the structure {string} should have {int} batch(es) on the server',
  async (page: Page, structure: string, n: number) => {
    await expect.poll(() => page.evaluate(async ([schema, s]) => {
      const ids = (await grok.dapi.domains.table(`${schema}.compound`).query(
        {filter: {property: 'structure', operator: '=', value: s}, columns: ['structure']})).map((r: any) => r.id);
      return ids.length === 0 ? 0 : grok.dapi.domains.table(`${schema}.compound_batch`)
        .count({filter: {property: 'compound_id', operator: '=', value: ids}});
    }, [SCHEMA, structure] as const), {message: `live batches of the compound ${structure}`}).toBe(n);
  }, {tier: 'api', description: 'reads the batch table on the server: the live batches of the compound of that standardized structure'});

/** Soft-deletes the batches the bulk example registers, and its new compounds by name. */
async function sweepBulk(page: Page): Promise<void> {
  await page.evaluate(async ([schema, names]) => {
    const batches = grok.dapi.domains.table(`${schema}.compound_batch`);
    const compounds = grok.dapi.domains.table(`${schema}.compound`);
    const ours = {property: 'epa_batch_id', operator: 'like', value: 'BDD-%'};
    if (await batches.count({filter: ours}) > 0)
      await batches.deleteWhere(ours);
    const named = {property: 'common_name', operator: '=', value: names};
    if (await compounds.count({filter: named}) > 0)
      await compounds.deleteWhere(named);
  }, [SCHEMA, BULK_NAMES] as const);
}

export const bulkExample = Given('the bulk registration example is open, and its compounds are not registered now and at feature end',
  async (page: Page) => {
    silent(page);
    await sweepBulk(page);
    atFeatureEnd(page, () => sweepBulk(page));
    await page.evaluate((csv) => {
      const df = (window as any).DG.DataFrame.fromCsv(csv);
      df.name = 'bulk-registration';
      grok.shell.addTableView(df);
    }, readFileSync(BULK_EXAMPLE, 'utf8'));
  }, {tier: 'api', description: 'opens bdd/fixtures/bulk-registration.csv (three new structures, a salt of one of them, a ' +
    'registered one, a broken one) as a table, and deletes what it registers now and when the feature ends'});

export const bulkRegistered = Then('the bulk registration example should have {int} new compounds and {int} batches on the server',
  async (page: Page, compounds: number, batches: number) => {
    await expect.poll(() => page.evaluate(([schema, names]) => Promise.all([
      grok.dapi.domains.table(`${schema}.compound`).count({filter: {property: 'common_name', operator: '=', value: names}}),
      grok.dapi.domains.table(`${schema}.compound_batch`).count(
        {filter: {property: 'epa_batch_id', operator: 'like', value: 'BDD-%'}}),
    ]), [SCHEMA, BULK_NAMES] as const), {message: 'compounds and batches of the bulk example'}).toEqual([compounds, batches]);
  }, {tier: 'api', description: 'reads the server: the example\'s new compounds by name, its batches by their BDD- ids'});
