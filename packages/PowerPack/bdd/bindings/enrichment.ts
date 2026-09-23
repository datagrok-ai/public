/* Column enrichment of DB Explorer (PowerPack db-explorer.ts): the Enrich pane a database-bound
   column shows in the context panel, its enrichment rows, and the editor dialog, which is the
   platform's visual query editor cut down to its Data and Join parts. An enrichment is a JSON file
   under System:AppData/PowerPack/enrichments/<connection>/<db>/<schema>/<table>/<column>/; the ones
   a feature makes are deleted when it ends, with the queries it saved. */
import {Page} from '@playwright/test';
import {Given, Then, When, element, kind} from '@datagrok-libraries/bdd';
import {type ElementRef, atFeatureEnd, expect, gestures, locate, pollMs} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

const ENRICHMENTS = 'System:AppData/PowerPack/enrichments/';

element('add enrichment button', {selector: '.power-pack-enrich-add',
  description: '"+ Add enrichment" under the enrichments of the Enrich pane'});
element('no enrichments note', {selector: '.power-pack-enrich-empty',
  description: '"No enrichments yet." — what the Enrich pane shows for a column with none'});
kind('enrichment', {
  selector: '.power-pack-enrichment-row',
  match: ['label'],
  labelSelector: 'a.ui-link',
  parts: {link: 'a.ui-link', 'edit icon': '.fa-pencil', 'delete icon': '.fa-times'},
  description: 'a row of the Enrich pane by the enrichment\'s name; its link applies it, the pencil edits it, the cross deletes it',
});

/** Enrichments whose name starts with the prefix, on every column of every connection. */
async function deleteEnrichments(page: Page, prefix: string): Promise<void> {
  await page.evaluate(async ([root, p]) => {
    if (!await grok.dapi.files.exists(root))
      return;
    for (const f of await grok.dapi.files.list(root, true))
      if (!f.isDirectory && f.name.startsWith(p) && f.name.endsWith('.json'))
        await grok.dapi.files.delete(f.fullPath);
  }, [ENRICHMENTS, prefix] as [string, string]);
}

export const noEnrichments = Given('no enrichment whose name starts with {string} is on the server', async (page: Page, prefix: string) => {
  await deleteEnrichments(page, prefix);
  atFeatureEnd(page, () => deleteEnrichments(page, prefix));
}, {tier: 'api', description: 'deletes every enrichment whose name starts with the text, now and when the feature ends'});

export const enrichmentsOnServer = Then('{int} enrichment(s) named {string} should be on the server', async (page: Page, count: number, name: string) => {
  await expect.poll(() => page.evaluate(async ([root, n]) => (await grok.dapi.files.list(root, true))
    .filter((f: any) => f.name === `${n}.json`).length, [ENRICHMENTS, name] as [string, string]),
  {message: `the enrichment files named ${name}.json`}).toBe(count);
}, {tier: 'api', description: 'the enrichment files the server keeps under System:AppData/PowerPack/enrichments'});

export const savedQuery = Given('a query {string} on {string} reads {string}', async (page: Page, name: string, connection: string, sql: string) => {
  const id: string = await page.evaluate(async ([n, c, s]) => {
    for (const q of await grok.dapi.queries.filter(`friendlyName = "${n}"`).list())
      await grok.dapi.queries.delete(q);
    const conn = (await grok.dapi.connections.list({pageSize: 1000})).find((x: any) => x.nqName === c);
    if (!conn)
      throw new Error(`no connection ${c}`);
    const q = conn.query(n, s);
    await grok.dapi.queries.save(q);
    return String(q.id);
  }, [name, connection, sql] as [string, string, string]);
  atFeatureEnd(page, () => page.evaluate(async (i) => {
    const q = await grok.dapi.queries.find(i).catch(() => null);
    if (q)
      await grok.dapi.queries.delete(q);
  }, id));
}, {tier: 'api', description: 'a SQL query saved on the connection (an earlier one of that name replaced), deleted when the feature ends'});

/** The Join part of the dialog: its first row names the two tables, its second the key pair. */
function joinRows(page: Page, dialog: ElementRef) {
  return locate(page, dialog).then((d) => d.first().locator('.grok-join-item'));
}

export const openJoinedColumns = When('user opens the columns of the joined {string} table in {element}', async (page: Page, table: string, dialog: ElementRef) => {
  const tag = (await joinRows(page, dialog)).locator('.grok-join-content .d4-tag > div').filter({hasText: new RegExp(`\\.${table}\\b`)}).last();
  await tag.click();
}, {tier: 'ui', description: 'a click on the joined table\'s tag, which opens the "Select columns..." dialog of its columns'});

export const joinedTableReads = Then('the joined {string} table in {element} should read {string}', async (page: Page, table: string, dialog: ElementRef, text: string) => {
  const tag = (await joinRows(page, dialog)).locator('.grok-join-content .d4-tag > div').filter({hasText: new RegExp(`\\.${table}\\b`)}).last();
  await expect(tag, `the tag of the joined ${table} table`).toHaveText(text);
}, {description: 'the tag the Join row shows for the table: "datagrok.public.users_sessions(4/12)" — selected of all columns'});

export const pickJoinKey = When('user picks {string} as the key of the main table in {element}', async (page: Page, column: string, dialog: ElementRef) => {
  const selector = (await joinRows(page, dialog)).locator('.grok-join-row').nth(1).locator('.d4-column-selector').first();
  await gestures.openColumnSelector(page, selector);
  await gestures.pickInColumnGrid(page, column, 'the key of the main table', selector);
}, {tier: 'ui', description: 'the left column selector of the Join\'s "on" row'});

export const joinReads = Then('the join key in {element} should read {string}', async (page: Page, dialog: ElementRef, text: string) => {
  const on = (await joinRows(page, dialog)).locator('.grok-join-row').nth(1).locator('.d4-tag').first();
  await expect.poll(async () => (await on.innerText()).replace(/\s+/g, ' ').trim(), {message: 'the "on" row of the Join'}).toBe(text);
}, {description: 'the key pair of the join as shown: "session_id = id" — the main table\'s column on the left'});

export const tableColumnsExactly = Then('the table should have the columns {string}', async (page: Page, list: string) => {
  const wanted = list.split(',').map((n) => n.trim()).filter(Boolean);
  await expect.poll(() => page.evaluate(() => (grok.shell.t?.columns.names() ?? ['no current table']) as string[]), {message: 'the columns of the current table',
    timeout: pollMs(30000)}).toEqual(wanted);
}, {description: 'exactly these columns, in this order'});

export const scrollToMiddle = When('user scrolls {element} to the middle of its list', async (page: Page, target: ElementRef) => {
  const loc = (await locate(page, target)).first();
  await loc.evaluate((e) => e.scrollIntoView({block: 'center'}));
  await expect.poll(async () => {
    const box = await loc.boundingBox();
    const size = page.viewportSize();
    return !!box && !!size && box.y > 100 && box.y + box.height < size.height - 100;
  }, {message: `${target.phrase} away from the edges of the window`}).toBe(true);
}, {tier: 'ui', description: 'the list scrolled until the element sits mid-window, clear of the status bar a node on the last line hides under'});

element('layouts pane', {selector: '.d4-toolbox .d4-pane-layouts',
  description: 'the Layouts section of the toolbox: its Save button and the cards of the layouts that fit the table'});
kind('layout card', {
  selector: '.d4-pane-layouts .grok-suggestions-chart-card',
  match: ['label'],
  labelSelector: '.grok-gallery-grid-item-title',
  description: 'a saved layout in the Layouts section of the toolbox, by its name (the view it was saved from); a click applies it',
});

/** Layouts saved under a name during the feature are deleted when it ends: the Save of the Layouts
 * section names a layout after the view it was saved from. */
export const layoutsDeleted = Given('the layouts named {string} are deleted when the feature ends', async (page: Page, name: string) => {
  atFeatureEnd(page, () => page.evaluate(async (n) => {
    for (const l of await grok.dapi.layouts.list({pageSize: 5000}))
      if ([l.friendlyName, l.name].some((x: any) => String(x).toLowerCase() === n.toLowerCase()))
        await grok.dapi.layouts.delete(l);
  }, name));
}, {tier: 'api', description: 'every layout of that name (as the Layouts section shows it) is deleted at feature end'});

export const enrichmentJoins = Then('the enrichment {string} on the server should select the column {string}', async (page: Page, name: string, column: string) => {
  await expect.poll(() => page.evaluate(async ([root, n]) => {
    const files = (await grok.dapi.files.list(root, true)).filter((f: any) => f.name === `${n}.json`);
    if (files.length !== 1)
      return [`${files.length} enrichment files named ${n}.json`];
    return (JSON.parse(await grok.dapi.files.readAsText(files[0].fullPath)).fields ?? []).map((x: string) => String(x));
  }, [ENRICHMENTS, name] as [string, string]), {message: `the fields the enrichment ${name} selects, as saved`})
    .toEqual(expect.arrayContaining([expect.stringMatching(new RegExp(`(^|\.)${column}$`))]));
}, {tier: 'api', description: 'the saved configuration, read back: what the next application of the enrichment will join'});
