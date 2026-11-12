/* Hit Design's campaigns, as the features use them. A campaign is two files on the AppData share
   (campaign.json and its table, enriched_table.csv), and every save also writes its row in the
   hitdesign database's update log (a molecule edit would register V-iDs as well, which these features
   never do). The campaign a feature uses is made for it through the JS API, the template embedded the
   way a saved campaign carries it, and removed at feature end: its folder, and its rows in the
   package's database, read back until they are gone. A run killed before its end is swept up when the
   next one starts. */
import {Page} from '@playwright/test';
import {Given, When} from '@datagrok-libraries/bdd';
import {type ElementRef, atFeatureEnd, expect, gestures, locate, pollMs} from '@datagrok-libraries/bdd/runtime';
import {keyOf, removedAndReadBack} from './campaign.js';

declare const grok: any;
declare const DG: any;

const CAMPAIGNS = 'System:AppData/HitTriage/Hit Design/campaigns';

export const campaignOnServer = Given('the Hit Design campaign {string} with a table of:',
  async (page: Page, name: string, table: string[][]) => {
    const key = keyOf(name);
    await removedAndReadBack(page, 'Hit Design', name, key);
    atFeatureEnd(page, () => removedAndReadBack(page, 'Hit Design', name, key));
    await page.evaluate(async ([root, n, k, rows]) => {
      const id = `${k}-1`;
      const csv = rows.map((r) => r.join(',')).join('\n') + '\n';
      const types: {[_: string]: string} = {};
      const semTypes: {[_: string]: string | null} = {};
      for (const c of DG.DataFrame.fromCsv(csv).columns.toList()) {
        types[c.name] = c.type;
        semTypes[c.name] = c.name === 'Molecule' ? 'Molecule' : null;
      }
      await grok.dapi.files.writeAsText(`${root}/${id}/enriched_table.csv`, csv);
      await grok.dapi.files.writeAsText(`${root}/${id}/campaign.json`, JSON.stringify({
        name: id, friendlyName: n, templateName: n, status: 'In Progress', createDate: '2999/12/31', campaignFields: {},
        columnSemTypes: semTypes, columnTypes: types, rowCount: rows.length - 1, filteredRowCount: rows.length - 1,
        savePath: `System.AppData/HitTriage/Hit Design/campaigns/${id}/enriched_table.csv`, version: 1,
        template: {name: n, key: k, campaignFields: [], stages: ['Design'],
          compute: {descriptors: {enabled: false, args: []}, functions: [], scripts: [], queries: []}}}));
    }, [CAMPAIGNS, name, key, table] as [string, string, string, string[][]]);
    await showFirstOnLanding(page);
  }, {tier: 'api', description: 'a campaign of that friendly name made through the JS API, its id "BDD" and the name\'s letters ("BDDBDDCOLUMNS-1"): its table the rows given (the first one the header), its template embedded (no compute, one "Design" stage), created 2999/12/31 so the landing sorted newest first lists it at the top; removed at feature end with its rows in the hitdesign database, read back, and swept the same way before it is made'});

/* The landing lists the campaigns the page read once (the package caches the list for the page) and keeps its
   columns, grouping and sorting in the browser: a campaign made for a feature is listed only once the page drops its
   copy, and is in view at the top of a table sorted by Created, newest first, with no grouping. The landing's
   settings are put back at feature end. */
const LANDING_SETTINGS = /^HDCampaign|^Accordion:Hit-Design-campaigns-group-/;

async function showFirstOnLanding(page: Page): Promise<void> {
  const saved: [string, string][] = await page.evaluate((re) => {
    const keys = Object.keys(localStorage).filter((k) => new RegExp(re).test(k));
    const out = keys.map((k) => [k, localStorage.getItem(k) ?? ''] as [string, string]);
    for (const k of keys)
      localStorage.removeItem(k);
    localStorage.setItem('HDCampaignsTableSorting', JSON.stringify({columnName: 'Created', ascending: false}));
    const cache = (window as any).hittriage?._package?.campaignsCache;
    if (cache)
      delete cache['Hit Design'];
    return out;
  }, LANDING_SETTINGS.source);
  atFeatureEnd(page, () => page.evaluate(([re, back]) => {
    for (const k of Object.keys(localStorage).filter((k) => new RegExp(re).test(k)))
      localStorage.removeItem(k);
    for (const [k, v] of back)
      localStorage.setItem(k, v);
    const cache = (window as any).hittriage?._package?.campaignsCache;
    if (cache)
      delete cache['Hit Design'];
  }, [LANDING_SETTINGS.source, saved] as [string, [string, string][]]));
}

/** The Tags input of the Column Properties dialog is a key and value per row, with an empty row only
 * when the column has no tags: a new tag goes into that row, else into the one the first row's "+"
 * inserts. */
export const addTag = When('user adds the {string} tag with value {string} to {element}',
  async (page: Page, tag: string, value: string, target: ElementRef) => {
    const root = (await locate(page, target)).first()
      .locator('xpath=ancestor-or-self::*[contains(concat(" ", normalize-space(@class), " "), " ui-input-root ")][1]');
    const rows = root.locator('tr');
    const emptyRow = () => rows.evaluateAll((trs) => trs.findIndex((tr) => {
      const inputs = Array.from(tr.querySelectorAll('input'));
      return inputs.length === 2 && inputs.every((i) => i.value === '');
    }));
    if (await emptyRow() < 0)
      await rows.first().locator('button').first().click();
    let row = -1;
    await expect.poll(async () => row = await emptyRow(), {message: `an empty row in ${target.phrase}`}).toBeGreaterThanOrEqual(0);
    const [keyInput, valueInput] = [rows.nth(row).locator('input').nth(0), rows.nth(row).locator('input').nth(1)];
    await keyInput.click();
    await gestures.typeVerified(keyInput, tag, `the key of a new row in ${target.phrase}`);
    await valueInput.click();
    await gestures.typeVerified(valueInput, value, `the value of the "${tag}" row in ${target.phrase}`);
  }, {tier: 'ui', description: 'types the tag and its value into an empty row of a Tags input (the Column Properties dialog\'s), a new one made with the first row\'s "+" when every row holds a tag'});
