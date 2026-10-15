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

declare const grok: any;
declare const DG: any;

const CAMPAIGNS = 'System:AppData/HitTriage/Hit Design/campaigns';

/** Removes every campaign of the feature's key (`<key>-<n>`) — its folder and its rows in the hitdesign
 * database — and returns what is still there: '' once all of it is gone. */
function removeCampaigns(page: Page, key: string): Promise<string> {
  return page.evaluate(async ([root, k]) => {
    const ids: string[] = (await grok.dapi.files.list(root))
      .filter((f: any) => f.isDirectory && f.name.startsWith(`${k}-`)).map((f: any) => f.name);
    for (const id of ids)
      await grok.dapi.files.delete(`${root}/${id}`);
    const conn = (await grok.dapi.connections.filter('name = "hitdesign"').list())
      .find((c: any) => c.nqName === 'HitTriage:hitdesign');
    const where = `app_name = 'Hit Design' and campaign_id like '${k}-%'`;
    for (const table of ['update_logs', 'campaign_locks', 'campaign_vids'])
      await conn.query('bdd-hit-design-cleanup', `delete from hitdesign.${table} where ${where} returning campaign_id`)
        .apply();
    const left: string[] = [];
    for (const id of ids)
      if (await grok.dapi.files.exists(`${root}/${id}`))
        left.push(`folder ${id}`);
    for (const table of ['update_logs', 'campaign_locks', 'campaign_vids']) {
      const rows = (await conn.query('bdd-hit-design-cleanup', `select campaign_id from hitdesign.${table} where ${where}`)
        .apply()).rowCount;
      if (rows > 0)
        left.push(`${rows} row(s) in hitdesign.${table}`);
    }
    return left.join(', ');
  }, [CAMPAIGNS, key] as [string, string]);
}

async function removedAndReadBack(page: Page, name: string, key: string): Promise<void> {
  let left = '';
  await expect.poll(async () => left = await removeCampaigns(page, key),
    {message: `the Hit Design campaign "${name}" removed from the share and the hitdesign database`, timeout: pollMs(30000)})
    .toBe('').catch((e) => {
      throw new Error(`the Hit Design campaign "${name}" is still on the server: ${left || e?.message}`);
    });
}

export const campaignOnServer = Given('the Hit Design campaign {string} with a table of:',
  async (page: Page, name: string, table: string[][]) => {
    const key = 'BDD' + name.replace(/[^A-Za-z]/g, '').toUpperCase();
    await removedAndReadBack(page, name, key);
    atFeatureEnd(page, () => removedAndReadBack(page, name, key));
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
        name: id, friendlyName: n, templateName: n, status: 'In Progress', createDate: '2026/10/15', campaignFields: {},
        columnSemTypes: semTypes, columnTypes: types, rowCount: rows.length - 1, filteredRowCount: rows.length - 1,
        savePath: `System.AppData/HitTriage/Hit Design/campaigns/${id}/enriched_table.csv`, version: 1,
        template: {name: n, key: k, campaignFields: [], stages: ['Design'],
          compute: {descriptors: {enabled: false, args: []}, functions: [], scripts: [], queries: []}}}));
    }, [CAMPAIGNS, name, key, table] as [string, string, string, string[][]]);
  }, {tier: 'api', description: 'a campaign of that friendly name made through the JS API, its id "BDD" and the name\'s letters ("BDDBDDCOLUMNS-1"): its table the rows given (the first one the header), its template embedded (no compute, one "Design" stage); removed at feature end with its rows in the hitdesign database, read back, and swept the same way before it is made'});

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
