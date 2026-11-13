/* Hit Design's campaigns, as the features use them. A campaign saved on the server cannot be taken back whole: a
   molecule written into it registers its V-iD, and every save logs the campaign's update, in the hitdesign database,
   which the package's queries add to and never delete from. So a feature works on a fixed fixture made once per stand
   and reused (the bdd library's rule for what cannot be undone): the campaign named here, of a template of the same
   name with no compute and no campaign fields, written through the JS API. The campaign's files are written as it
   starts when the feature starts, since a killed run never reached its end, and again at feature end, read back. A
   molecule registered in it again is the same V-iD (the database keys its rows by campaign and molecule), so a run adds
   no row there. */
import {Page} from '@playwright/test';
import {Given, Then} from '@datagrok-libraries/bdd';
import {atFeatureEnd, expect, pollMs} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

const ROOT = 'System:AppData/HitTriage/Hit Design';

interface CampaignFile {
  path: string;
  text: string;
}

/** The campaign of that name (its friendly name) as the feature starts on it, its template made when missing; its files. */
function campaignFixture(page: Page, name: string): Promise<CampaignFile[]> {
  return page.evaluate(async ([root, n]) => {
    const key = 'BDD' + n.replace(/[^A-Za-z]/g, '').toUpperCase().slice(0, 6);
    const template = `${root}/templates/${n}.json`;
    if (!await grok.dapi.files.exists(template)) {
      await grok.dapi.files.writeAsText(template, JSON.stringify({name: n, key, campaignFields: [], stages: ['Design'],
        compute: {descriptors: {enabled: false, args: []}, functions: []}}));
    }
    const id = `${key}-1`;
    const files = [
      {path: `${root}/campaigns/${id}/enriched_table.csv`, text: 'Molecule,Stage,V-iD\nCCO,Design,\n'},
      {path: `${root}/campaigns/${id}/campaign.json`, text: JSON.stringify({
        name: id, friendlyName: n, templateName: n, status: 'In Progress', createDate: '2026/10/07', campaignFields: {},
        columnSemTypes: {'Molecule': 'Molecule', 'Stage': null, 'V-iD': null},
        columnTypes: {'Molecule': 'string', 'Stage': 'string', 'V-iD': 'string'},
        rowCount: 1, filteredRowCount: 1, savePath: `System.AppData/HitTriage/Hit Design/campaigns/${id}/enriched_table.csv`,
        version: 1})},
    ];
    for (const f of files)
      await grok.dapi.files.writeAsText(f.path, f.text);
    return files;
  }, [ROOT, name] as [string, string]);
}

export const campaignOnServer = Given('the Hit Design campaign {string} is on the server, and is as it was again when the feature ends',
  async (page: Page, name: string) => {
    await campaignFixture(page, name);
    atFeatureEnd(page, async () => {
      const files = await campaignFixture(page, name);
      await expect.poll(() => page.evaluate(async (fs) => {
        for (const f of fs) {
          if (await grok.dapi.files.readAsText(f.path) !== f.text)
            return f.path;
        }
        return '';
      }, files), {message: `the files of the Hit Design campaign "${name}", put back`, timeout: pollMs(30000)}).toBe('');
    });
  }, {tier: 'api', description: 'a fixture made once per stand (a campaign cannot be deleted from the database): the campaign of that friendly name and its template of the same name (no compute, no fields), written through the JS API; its files written as it starts, at the start, since a killed run never reached its end, and again at feature end, read back'});

export const campaignSaved = Then('the Hit Design campaign {string} should be saved with {int} rows, row {int} the molecule {string} with its V-iD',
  async (page: Page, name: string, rows: number, row: number, smiles: string) => {
    let seen = '';
    await expect.poll(async () => {
      const r = await page.evaluate(async ([root, n, i, want]) => {
        const campaigns = `${root}/campaigns`;
        for (const folder of await grok.dapi.files.list(campaigns)) {
          const json = `${campaigns}/${folder.name}/campaign.json`;
          if (!folder.isDirectory || !await grok.dapi.files.exists(json))
            continue;
          const c = JSON.parse(await grok.dapi.files.readAsText(json));
          if (c.friendlyName !== n)
            continue;
          const df = await grok.dapi.files.readCsv(c.savePath.replace(/^System\.AppData\//, 'System:AppData/'));
          const rdkit = await grok.functions.call('Chem:getRdKitModule');
          const canonical = (v: string) => {
            const mol = rdkit.get_mol(v);
            try {
              return mol.get_smiles();
            }
            finally {
              mol.delete();
            }
          };
          // the save the edit started may not have come yet: the table holds no such row then
          const has = i <= df.rowCount;
          const molecule = has ? String(df.col('Molecule')?.get(i - 1) ?? '') : '';
          return {rowCount: c.rowCount, rows: df.rowCount, molecule: molecule ? canonical(molecule) : '', want: canonical(want),
            vid: has ? String(df.col('V-iD')?.get(i - 1) ?? '') : ''};
        }
        return null;
      }, [ROOT, name, row, smiles] as [string, string, number, string]);
      seen = JSON.stringify(r);
      return r !== null && r.rowCount === rows && r.rows === rows && r.molecule === r.want && /^V\d+$/.test(r.vid);
    }, {message: `the campaign "${name}" as saved on the server`, timeout: pollMs(30000)}).toBe(true).catch(() => {
      throw new Error(`the campaign "${name}" on the server is not ${rows} rows with ${smiles} and a V-iD in row ${row}: ${seen}`);
    });
  }, {description: 'its files on the server, read back until the save the edit started has written both: campaign.json\'s row count, and the table\'s row read by RDKit with a V-iD (V and digits)'});
