/* The campaign view of Hit Design and PeptiHit, as its features use it. A campaign is two files on the AppData share
   (campaign.json and its table) and rows in the hitdesign database: every save logs the campaign's update, a save
   holds its lock, and a molecule written into the table registers its V-iD (a dictionary row, and a row tying it to
   the campaign). A feature's campaign is made through the JS API, its template embedded the way a saved campaign
   carries it (no template file), and removed at feature end and swept before it is made: its folder, its rows, and
   the dictionary rows it registered that nothing else uses. Which dictionary rows are its own is told by the highest
   dictionary id when the campaign was made, kept in the campaign's folder for a sweep after a killed run.

   A campaign is opened by its link (`?campaignId=`), not from the landing's list: the package lists the campaigns
   once per page, so a list read by an earlier feature on the same page does not hold a campaign made since. */
import {resolve} from 'path';
import {Page} from '@playwright/test';
import {Given, Then, When} from '@datagrok-libraries/bdd';
import {type ElementRef, atFeatureEnd, expect, gestures, locate, pollMs} from '@datagrok-libraries/bdd/runtime';
import {sharingLogin, sharingShownName} from '@datagrok-libraries/bdd/bindings/platform/steps';

declare const grok: any;
declare const DG: any;

export type App = 'Hit Design' | 'PeptiHit';
const ROOTS: Record<App, string> = {
  'Hit Design': 'System:AppData/HitTriage/Hit Design',
  'PeptiHit': 'System:AppData/HitTriage/PeptiHit',
};
const APP_PATHS: Record<App, string> = {'Hit Design': 'HitDesign', 'PeptiHit': 'PeptiHit'};
const FLOOR_FILE = 'bdd-vid-floor.txt';
const LAYOUT = 'System:AppData/HitTriage/Hit Design/HD-template-layout.layout';

interface Campaign {app: App; name: string; key: string; id: string}
const campaigns = new WeakMap<Page, Map<string, Campaign>>();

export const keyOf = (name: string): string => 'BDD' + name.replace(/[^A-Za-z]/g, '').toUpperCase();

export function campaignOf(page: Page, name: string): Campaign {
  const c = campaigns.get(page)?.get(name);
  if (!c)
    throw new Error(`no campaign "${name}" was made by this feature`);
  return c;
}

/** Removes every campaign of the key in the app — its folder, its rows in the hitdesign database, and the dictionary
 * rows above the floor its folder kept that no campaign refers to any more — and returns what is still there. */
function removeCampaigns(page: Page, app: App, key: string): Promise<string> {
  return page.evaluate(async ([root, a, k, floorFile]) => {
    const folders: string[] = (await grok.dapi.files.list(`${root}/campaigns`))
      .filter((f: any) => f.isDirectory && f.name.startsWith(`${k}-`)).map((f: any) => f.name);
    let floor: number | null = null;
    for (const id of folders) {
      const path = `${root}/campaigns/${id}/${floorFile}`;
      if (await grok.dapi.files.exists(path)) {
        const n = parseInt(await grok.dapi.files.readAsText(path), 10);
        if (Number.isFinite(n))
          floor = floor === null ? n : Math.min(floor, n);
      }
    }
    const conn = (await grok.dapi.connections.filter('name = "hitdesign"').list())
      .find((c: any) => c.nqName === 'HitTriage:hitdesign');
    const where = `app_name = '${a}' and campaign_id like '${k}-%'`;
    const vids: string[] = ((await conn.query('bdd-campaign-cleanup', `select vid from hitdesign.campaign_vids where ${where}`)
      .apply()).col('vid')?.toList() ?? []).filter((v: string) => /^V\d+$/.test(v));
    for (const table of ['update_logs', 'campaign_locks', 'campaign_vids'])
      await conn.query('bdd-campaign-cleanup', `delete from hitdesign.${table} where ${where} returning campaign_id`).apply();
    const own = floor !== null && vids.length > 0 ?
      `vid in (${vids.map((v) => `'${v}'`).join(', ')}) and id > ${floor} and ` +
        'not exists (select 1 from hitdesign.campaign_vids c where c.vid = vid_dictionary.vid)' : null;
    if (own)
      await conn.query('bdd-campaign-cleanup', `delete from hitdesign.vid_dictionary where ${own} returning vid`).apply();
    for (const id of folders)
      await grok.dapi.files.delete(`${root}/campaigns/${id}`);
    const left: string[] = [];
    for (const id of folders) {
      if (await grok.dapi.files.exists(`${root}/campaigns/${id}`))
        left.push(`folder ${id}`);
    }
    for (const table of ['update_logs', 'campaign_locks', 'campaign_vids']) {
      const rows = (await conn.query('bdd-campaign-cleanup', `select campaign_id from hitdesign.${table} where ${where}`)
        .apply()).rowCount;
      if (rows > 0)
        left.push(`${rows} row(s) in hitdesign.${table}`);
    }
    if (own) {
      const rows = (await conn.query('bdd-campaign-cleanup', `select vid from hitdesign.vid_dictionary where ${own}`).apply()).rowCount;
      if (rows > 0)
        left.push(`${rows} row(s) in hitdesign.vid_dictionary`);
    }
    return left.join(', ');
  }, [ROOTS[app], app, key, FLOOR_FILE] as [string, App, string, string]);
}

export async function removedAndReadBack(page: Page, app: App, name: string, key: string): Promise<void> {
  // a publish of the package re-creates its connection: a feature running meanwhile waits for it
  await expect.poll(() => page.evaluate(async () => (await grok.dapi.connections.filter('name = "hitdesign"').list())
    .some((c: any) => c.nqName === 'HitTriage:hitdesign')), {message: 'the package\'s hitdesign connection', timeout: pollMs(120000)})
    .toBe(true);
  let left = '';
  await expect.poll(async () => left = await removeCampaigns(page, app, key),
    {message: `the ${app} campaign "${name}" removed from the share and the hitdesign database`, timeout: pollMs(30000)})
    .toBe('').catch((e) => {
      throw new Error(`the ${app} campaign "${name}" is still on the server: ${left || e?.message}`);
    });
}

export interface CampaignOptions {
  /** the second account is the campaign's author and its only editor; everyone may view it */
  sharingUserEdits?: boolean;
  /** a Molecule field after Chemist, set to the SMILES */
  molecule?: {field: string, smiles: string};
}

/** A campaign of that friendly name, its id the key and "-1": its table the rows given (the first one the header), its
 * template embedded with the stages given, a "Chemist" field (and the molecule field asked for), the JS-only Chemical
 * Properties compute (MW), the Demo File Submit preset and, for Hit Design, the package's shipped layout file. */
export async function makeCampaign(page: Page, app: App, name: string, stages: string[], table: string[][],
  options: CampaignOptions = {}): Promise<Campaign> {
  const key = keyOf(name);
  await removedAndReadBack(page, app, name, key);
  atFeatureEnd(page, () => removedAndReadBack(page, app, name, key));
  const sharing = options.sharingUserEdits ? sharingLogin() : '';
  await page.evaluate(async ([root, a, n, k, rows, st, layout, floorFile, second, [molField, smiles]]) => {
    const id = `${k}-1`;
    const csv = rows.map((r) => r.map((v) => /[",\n]/.test(v) ? `"${v.replace(/"/g, '""')}"` : v).join(',')).join('\n') + '\n';
    const types: {[_: string]: string} = {};
    const semTypes: {[_: string]: string | null} = {};
    const df = DG.DataFrame.fromCsv(csv);
    for (const c of df.columns.toList()) {
      types[c.name] = c.type;
      semTypes[c.name] = c.name === 'Molecule' ? 'Molecule' : c.name === 'Helm' ? 'Macromolecule' : null;
    }
    const conn = (await grok.dapi.connections.filter('name = "hitdesign"').list())
      .find((c: any) => c.nqName === 'HitTriage:hitdesign');
    const floor = (await conn.query('bdd-campaign-floor', 'select coalesce(max(id), 0) as floor from hitdesign.vid_dictionary')
      .apply()).get('floor', 0);
    const allUsers = DG.Group.defaultGroupsIds['All users'];
    let owner: {[_: string]: any} = {};
    if (second) {
      const user = await grok.dapi.users.filter(`login = ${JSON.stringify(second)}`).first();
      if (!user)
        throw new Error(`no user "${second}" on the server`);
      owner = {authorUserId: user.id, authorUserFriendlyName: user.friendlyName,
        permissions: {view: [allUsers], edit: [user.group.id]}};
    }
    const folder = `${root}/campaigns/${id}`;
    await grok.dapi.files.writeAsText(`${folder}/${floorFile}`, String(floor));
    await grok.dapi.files.writeAsText(`${folder}/enriched_table.csv`, csv);
    await grok.dapi.files.writeAsText(`${folder}/campaign.json`, JSON.stringify({
      name: id, friendlyName: n, templateName: n, status: 'In Progress', createDate: '2026/10/15',
      campaignFields: {Chemist: 'BDD chemist', ...molField ? {[molField]: smiles} : {}}, columnSemTypes: semTypes, columnTypes: types, rowCount: df.rowCount,
      filteredRowCount: df.rowCount, savePath: `${root.replace('System:AppData', 'System.AppData')}/campaigns/${id}/enriched_table.csv`,
      version: 1, ...owner,
      template: {name: n, key: k, campaignFields: [{name: 'Chemist', type: 'String', required: false},
        ...molField ? [{name: molField, type: 'Molecule', required: false}] : []], stages: st,
        submit: {fName: 'demoFileSubmit', package: 'HitTriage'}, ...(a === 'Hit Design' ? {localLayoutPath: layout} : {}),
        compute: {descriptors: {enabled: false, args: []}, scripts: [], queries: [],
          functions: [{name: 'addChemPropertiesColumns', package: 'Chem', args: {MW: true, HBA: false, HBD: false, logP: false,
            logS: false, PSA: false, rotatableBonds: false, stereoCenters: false, moleculeCharge: false, molecularFormula: false}}]}},
    }));
  }, [ROOTS[app], app, name, key, table, stages, LAYOUT, FLOOR_FILE, sharing,
    [options.molecule?.field ?? '', options.molecule?.smiles ?? '']] as
    [string, App, string, string, string[][], string[], string, string, string, [string, string]]);
  const campaign = {app, name, key, id: `${key}-1`};
  if (!campaigns.has(page))
    campaigns.set(page, new Map());
  campaigns.get(page)!.set(name, campaign);
  return campaign;
}

const stagesOf = (list: string): string[] => list.split(',').map((s) => s.trim()).filter((s) => s.length > 0);

export const hitDesignCampaign = Given('the Hit Design campaign {string} with the stages {string} and a table of:',
  async (page: Page, name: string, stages: string, table: string[][]) => {
    await makeCampaign(page, 'Hit Design', name, stagesOf(stages), table);
  }, {tier: 'api', description: 'made through the JS API, its id "BDD" + the name\'s letters + "-1": the rows given (the first one the header), the template embedded (the stages given, a Chemist field set to "BDD chemist", the Chemical Properties compute with MW, the Demo File Submit preset, the package\'s HD-template-layout.layout); removed at feature end with its hitdesign rows and the V-iDs it registered, read back, and swept before it is made'});

export const hitDesignCampaignWithMolecule = Given(
  'the Hit Design campaign {string} with the stages {string}, the molecule field {string} set to {string} and a table of:',
  async (page: Page, name: string, stages: string, field: string, smiles: string, table: string[][]) => {
    await makeCampaign(page, 'Hit Design', name, stagesOf(stages), table, {molecule: {field, smiles}});
  }, {tier: 'api', description: 'as "the Hit Design campaign … with the stages …", its template also holding a Molecule field of that name after Chemist, set to the SMILES'});

export const viewOnlyCampaign = Given('the Hit Design campaign {string}, which only the sharing user may edit, with a table of:',
  async (page: Page, name: string, table: string[][]) => {
    await makeCampaign(page, 'Hit Design', name, ['Idea', 'Design'], table, {sharingUserEdits: true});
  }, {tier: 'api', description: 'as "the Hit Design campaign … with the stages …", the stages Idea and Design; its author is the second account (DATAGROK_SHARING_LOGIN, bddsecond), whose personal group alone may edit it, and All users may view it'});

/* The card layout a closed tracker keeps with the campaign names the table it was made over — "Table" while the
   campaign was the only one open. */
export const cardLayoutOver = Given('the card layout of the campaign {string} was made over a table named {string}',
  async (page: Page, name: string, table: string) => {
    const c = campaignOf(page, name);
    await page.evaluate(async ([root, id, t]) => {
      const path = `${root}/campaigns/${id}/campaign.json`;
      const campaign = JSON.parse(await grok.dapi.files.readAsText(path));
      const field = (column: string, inputType: string, top: number) => ({left: 83, top, width: 100, height: 20,
        type: 'field', viewerSettings: {column, format: null, inputType}});
      campaign.tilesViewerFormSketch = JSON.stringify({'#type': 'SketchState', table: t, formDesigned: false,
        elementStates: [field('V-iD', 'Text', 20), field('Stage', 'Text', 44), field('MW', 'Float', 68)]});
      await grok.dapi.files.writeAsText(path, JSON.stringify(campaign));
    }, [ROOTS[c.app], c.id, table] as [string, string, string]);
  }, {tier: 'api', description: 'writes into the campaign\'s campaign.json the card layout a tracker keeps (`tilesViewerFormSketch`): the V-iD, Stage and MW fields, over a table of the name given'});

export const openByLink = When('user opens the campaign {string} by its link', async (page: Page, name: string) => {
  const c = campaignOf(page, name);
  await page.evaluate((path) => {
    grok.shell.route(path);
  }, `/apps/HitTriage/${APP_PATHS[c.app]}?campaignId=${c.id}`);
  const wanted = `TableView: ${c.name}`;
  let current = '';
  await expect.poll(async () => (current = await page.evaluate(() => `${grok.shell.v?.type}: ${grok.shell.v?.name}`)) === wanted,
    {message: `the campaign view of "${name}"`, timeout: pollMs(60000)}).toBe(true).catch(() => {
    throw new Error(`the campaign view of "${name}" (${wanted}) is not current; current: ${current}`);
  });
}, {tier: 'api', description: 'the app\'s address with ?campaignId=<id> (a deep link); done when the campaign\'s table view, named after the campaign, is current'});

/* The shell renders the crumbs of the current view ("Hit Design › <campaign>") from its parent view: the parent's
   crumb is there with the view, so the crumb of this campaign is the one shown while its view is current. */
export const backByCrumb = When('user goes back by the {string} crumb of the campaign {string}', async (page: Page, crumb: string, name: string) => {
  await expect.poll(() => page.evaluate(() => grok.shell.v?.name), {message: `the campaign view of "${name}" current`})
    .toBe(name);
  const link = page.locator('.d4-view-ribbon > .d4-ribbon-name:not(.grok-ns-entity-name)', {hasText: crumb})
    .filter({visible: true}).first();
  await expect(link, `the "${crumb}" crumb of "${name}"`).toBeVisible();
  await link.click();
}, {tier: 'ui', description: 'clicks the parent crumb the shell shows before the campaign\'s name in its view\'s ribbon, once that view is current'});

/* ---------------------------------------------------------------- what a save put on the server */

interface Saved {rows: number; version: number; status: string; name: string; edit: string[];
  values: {[col: string]: string[]}}

function savedCampaign(page: Page, c: Campaign): Promise<Saved | string> {
  return page.evaluate(async ([root, id]) => {
    const folder = `${root}/campaigns/${id}`;
    if (!await grok.dapi.files.exists(`${folder}/campaign.json`))
      return `no ${folder}/campaign.json`;
    const text = await grok.dapi.files.readAsText(`${folder}/campaign.json`);
    let json: any;
    try {
      json = JSON.parse(text);
    }
    catch {
      // read while a save rewrites it: a throw would end the poll
      return `${folder}/campaign.json is being written (${text.length} characters read)`;
    }
    const df = DG.DataFrame.fromCsv(await grok.dapi.files.readAsText(`${folder}/enriched_table.csv`));
    const values: {[col: string]: string[]} = {};
    for (const col of df.columns.toList())
      values[col.name] = Array.from({length: df.rowCount}, (_, i) => col.isNone(i) ? '' : String(col.get(i)));
    return {rows: df.rowCount, version: json.version ?? 0, status: json.status ?? '', name: json.friendlyName ?? '',
      edit: json.permissions?.edit ?? [], values};
  }, [ROOTS[c.app], c.id] as [string, string]);
}

async function pollSaved(page: Page, name: string, what: string, test: (s: Saved) => string): Promise<void> {
  const c = campaignOf(page, name);
  let seen = '';
  await expect.poll(async () => {
    const s = await savedCampaign(page, c);
    seen = typeof s === 'string' ? s : test(s);
    return seen;
  }, {message: `the campaign "${name}" as saved on the server: ${what}`, timeout: pollMs(30000)}).toBe('').catch(() => {
    throw new Error(`the campaign "${name}" on the server is not ${what}: ${seen}`);
  });
}

const cellOf = (s: Saved, column: string, row: number): string | null =>
  s.values[column] === undefined ? null : s.values[column][row - 1] ?? null;

export const savedRows = Then('the campaign {string} should be saved with {int} rows', (page: Page, name: string, rows: number) =>
  pollSaved(page, name, `${rows} rows`, (s) => s.rows === rows ? '' : `${s.rows} rows`),
{tier: 'api', description: 'its table on the share, read back until the save has written it'});

export const savedCell = Then('the campaign {string} should be saved with {string} in row {int} of the {string} column',
  (page: Page, name: string, value: string, row: number, column: string) =>
    pollSaved(page, name, `"${value}" in row ${row} of ${column}`, (s) => {
      const v = cellOf(s, column, row);
      return v === value ? '' : v === null ? `no ${column} column or no row ${row} (${s.rows} rows; ${Object.keys(s.values).join(', ')})` : `"${v}"`;
    }),
{tier: 'api', description: 'its table on the share, row counted from 1, the value as the CSV holds it; read back until the save has written it'});

export const savedCellMatching = Then('the campaign {string} should be saved with a value matching {string} in row {int} of the {string} column',
  (page: Page, name: string, pattern: string, row: number, column: string) =>
    pollSaved(page, name, `/${pattern}/ in row ${row} of ${column}`, (s) => {
      const v = cellOf(s, column, row);
      return v !== null && new RegExp(pattern).test(v) ? '' : v === null ? `no ${column} column or no row ${row} (${s.rows} rows)` : `"${v}"`;
    }),
{tier: 'api', description: 'as the step above, the value matched by the regular expression'});

export const savedStatus = Then('the campaign {string} should be saved with the status {string}', (page: Page, name: string, status: string) =>
  pollSaved(page, name, `status "${status}"`, (s) => s.status === status ? '' : `status "${s.status}"`),
{tier: 'api', description: 'campaign.json on the share, read back until the save has written it'});

export const savedName = Then('the campaign {string} should be saved with the name {string}', (page: Page, name: string, shown: string) =>
  pollSaved(page, name, `named "${shown}"`, (s) => s.name === shown ? '' : `named "${s.name}"`),
{tier: 'api', description: 'the friendly name campaign.json carries, read back until the save has written it'});

export const savedVersion = Then('the campaign {string} should be at version {int} on the server', (page: Page, name: string, version: number) =>
  pollSaved(page, name, `version ${version}`, (s) => s.version === version ? '' : `version ${s.version}`),
{tier: 'api', description: 'the version campaign.json carries (a save adds one)'});

export const sharingUserMayEdit = Then('the campaign {string} should let the sharing user edit it on the server', async (page: Page, name: string) => {
  const group = await page.evaluate(async (login) => (await grok.dapi.users.filter(`login = ${JSON.stringify(login)}`).first())?.group?.id ?? '',
    sharingLogin());
  if (!group)
    throw new Error(`no personal group of "${sharingLogin()}" on the server`);
  await pollSaved(page, name, `editable by ${sharingLogin()}`, (s) => s.edit.includes(group) ? '' : `edit groups ${s.edit.join(', ')}`);
}, {tier: 'api', description: 'the edit permissions of campaign.json hold the personal group of DATAGROK_SHARING_LOGIN (bddsecond)'});

/* ---------------------------------------------------------------- the update log (what the refresh poll reads) */

const updateLogs = new WeakMap<Page, Map<string, string>>();

function updateLogOf(page: Page, c: Campaign): Promise<string> {
  return page.evaluate(async ([a, id]) => {
    const conn = (await grok.dapi.connections.filter('name = "hitdesign"').list())
      .find((x: any) => x.nqName === 'HitTriage:hitdesign');
    const df = await conn.query('bdd-campaign-log', `select updated_at from hitdesign.update_logs where app_name = '${a}' and campaign_id = '${id}'`).apply();
    return df.rowCount > 0 ? df.col('updated_at').getString(0) : '';
  }, [c.app, c.id] as [string, string]);
}

export const rememberUpdateLog = When('user remembers the update log of the campaign {string}', async (page: Page, name: string) => {
  if (!updateLogs.has(page))
    updateLogs.set(page, new Map());
  updateLogs.get(page)!.set(name, await updateLogOf(page, campaignOf(page, name)));
}, {tier: 'api', description: 'the last update the hitdesign database logged for it (none before its first save), which the campaign view polls every 3 s'});

/* A save takes the campaign's lock before anything else and lets it go last (the release is what logs an update), so
   a lock gone after the save's balloon is the end of that save: the update log is read after it, never before. */
export const noLockHeld = Then('the campaign {string} should hold no lock on the server', async (page: Page, name: string) => {
  const c = campaignOf(page, name);
  await expect.poll(() => page.evaluate(async ([a, id]) => {
    const conn = (await grok.dapi.connections.filter('name = "hitdesign"').list())
      .find((x: any) => x.nqName === 'HitTriage:hitdesign');
    return (await conn.query('bdd-campaign-lock', `select campaign_id from hitdesign.campaign_locks where app_name = '${a}' and campaign_id = '${id}'`)
      .apply()).rowCount;
  }, [c.app, c.id] as [string, string]), {message: `the save lock of "${name}" let go`, timeout: pollMs(35000)}).toBe(0);
}, {tier: 'api', description: 'no row of the campaign in hitdesign.campaign_locks, polled: the save that took the lock has ended'});

export const updateLogAsRemembered = Then('the update log of the campaign {string} should be as remembered', async (page: Page, name: string) => {
  const before = updateLogs.get(page)?.get(name);
  if (before === undefined)
    throw new Error(`the update log of "${name}" was not remembered`);
  expect(await updateLogOf(page, campaignOf(page, name)), `the update log of "${name}" (remembered: "${before}")`).toBe(before);
}, {tier: 'api', description: 'read once, after an end signal the step before gave (the refused save\'s balloon)'});

/* ---------------------------------------------------------------- a save from another session */

export const otherSessionSaves = When('another session saves the campaign {string} with a new row of {string}', async (page: Page, name: string, smiles: string) => {
  const c = campaignOf(page, name);
  const error = await page.evaluate(async ([root, a, id, mol]) => {
    const lock = await grok.data.query('HitTriage:AcquireCampaignLock', {appName: a, campaignId: id, lockedBy: 'bdd second session'});
    if (!lock || lock.rowCount === 0)
      return 'the campaign lock was not acquired';
    const folder = `${root}/campaigns/${id}`;
    const json = JSON.parse(await grok.dapi.files.readAsText(`${folder}/campaign.json`));
    const df = DG.DataFrame.fromCsv(await grok.dapi.files.readAsText(`${folder}/enriched_table.csv`));
    df.rows.addNew(null, false);
    df.set('Molecule', df.rowCount - 1, mol);
    if (df.col('Stage') && json.template?.stages?.length)
      df.set('Stage', df.rowCount - 1, json.template.stages[0]);
    await grok.dapi.files.writeAsText(`${folder}/enriched_table.csv`, df.toCsv());
    json.version = (json.version ?? 0) + 1;
    json.rowCount = df.rowCount;
    json.filteredRowCount = df.rowCount;
    await grok.dapi.files.writeAsText(`${folder}/campaign.json`, JSON.stringify(json));
    const released = await grok.data.query('HitTriage:ReleaseCampaignLock', {appName: a, campaignId: id});
    return released?.rowCount > 0 ? '' : 'the campaign lock was not released';
  }, [ROOTS[c.app], c.app, c.id, smiles] as [string, App, string, string]);
  if (error)
    throw new Error(`another session's save of "${name}" failed: ${error}`);
}, {tier: 'api', description: 'what a save in another browser does, through the package\'s own lock queries and the files API: the lock taken, the row appended (the first stage) to the table, the version raised, the lock released, which logs the update the campaign view polls for'});

/* ---------------------------------------------------------------- files in the user's home */

const homePath = (page: Page, file: string): Promise<string> =>
  page.evaluate((f) => `${grok.shell.user.project.name}:Home/${f}`, file);

export const userFileHolds = Given('the user\'s file {string} holds the table:', async (page: Page, file: string, table: string[][]) => {
  const path = await homePath(page, file);
  const remove = async () => {
    await page.evaluate(async (p) => {
      if (await grok.dapi.files.exists(p))
        await grok.dapi.files.delete(p);
    }, path);
    await expect.poll(() => page.evaluate((p) => grok.dapi.files.exists(p), path),
      {message: `${path} deleted`, timeout: pollMs(30000)}).toBe(false);
  };
  await remove();
  atFeatureEnd(page, remove);
  await page.evaluate(([p, rows]) => grok.dapi.files.writeAsText(p, rows.map((r) => r.join(',')).join('\n') + '\n'),
    [path, table] as [string, string[][]]);
}, {tier: 'api', description: 'a CSV of the rows given (the first one the header) in the account\'s home folder, written through the files API; deleted at feature end, read back'});

export const holdsUserFilePath = Then('{element} should hold the path of the user\'s file {string}',
  async (page: Page, target: ElementRef, file: string) => {
    const path = await homePath(page, file);
    await expect((await locate(page, target)).locator('input[type="text"]').first(), `the share path ${path}`).toHaveValue(path);
  }, {tier: 'ui', description: 'the text of the file input\'s box: the file\'s share path (<home project>:Home/<file>)'});

export const chooseProjectFile = When('user clicks on {element} and chooses the {string} file of the project',
  async (page: Page, target: ElementRef, file: string) => {
    const chooser = page.waitForEvent('filechooser', {timeout: pollMs(15000)});
    await gestures.click(page, target);
    await (await chooser).setFiles(resolve(process.env.BDD_ROOT ?? process.cwd(), file));
  }, {tier: 'ui', description: 'the click opens the browser\'s file chooser, which is answered with a file of the bdd project'});

type MergeKept = {filePath?: string; autoMergeOnOpen?: boolean} | null;

const mergeConfigOf = (page: Page, c: Campaign): Promise<MergeKept> =>
  page.evaluate(async ([root, id]) => JSON.parse(await grok.dapi.files.readAsText(`${root}/campaigns/${id}/campaign.json`))
    .mergeConfig ?? null, [ROOTS[c.app], c.id] as [string, string]);

async function pollMergeConfig(page: Page, name: string, what: string, test: (cfg: MergeKept) => boolean): Promise<void> {
  const c = campaignOf(page, name);
  let seen: MergeKept = null;
  await expect.poll(async () => {
    seen = await mergeConfigOf(page, c);
    return test(seen);
  }, {message: `the merge kept with the campaign "${name}": ${what}`, timeout: pollMs(30000)}).toBe(true).catch(() => {
    throw new Error(`the campaign "${name}" on the server does not ${what}: its mergeConfig is ${JSON.stringify(seen)}`);
  });
}

export const savedToMergeUserFile = Then('the campaign {string} should be saved to merge the user\'s file {string} on every open',
  async (page: Page, name: string, file: string) => {
    const path = await homePath(page, file);
    await pollMergeConfig(page, name, `merge ${path} on every open`, (cfg) => cfg?.filePath === path && cfg.autoMergeOnOpen === true);
  }, {tier: 'api', description: 'its campaign.json\'s mergeConfig names the file\'s share path with Auto-merge on campaign open, read back until the save wrote it'});

export const savedToMergeNoFile = Then('the campaign {string} should be saved to merge no file on open', (page: Page, name: string) =>
  pollMergeConfig(page, name, 'keep a merge with no file to read on open',
    (cfg) => cfg !== null && cfg.filePath === undefined && cfg.autoMergeOnOpen === false),
{tier: 'api', description: 'its campaign.json keeps the merge\'s choices but no file path, and Auto-merge on campaign open off'});

export const enterUserFilePath = When('user enters the path of the user\'s file {string} into {element}',
  async (page: Page, file: string, target: ElementRef) => gestures.typeInto(page, target, await homePath(page, file), true),
  {tier: 'ui', description: 'types the file\'s share path (<home project>:Home/<file>) and presses Enter, as a person pastes it'});

/* ---------------------------------------------------------------- reads of the page */

export const valueMatching = Then('the value of {string} column in row {int} should match {string}',
  async (page: Page, column: string, row: number, pattern: string) => {
    await expect.poll(() => page.evaluate(([c, r]) => {
      const col = grok.shell.t?.col(c);
      return col ? r <= col.length ? col.isNone(r - 1) ? '' : String(col.get(r - 1)) : `(the table has ${col.length} rows)` : `(no ${c} column)`;
    }, [column, row] as [string, number]), {message: `"${column}" in row ${row} of the current table`}).toMatch(new RegExp(pattern));
  }, {description: 'the current table\'s cell, rows from 1, matched by the regular expression; polled'});

/* The sketcher applies a typed molecule on its own time (OpenChemLib loads it into its editor after the Enter), and
   OK writes what it holds then: an OK right after the Enter wrote the old molecule. The step ends on the sketcher's
   own state. */
export const enterMoleculeInSketcher = When('user enters the molecule {string} into sketcher dialog', async (page: Page, smiles: string) => {
  const dialog = page.locator('.d4-dialog:has(input[placeholder^="SMILES"])').filter({visible: true}).last();
  const input = dialog.locator('input[placeholder^="SMILES"]');
  await input.click();
  await input.fill(smiles);
  await input.press('Enter');
  await expect.poll(() => dialog.evaluate((d, s) => {
    for (const e of Array.from(d.querySelectorAll('[data-widget]'))) {
      const w = DG.Widget.find(e);
      if (w && typeof w.getSmiles === 'function') {
        const now = String(w.getSmiles() ?? '');
        return now === '' ? '(empty)' : grok.chem.convert(now, DG.chem.Notation.Unknown, DG.chem.Notation.Smiles) ===
          grok.chem.convert(s, DG.chem.Notation.Unknown, DG.chem.Notation.Smiles) ? '' : now;
      }
    }
    return '(no sketcher in the dialog)';
  }, smiles), {message: `the sketcher holding ${smiles}`}).toBe('');
}, {tier: 'ui', description: 'types the SMILES into the molecule field of the open sketcher dialog and presses Enter; done when the sketcher holds that molecule (DG.chem.Sketcher.getSmiles, compared canonically), so OK writes it'});

/* The summary, read-only form rows or a key-value table (ui.tableFromMap), is read as label → shown text, the editable
   inputs below it left out. */
export const summaryIs = Then('the summary in {element} should be:', async (page: Page, target: ElementRef, rows: string[][]) => {
  const root = (await locate(page, target)).filter({visible: true}).first();
  const read = (): Promise<[string, string][]> => root.evaluate((r) => {
    const table = r.querySelector('table.d4-item-table');
    if (table) {
      return Array.from(table.querySelectorAll(':scope > tbody > tr'))
        .map((tr) => Array.from(tr.children).map((td) => (td.textContent ?? '').trim()))
        .map((cells) => [cells[0] ?? '', cells[1] ?? ''] as [string, string]);
    }
    return Array.from(r.querySelectorAll('label')).filter((l) => (l as HTMLElement).offsetWidth > 0)
      .map((l) => {
        const field = (l.htmlFor ? document.getElementById(l.htmlFor) : null) as HTMLInputElement | null;
        const text = field && !field.readOnly && 'value' in field ? null :
          field ? String(field.value ?? field.textContent ?? '') : (l.parentElement?.textContent ?? '').replace(l.textContent ?? '', '');
        return text === null ? null : [(l.textContent ?? '').trim(), text.trim()] as [string, string];
      }).filter((x): x is [string, string] => x !== null);
  });
  // a value "*" is not compared (it changes by design), one ending in "…" is a beginning
  const holds = (got: [string, string][]): boolean => got.length === rows.length && rows.every(([k, v], i) =>
    got[i][0] === k && (v === '*' || (v.endsWith('…') ? got[i][1].startsWith(v.slice(0, -1)) : got[i][1] === v)));
  let got: [string, string][] = [];
  await expect.poll(async () => holds(got = await read()), {message: `the summary of ${target.phrase}`}).toBe(true).catch(() => {
    throw new Error(`the summary of ${target.phrase} is not ${rows.map(([k, v]) => `${k} = ${v}`).join('; ')}: ` +
      `it shows ${got.map(([k, v]) => `${k} = ${v}`).join('; ')}`);
  });
}, {tier: 'ui', description: 'the read-only summary in the element (a key-value table, or read-only form rows): exactly these keys, in this order, each showing that text — "*" any text, a value ending in "…" the text\'s beginning'});

/* A Molecule field's structure is the row's value: inside its row, which is as tall as the picture, starting where the
   text values start, clear of the next row (it hung over the next row from a 28 px row before). */
export const structureInRow = Then('the structure of the {string} row in {element} should sit in the value column of its row',
  async (page: Page, caption: string, target: ElementRef) => {
    const root = (await locate(page, target)).filter({visible: true}).first();
    let seen = '';
    await expect.poll(async () => seen = await root.evaluate((r, c) => {
      const rows = Array.from(r.querySelectorAll<HTMLElement>('.u2-input-root')).filter((x) => x.getClientRects().length > 0);
      const at = rows.findIndex((x) => x.querySelector(':scope > .u2-input-label')?.textContent?.trim() === c);
      if (at < 0)
        return `no "${c}" row`;
      const box = (el: Element | null | undefined) => el?.getBoundingClientRect() ?? null;
      const row = box(rows[at])!;
      const picture = box(rows[at].querySelector('canvas, svg, img'));
      const next = box(rows[at + 1]);
      const value = box(rows.find((x) => x.dataset.u2 === 'text-input')?.querySelector('[data-u2-part="editor"]'));
      if (picture === null || picture.height < 2)
        return 'no structure drawn';
      const problems: string[] = [];
      if (picture.top < row.top - 1 || picture.bottom > row.bottom + 1 || picture.right > row.right + 1)
        problems.push(`the structure ${Math.round(picture.top)}–${Math.round(picture.bottom)} is outside its row ${Math.round(row.top)}–${Math.round(row.bottom)}`);
      if (next !== null && picture.bottom > next.top + 1)
        problems.push('it covers the next row');
      if (value !== null && Math.abs(picture.left - value.left) > 8)
        problems.push(`it starts at x ${Math.round(picture.left)}, the values at ${Math.round(value.left)}`);
      return problems.join('; ');
    }, caption), {message: `the structure of the "${caption}" row in ${target.phrase}`, timeout: pollMs(15000)}).toBe('')
      .catch(() => {
        throw new Error(`the structure of the "${caption}" row in ${target.phrase} is misplaced: ${seen}`);
      });
  }, {tier: 'ui', description: 'the drawn structure lies inside its form row (the row grew to it), clear of the next row, and starts within 8 px of where the text values start'});

/* Status suggestions, a u2 suggest input's list or a platform popup menu, are read and picked by the text of their
   items. */
const SUGGESTIONS = '.d4-menu-popup .d4-menu-item-label, .u2-overlay [data-u2-part="popup"] [role="option"]';

export const suggests = Then('{element} should suggest {string}', async (page: Page, target: ElementRef, text: string) => {
  await locate(page, target);
  await expect.poll(() => page.evaluate((sel) => Array.from(document.querySelectorAll(sel))
    .filter((e) => (e as HTMLElement).offsetWidth > 0).map((e) => (e.textContent ?? '').trim()), SUGGESTIONS),
  {message: `the suggestions shown for ${target.phrase}`}).toContain(text);
}, {tier: 'ui', description: 'the suggestions open under the input (a popup menu, or a suggest input\'s list) hold an item of exactly that text'});

export const takeSuggestion = When('user takes the suggestion {string} of {element}', async (page: Page, text: string, target: ElementRef) => {
  await locate(page, target);
  await page.locator(SUGGESTIONS).filter({visible: true}).filter({hasText: new RegExp(`^\\s*${text.replace(/[.*+?^${}()|[\]\\]/g, '\\$&')}\\s*$`)})
    .first().click();
}, {tier: 'ui', description: 'clicks the suggestion of exactly that text shown for the input'});

/* Edit Permissions is an access grid: a user or group is added from the picker under its rows, then
   given its capabilities in its row. The feature names the outcome; this step is the one place that knows the form. */
export const grantSharingUserEdit = When('user lets the sharing user edit in {element}', async (page: Page, target: ElementRef) => {
  const dialog = (await locate(page, target)).filter({visible: true}).first();
  const grid = dialog.locator('[data-u2="access-grid"]');
  await expect(grid, `the access grid of ${target.phrase}`).toBeVisible({timeout: pollMs(15000)});
  const shown = new RegExp(sharingShownName(), 'i');
  await gestures.typeVerified(grid.locator('input.u2-typeahead-input'), sharingLogin().split('@')[0],
    `the add picker of ${target.phrase}`);
  const option = page.locator('[role="option"]').filter({has: page.locator('.u2-typeahead-user-name', {hasText: shown})}).first();
  await expect(option, `the "${sharingLogin()}" suggestion of the add picker`).toBeVisible({timeout: pollMs(15000)});
  await option.click();
  const canEdit = grid.locator('tr.u2-access-grid-row').filter({hasText: shown}).locator('input[data-u2-key="edit"]');
  await expect(canEdit, 'the Can Edit box of the second account\'s row').toBeVisible({timeout: pollMs(15000)});
  await canEdit.check();
  await expect(canEdit, 'the Can Edit box of the second account\'s row').toBeChecked();
}, {tier: 'ui', description: 'gives the second account (DATAGROK_SHARING_LOGIN) edit access in the Edit Permissions grid: adds it from the picker\'s suggestions, then checks Can Edit in its row'});

/* Modify Stages is a records input: a row per stage, "Add stage" under them. */
export const addStage = When('user adds the stage {string} in {element}', async (page: Page, stage: string, target: ElementRef) => {
  const dialog = (await locate(page, target)).filter({visible: true}).first();
  const records = dialog.locator('[data-u2="records-input"]');
  const cells = records.locator('[role="row"] input');
  const before = await cells.count();
  await records.getByText('Add stage', {exact: true}).click();
  await expect(cells, `a new row in ${target.phrase}`).toHaveCount(before + 1);
  await gestures.typeVerified(cells.last(), stage, `the new row of ${target.phrase}`);
}, {tier: 'ui', description: 'adds a row with "Add stage" and types the stage into it'});

export const namesSharingUser = Then('{element} should name the sharing user', async (page: Page, target: ElementRef) => {
  await expect((await locate(page, target)).filter({visible: true}).first(), `${target.phrase}, naming the second account`)
    .toContainText(new RegExp(sharingShownName(), 'i'));
}, {tier: 'ui', description: 'the element shows the second account (DATAGROK_SHARING_LOGIN) as the platform names it'});

/* A row of the Edit Permissions grid starts with who it is: a group's two-person icon, a user's picture or initials. */
async function principalPicture(page: Page, target: ElementRef, who: RegExp, kind: 'group' | 'user'): Promise<void> {
  const row = (await locate(page, target)).filter({visible: true}).first()
    .locator('[data-u2="access-grid"] tr.u2-access-grid-row')
    .filter({has: page.locator('.u2-access-grid-principal', {hasText: who})}).first();
  await expect(row.locator(`.u2-access-grid-icon > .u2-avatar${kind === 'group' ? '.u2-avatar-group' : ':not(.u2-avatar-group)'}`),
    `the ${kind} picture in the ${who} row of ${target.phrase}`).toBeVisible({timeout: pollMs(15000)});
}

export const groupPicture = Then('the {string} row of {element} should show a group icon',
  (page: Page, name: string, target: ElementRef) =>
    principalPicture(page, target, new RegExp(`^${name.replace(/[.*+?^${}()|[\]\\]/g, '\\$&')}$`), 'group'),
  {tier: 'ui', description: 'the access grid row of that group starts with the group icon'});

export const sharingUserPicture = Then('the sharing user\'s row of {element} should show the user\'s picture',
  (page: Page, target: ElementRef) => principalPicture(page, target, new RegExp(sharingShownName(), 'i'), 'user'),
  {tier: 'ui', description: 'the access grid row of the second account (DATAGROK_SHARING_LOGIN) starts with its picture or initials'});

export const ribbonLastGroup = Then('the last group of the ribbon should list {string}', async (page: Page, list: string) => {
  const want = list.split(/\s*,\s*/).filter((s) => s.length > 0);
  await expect.poll(() => page.evaluate(() => {
    const panels = Array.from(document.querySelectorAll('.d4-ribbon-panel')).filter((p) => (p as HTMLElement).offsetWidth > 0);
    const last = panels[panels.length - 1];
    const group = last?.querySelector('[data-u2="ribbon-group"]');
    if (group)
      return (DG.U2.Control.forElement(group).getItems().list as {key: string}[]).map((i) => i.key);
    return last ? Array.from(last.querySelectorAll(':scope > .d4-ribbon-item > *'))
      .map((e) => (e.getAttribute('aria-label') ?? e.textContent ?? '').trim()).filter((s) => s.length > 0) : [];
  }), {message: 'the items of the ribbon\'s last group, by their names'}).toEqual(want);
}, {tier: 'ui', description: 'the shown ribbon\'s last panel, item by item in order: a u2 ribbon group\'s items by their keys (an icon\'s name or tooltip, a button\'s text, folded ones included), else each element by its aria-label or text; an item with no name (a separator) is left out'});

export const closeViewByName = When('user closes the {string} view', async (page: Page, name: string) => {
  await page.evaluate((n) => {
    const view = Array.from(grok.shell.views as Iterable<any>).find((v) => v.name === n);
    if (!view)
      throw new Error(`no "${n}" view; open: ${Array.from(grok.shell.views as Iterable<any>).map((v) => v.name).join(', ')}`);
    view.close();
  }, name);
  await expect.poll(() => page.evaluate((n) => Array.from(grok.shell.views as Iterable<any>).some((v) => v.name === n), name),
    {message: `the "${name}" view closed`}).toBe(false);
}, {tier: 'api', description: 'an open view by its name, current or not — what the ✕ of its tab does where the tabs are shown'});
