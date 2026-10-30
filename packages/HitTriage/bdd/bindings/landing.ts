/* The landings of the four HitTriage apps (Hit Design, PeptiHit, Hit Triage, PepTriage) as the landing and template
   features use them: the campaigns table (its rows, the icons above it, its groups, the row actions), the template
   editor's own inputs, and the fixtures these features stand on.

   A campaign is two files on the AppData share, campaign.json and enriched_table.csv, written here in the shape the
   app writes them; Hit Design and PeptiHit also keep the time of the last save in the hitdesign database
   (update_logs), which a fixture writes when its row says how long ago it was modified. A template is one file,
   <app>/templates/<name>.json.

   Every fixture carries the feature's {time}: a campaign code is a family in capitals, the time and a number
   ("BDDLND1761512345678-1"), a template name ends in the time, a template key is the family and the time. It goes
   at feature end with what the feature made from it (the campaigns a template's key started through the UI, their
   rows in the hitdesign database), read back until gone; the same family left by a run older than an hour is swept
   when the next one starts, so a second run of the feature at the same time keeps its own.

   The app keeps the campaign list it read for the life of the page (HTPackage.campaignsCache): campaigns
   put on the server before the app opens drop that copy, as a fresh page would not have it; campaigns saved "from
   elsewhere" leave it alone, which is what the refresh claims are about. */
import {Page} from '@playwright/test';
import {Given, Then, When, element, kind} from '@datagrok-libraries/bdd';
import {type ElementRef, atFeatureEnd, exactText, expect, gestures, locate, pollMs, STALE_AFTER_MS}
  from '@datagrok-libraries/bdd/runtime';

declare const grok: any;
declare const DG: any;

type App = 'Hit Design' | 'PeptiHit' | 'Hit Triage' | 'PepTriage';
const APPS: App[] = ['Hit Design', 'PeptiHit', 'Hit Triage', 'PepTriage'];
const isDesign = (app: App): boolean => app === 'Hit Design' || app === 'PeptiHit';
const share = (app: App): string => `System:AppData/HitTriage/${app}`;

function appOf(name: string): App {
  if (!APPS.includes(name as App))
    throw new Error(`"${name}" is not a HitTriage app; the apps are ${APPS.join(', ')}`);
  return name as App;
}

// ---- the nouns -------------------------------------------------------------------------------------------------

const TABLE = '[data-u2="catalog-table"][data-u2-name="campaigns"]';

element('campaigns table', {selector: TABLE,
  description: 'the table under "Continue Campaign" on an app\'s landing'});
element('campaign columns icon', {selector: '[data-u2-name="campaigns-columns"]',
  description: 'the eye above the campaigns table ("Edit visible columns in campaigns table"): a menu of the columns to show'});
element('campaign grouping icon', {selector: '[data-u2-name="campaigns-group"]',
  description: 'the layer-group icon above the campaigns table ("Group Campaigns. Current: …"): a menu of the groupings'});
element('campaigns refresh icon', {selector: '[data-u2-name="campaigns-refresh"]',
  description: 'the sync icon above the campaigns table ("Refresh campaigns table")'});
kind('campaign row', {
  selector: `${TABLE} tr.u2-table-row`,
  match: ['label'],
  labelSelector: ':scope > td:first-child',
  parts: {'permissions icon': '[aria-label="Manage campaign permissions"]', 'delete icon': '[aria-label="Delete campaign"]',
    'name link': 'a'},
  description: 'a row of the campaigns table by the campaign\'s name ("BDD Alpha" campaign row); "permissions icon of …" and "delete icon of …" are its row actions',
});
kind('campaign group', {
  selector: `${TABLE} tr.u2-table-section`,
  match: ['label'],
  labelSelector: '[data-u2-part="item-label"]',
  description: 'the header of a group of the campaigns table by its label ("In Progress" campaign group); expanded or collapsed',
});
element('HELM editor dialog', {selector: '.d4-dialog-full-screen:has([data-testid="app-root"])',
  description: 'the full-screen HELM editor (Helm:editMoleculeCell) a new PeptiHit campaign opens on its first row; no title, CANCEL and OK'});
const TEMPLATE_EDITOR = '[data-u2-name="template-editor"]';
element('template name input', {selector: `${TEMPLATE_EDITOR} [data-u2-name="name"], .ui-form [name="input-host-Name"]:not(.ui-items-grid *)`,
  description: 'the Name input of the template editor\'s Details (Hit Design and PeptiHit: u2; Hit Triage and PepTriage: the platform form, where a field row has a Name input too)'});
element('template key input', {selector: `${TEMPLATE_EDITOR} [data-u2-name="key"], .ui-form [name="input-host-Key"]`,
  description: 'the Key input of the template editor\'s Details: the prefix of the campaign codes'});

// ---- the accounts ----------------------------------------------------------------------------------------------

interface Account {id: string; name: string; group: string}

function accounts(page: Page): Promise<{me: Account, second: Account | null}> {
  return page.evaluate(async (login) => {
    const me = DG.User.current();
    const second = login ? await grok.dapi.users.filter(`login = "${login}"`).first() : null;
    return {me: {id: me.id, name: me.friendlyName, group: me.group.id},
      second: second ? {id: second.id, name: second.friendlyName, group: second.group.id} : null};
  }, process.env.DATAGROK_SHARING_LOGIN ?? null);
}

// ---- the fixtures: names, families, sweeps ---------------------------------------------------------------------

const CODE = /^([A-Z]+)(\d{13})-(\d+)$/;
const KEY = /^([A-Z]+)(\d{13})$/;
const TEMPLATE_NAME = /^(.*\S)\s+(\d{13})$/;
const stale = (time: string | number): boolean => Date.now() - Number(time) > STALE_AFTER_MS;

function codeParts(code: string): {family: string, prefix: string} {
  const m = CODE.exec(code);
  if (!m)
    throw new Error(`a fixture campaign's code is a family in capitals, the feature's {time} and a number ("BDDLND{time}-1"), not "${code}"`);
  return {family: m[1], prefix: m[1] + m[2]};
}

function keyParts(key: string): {family: string, prefix: string} {
  const m = KEY.exec(key);
  if (!m)
    throw new Error(`a fixture template's key is a family in capitals and the feature's {time} ("BDDNEW{time}"), not "${key}"`);
  return {family: m[1], prefix: key};
}

function templateFamily(name: string): string {
  const m = TEMPLATE_NAME.exec(name);
  if (!m)
    throw new Error(`a fixture template's name ends in the feature's {time} ("BDD Picker {time}"), not "${name}"`);
  return m[1];
}

/** What a sweep removes: the campaign folders of this run's prefixes, the stale ones of the families, and their rows in
 * the hitdesign database. Returns what is left ('' when all of it is gone). */
async function sweepCampaigns(page: Page, app: App, families: string[], prefixes: string[]): Promise<string> {
  const root = `${share(app)}/campaigns`;
  const names: string[] = await page.evaluate(async (r) => (await grok.dapi.files.list(r))
    .filter((f: any) => f.isDirectory).map((f: any) => f.name), root);
  const doomed = names.filter((n) => prefixes.some((p) => n.startsWith(`${p}-`)) ||
    (CODE.test(n) && families.includes(CODE.exec(n)![1]) && stale(CODE.exec(n)![2])));
  const cutoff = Date.now() - STALE_AFTER_MS;
  const where = [...prefixes.map((p) => `campaign_id like '${p}-%'`),
    ...families.map((f) => `(campaign_id ~ '^${f}[0-9]{13}-' and substring(campaign_id from '^${f}([0-9]{13})')::bigint < ${cutoff})`)];
  return page.evaluate(async ([r, ids, a, w, design]) => {
    for (const id of ids)
      await grok.dapi.files.delete(`${r}/${id}`);
    const left: string[] = [];
    for (const id of ids) {
      if (await grok.dapi.files.exists(`${r}/${id}`))
        left.push(`folder ${id}`);
    }
    if (design && w.length > 0) {
      const conn = (await grok.dapi.connections.filter('name = "hitdesign"').list())
        .find((c: any) => c.nqName === 'HitTriage:hitdesign');
      const cond = `app_name = '${a}' and (${w.join(' or ')})`;
      for (const table of ['update_logs', 'campaign_locks', 'campaign_vids'])
        await conn.query('bdd-landing-cleanup', `delete from hitdesign.${table} where ${cond} returning campaign_id`).apply();
      for (const table of ['update_logs', 'campaign_locks', 'campaign_vids']) {
        const rows = (await conn.query('bdd-landing-cleanup', `select campaign_id from hitdesign.${table} where ${cond}`).apply()).rowCount;
        if (rows > 0)
          left.push(`${rows} row(s) in hitdesign.${table}`);
      }
    }
    return left.join(', ');
  }, [root, doomed, app, where, isDesign(app)] as [string, string[], string, string[], boolean]);
}

async function sweepTemplates(page: Page, app: App, families: string[], names: string[]): Promise<{left: string, keys: string[]}> {
  const root = `${share(app)}/templates`;
  const files: string[] = await page.evaluate(async (r) => (await grok.dapi.files.list(r))
    .filter((f: any) => !f.isDirectory && f.name.endsWith('.json')).map((f: any) => f.name.slice(0, -5)), root);
  const doomed = files.filter((n) => names.includes(n) ||
    (TEMPLATE_NAME.test(n) && families.includes(TEMPLATE_NAME.exec(n)![1]) && stale(TEMPLATE_NAME.exec(n)![2])));
  return page.evaluate(async ([r, ds]) => {
    const keys: string[] = [];
    for (const n of ds) {
      try {
        keys.push(String(JSON.parse(await grok.dapi.files.readAsText(`${r}/${n}.json`)).key ?? ''));
      }
      catch {
        // a file that no longer parses goes all the same
      }
      await grok.dapi.files.delete(`${r}/${n}.json`);
    }
    const left: string[] = [];
    for (const n of ds) {
      if (await grok.dapi.files.exists(`${r}/${n}.json`))
        left.push(`template ${n}`);
    }
    return {left: left.join(', '), keys: keys.filter((k) => k !== '')};
  }, [root, doomed] as [string, string[]]);
}

/** The list of campaigns the page read is the app's own cache for the life of the page: after a feature removes its
 * campaigns, a later feature on the same page would still see them listed. */
function dropPageCopy(page: Page, app: App): Promise<void> {
  return page.evaluate((a) => {
    const cache = (window as any).hittriage?._package?.campaignsCache;
    if (cache)
      delete cache[a];
  }, app);
}

async function readBack(what: string, sweep: () => Promise<string>): Promise<void> {
  let left = '';
  await expect.poll(async () => left = await sweep(), {message: `${what} removed from the server`, timeout: pollMs(30000)})
    .toBe('').catch((e) => {
      throw new Error(`${what} is still on the server: ${left || e?.message}`);
    });
}

/** The campaigns of a template key that a feature started through the UI go with the template. */
function campaignPrefixesOf(keys: string[]): {families: string[], prefixes: string[]} {
  const parts = keys.filter((k) => KEY.test(k)).map(keyParts);
  return {families: [...new Set(parts.map((p) => p.family))], prefixes: parts.map((p) => p.prefix)};
}

// ---- the fixtures: campaigns -----------------------------------------------------------------------------------

const SMILES = ['CCO', 'c1ccccc1', 'CC(=O)O', 'CCN', 'CCCC', 'c1ccncc1', 'CC(C)O', 'OCCO', 'CS(C)=O', 'C1CCCCC1'];
const HELMS = ['PEPTIDE1{A.G}$$$$V2.0', 'PEPTIDE1{G.A.K}$$$$V2.0', 'PEPTIDE1{L.A}$$$$V2.0', 'PEPTIDE1{K.K.G}$$$$V2.0'];
const NO_COMPUTE = {descriptors: {enabled: false, args: []}, functions: [], scripts: [], queries: []};
const COLUMNS = ['code', 'name', 'status', 'created', 'molecules', 'author', 'can view', 'can edit', 'template', 'fields', 'modified'];

interface CampaignFile {code: string, json: any, csv: string, modified: number | null}

/** "Owner = Ann; Due = date 2026-11-15": a Date field is stored the way the new-campaign form stores it. */
function parseFields(text: string): {values: {[k: string]: any}, defs: {name: string, type: string, required: boolean}[]} {
  const values: {[k: string]: any} = {};
  const defs: {name: string, type: string, required: boolean}[] = [];
  for (const part of text.split(';').map((p) => p.trim()).filter((p) => p !== '')) {
    const [name, ...rest] = part.split('=');
    const value = rest.join('=').trim();
    const date = /^date\s+(\d{4}-\d{2}-\d{2})$/.exec(value);
    if (date) {
      const d = new Date(`${date[1]}T12:00:00Z`);
      values[name.trim()] = {_value: d.getTime(), isUtc: false, date: d.toISOString()};
    }
    else
      values[name.trim()] = value;
    defs.push({name: name.trim(), type: date ? 'Date' : 'String', required: false});
  }
  return {values, defs};
}

function campaignFiles(app: App, table: string[][], who: {me: Account, second: Account | null}): CampaignFile[] {
  const [header, ...rows] = table;
  const unknown = header.filter((h) => !COLUMNS.includes(h));
  if (unknown.length > 0 || !header.includes('code'))
    throw new Error(`the campaigns table takes the columns ${COLUMNS.join(', ')} (code required); it has ${unknown.join(', ')}`);
  const today = new Date();
  return rows.map((cells) => {
    const r: {[k: string]: string} = {};
    header.forEach((h, i) => r[h] = cells[i]?.trim() ?? '');
    const code = r['code'];
    const {prefix} = codeParts(code);
    const author = (r['author'] || 'me') === 'me' ? who.me : who.second;
    if (!author)
      throw new Error(`"${code}" is authored by the sharing user, and the run has no second account (DATAGROK_SHARING_LOGIN)`);
    const groups = (v: string) => (v || 'everyone') === 'everyone' ? [DG_ALL_USERS] : [author.group];
    const n = Number(r['molecules'] || '1');
    const gone = /^gone:\s*/.test(r['template'] ?? '');
    const templateName = (r['template'] ?? '').replace(/^gone:\s*/, '') || prefix;
    const {values, defs} = parseFields(r['fields'] ?? '');
    const template: any = {name: templateName, key: prefix, campaignFields: defs, compute: NO_COMPUTE};
    const createDate = r['created'] || `${today.getFullYear()}/${today.getMonth() + 1}/${today.getDate()}`;
    const common: any = {
      name: code, ...(r['name'] ? {friendlyName: r['name']} : {}), templateName, status: r['status'] || 'In Progress',
      createDate, campaignFields: values, rowCount: n, filteredRowCount: n,
      authorUserId: author.id, authorUserFriendlyName: author.name, lastModifiedUserName: author.name,
      permissions: {view: groups(r['can view']), edit: groups(r['can edit'])},
    };
    const modified = r['modified'] ? Number(r['modified']) : null;
    if (isDesign(app)) {
      const helm = app === 'PeptiHit';
      template.stages = ['Design'];
      const cols = [...(helm ? ['Helm'] : []), 'Molecule', 'Stage', 'V-iD'];
      const csv = [cols.join(','), ...Array.from({length: n}, (_, i) =>
        [...(helm ? [HELMS[i % HELMS.length]] : []), SMILES[i % SMILES.length], 'Design', ''].join(','))].join('\n') + '\n';
      const json = {...common,
        columnSemTypes: {...(helm ? {Helm: 'Macromolecule'} : {}), Molecule: 'Molecule', Stage: null, 'V-iD': null},
        savePath: `System.AppData/HitTriage/${app}/campaigns/${code}/enriched_table.csv`,
        columnTypes: {...(helm ? {Helm: 'string'} : {}), Molecule: 'string', Stage: 'string', 'V-iD': 'string'},
        columnTags: helm ? {Helm: {'units': 'helm', 'cell.renderer': 'helm', '.alphabetIsMultichar': 'true'}} : {},
        version: 1, ...(gone ? {} : {template})};
      return {code, json, csv, modified};
    }
    const pep = app === 'PepTriage';
    const cols = [...(pep ? ['Helm'] : []), 'Molecule', 'Selected hits'];
    const csv = [cols.join(','), ...Array.from({length: n}, (_, i) =>
      [...(pep ? [HELMS[i % HELMS.length]] : []), SMILES[i % SMILES.length], 'false'].join(','))].join('\n') + '\n';
    const json = {...common, filters: [],
      ingest: {type: 'File', query: `${share(app)}/campaigns/${code}/enriched_table.csv`, molColName: 'Molecule'},
      columnSemTypes: {...(pep ? {Helm: 'Macromolecule'} : {}), Molecule: 'Molecule', 'Selected hits': null},
      columnTypes: {...(pep ? {Helm: 'string'} : {}), Molecule: 'string', 'Selected hits': 'bool'},
      ...(gone ? {} : {template: {...template, dataSourceType: 'File', ...(pep ? {sequenceColumnName: 'Helm'} : {})}})};
    return {code, json, csv, modified};
  });
}

const DG_ALL_USERS = 'a4b45840-9a50-11e6-9cc9-8546b8bf62e6';

async function writeCampaigns(page: Page, app: App, files: CampaignFile[], freshList: boolean): Promise<void> {
  // written oldest first, so the share's own times agree with the "modified" column
  const ordered = [...files].sort((a, b) => (b.modified ?? 0) - (a.modified ?? 0));
  await page.evaluate(async ([r, a, fs, design, drop]) => {
    let conn: any = null;
    for (const f of fs) {
      await grok.dapi.files.writeAsText(`${r}/${f.code}/enriched_table.csv`, f.csv);
      await grok.dapi.files.writeAsText(`${r}/${f.code}/campaign.json`, JSON.stringify(f.json));
      if (design && f.modified !== null) {
        conn ??= (await grok.dapi.connections.filter('name = "hitdesign"').list()).find((c: any) => c.nqName === 'HitTriage:hitdesign');
        await conn.query('bdd-landing-fixture', `insert into hitdesign.update_logs (app_name, campaign_id, updated_at) ` +
          `values ('${a}', '${f.code}', now() - interval '${Number(f.modified)} minutes') on conflict (app_name, campaign_id) ` +
          `do update set updated_at = excluded.updated_at returning campaign_id`).apply();
      }
    }
    const cache = (window as any).hittriage?._package?.campaignsCache;
    if (drop && cache)
      delete cache[a];
  }, [`${share(app)}/campaigns`, app, ordered, isDesign(app), freshList] as [string, string, CampaignFile[], boolean, boolean]);
  const missing: string[] = await page.evaluate(async ([r, codes]) => {
    const out: string[] = [];
    for (const c of codes) {
      if (!await grok.dapi.files.exists(`${r}/${c}/campaign.json`))
        out.push(c);
    }
    return out;
  }, [`${share(app)}/campaigns`, files.map((f) => f.code)] as [string, string[]]);
  if (missing.length > 0)
    throw new Error(`the ${app} campaigns ${missing.join(', ')} did not reach the server`);
}

async function putCampaigns(page: Page, appName: string, table: string[][], freshList: boolean): Promise<void> {
  const app = appOf(appName);
  const files = campaignFiles(app, table, await accounts(page));
  const parts = files.map((f) => codeParts(f.code));
  const families = [...new Set(parts.map((p) => p.family))];
  const prefixes = [...new Set(parts.map((p) => p.prefix))];
  // a campaign saved from elsewhere joins the ones the feature put there: a sweep of the prefix would remove them
  if (freshList)
    await readBack(`the ${app} campaigns of ${prefixes.join(', ')}`, () => sweepCampaigns(page, app, families, prefixes));
  atFeatureEnd(page, async () => {
    await readBack(`the ${app} campaigns of ${prefixes.join(', ')}`, () => sweepCampaigns(page, app, families, prefixes));
    await dropPageCopy(page, app);
  });
  await writeCampaigns(page, app, files, freshList);
}

export const campaignsOnServer = Given('the {string} campaigns are on the server:',
  (page: Page, app: string, table: string[][]) => putCampaigns(page, app, table, true),
  {tier: 'api', description: 'campaign folders written through the JS API in the shape the app writes (columns: code, name, status, created, molecules, author me|sharing user, can view / can edit everyone|author, template <name> or "gone: <name>" for a legacy campaign with no snapshot, fields "Owner = Ann; Due = date 2026-11-15", modified = minutes ago, also written to hitdesign.update_logs); the page drops its cached list, as a fresh page would; removed at feature end with every campaign of their prefixes and their hitdesign rows, read back'});

export const campaignsFromElsewhere = When('the {string} campaigns are saved on the server from elsewhere:',
  (page: Page, app: string, table: string[][]) => putCampaigns(page, app, table, false),
  {tier: 'api', description: 'the same files as "… campaigns are on the server", written while the landing is open, as another session would save them, beside the campaigns already there: the page keeps the list it read'});

// ---- the fixtures: templates ------------------------------------------------------------------------------------

const TEMPLATE_COLUMNS = ['name', 'key', 'fields', 'stages', 'submit', 'layout'];
const SUBMIT_FUNCTIONS: {[friendly: string]: {fName: string, package: string}} =
  {'Demo File Submit': {fName: 'demoFileSubmit', package: 'HitTriage'}};
const FIELD_TYPES = ['String', 'Number', 'Boolean', 'Date', 'Molecule'];

/** "Owner: String required; Budget: Number" */
function parseFieldDefs(text: string): {name: string, type: string, required: boolean}[] {
  return text.split(';').map((p) => p.trim()).filter((p) => p !== '').map((p) => {
    const m = /^(.+?):\s*(\w+)(\s+required)?$/.exec(p);
    if (!m || !FIELD_TYPES.includes(m[2]))
      throw new Error(`a template field is "<name>: <${FIELD_TYPES.join('|')}>[ required]", not "${p}"`);
    return {name: m[1].trim(), type: m[2], required: !!m[3]};
  });
}

/** "bar chart of Stage" or "scatter plot of MolWt": the layout a table view saves with that one viewer next to its grid. */
function layoutViewState(page: Page, text: string): Promise<string> {
  const m = /^(bar chart|scatter plot) of (\S+)$/.exec(text);
  if (!m)
    throw new Error(`a fixture template's layout is "bar chart of <column>" or "scatter plot of <column>", not "${text}"`);
  return page.evaluate(async ([viewer, col]) => {
    const extra = ['Molecule', 'Stage', 'V-iD'].includes(col) ? '' : `,${col}`;
    const df = DG.DataFrame.fromCsv(`Molecule,Stage,V-iD${extra}\nCCO,Design,${extra ? ',1.5' : ''}\nc1ccccc1,Design,${extra ? ',2.5' : ''}\n`);
    df.name = 'bdd-layout-fixture';
    const tv = grok.shell.addTableView(df);
    tv.addViewer(viewer === 'bar chart' ? 'Bar chart' : 'Scatter plot',
      viewer === 'bar chart' ? {splitColumnName: col} : {xColumnName: col, yColumnName: col});
    const state = tv.saveLayout().viewState;
    tv.close();
    grok.shell.closeTable(df);
    return state;
  }, [m[1], m[2]] as [string, string]);
}

async function templateJson(page: Page, app: App, r: {[k: string]: string}): Promise<any> {
  keyParts(r['key']);
  templateFamily(r['name']);
  const submit = r['submit'] ? SUBMIT_FUNCTIONS[r['submit']] : undefined;
  if (r['submit'] && !submit)
    throw new Error(`the fixture submit functions are ${Object.keys(SUBMIT_FUNCTIONS).join(', ')}, not "${r['submit']}"`);
  const out: any = {name: r['name'], key: r['key'], campaignFields: parseFieldDefs(r['fields'] ?? '')};
  if (isDesign(app))
    out.stages = (r['stages'] ?? '').split(',').map((s) => s.trim()).filter((s) => s !== '');
  else
    out.dataSourceType = 'File';
  if (app === 'PepTriage')
    out.sequenceColumnName = 'Helm';
  if (r['layout'])
    out.layoutViewState = await layoutViewState(page, r['layout']);
  out.compute = NO_COMPUTE;
  if (submit)
    out.submit = submit;
  return out;
}

async function templatesAtEnd(page: Page, app: App, names: string[], keys: string[]): Promise<void> {
  const families = [...new Set(names.map(templateFamily))];
  const what = `the ${app} templates ${names.join(', ')}`;
  const sweep = async (): Promise<string> => {
    const {left, keys: found} = await sweepTemplates(page, app, families, names);
    const c = campaignPrefixesOf([...keys, ...found]);
    const campaignsLeft = c.prefixes.length > 0 ? await sweepCampaigns(page, app, c.families, c.prefixes) : '';
    return [left, campaignsLeft].filter((s) => s !== '').join(', ');
  };
  await readBack(what, sweep);
  atFeatureEnd(page, async () => {
    await readBack(what, sweep);
    await dropPageCopy(page, app);
  });
}

export const templatesOnServer = Given('the {string} templates are on the server:',
  async (page: Page, appName: string, table: string[][]) => {
    const app = appOf(appName);
    const [header, ...rows] = table;
    const unknown = header.filter((h) => !TEMPLATE_COLUMNS.includes(h));
    if (unknown.length > 0 || !header.includes('name') || !header.includes('key'))
      throw new Error(`the templates table takes the columns ${TEMPLATE_COLUMNS.join(', ')} (name and key required); it has ${unknown.join(', ')}`);
    const specs = rows.map((cells) => {
      const r: {[k: string]: string} = {};
      header.forEach((h, i) => r[h] = cells[i]?.trim() ?? '');
      return r;
    });
    await templatesAtEnd(page, app, specs.map((r) => r['name']), specs.map((r) => r['key']));
    const jsons = [];
    for (const r of specs)
      jsons.push(await templateJson(page, app, r));
    await page.evaluate(async ([root, ts]) => {
      for (const t of ts)
        await grok.dapi.files.writeAsText(`${root}/${t.name}.json`, JSON.stringify(t));
    }, [`${share(app)}/templates`, jsons] as [string, any[]]);
  }, {tier: 'api', description: 'template files written through the JS API as the template editor writes them (columns: name, key, fields "Owner: String required; Due: Date", stages "Design, Make", submit "Demo File Submit", layout "bar chart of Stage" / "scatter plot of MolWt"; no compute); removed at feature end with the campaigns their keys started and those campaigns\' hitdesign rows, read back'});

export const templatesRemoved = Given('the {string} templates named {string} are removed now and at feature end',
  (page: Page, appName: string, names: string) =>
    templatesAtEnd(page, appOf(appName), names.split(',').map((n) => n.trim()).filter((n) => n !== ''), []),
  {tier: 'api', description: 'for the templates a feature makes through the UI (names ending in {time}): their files, and the campaigns their keys started, go now and at feature end'});

/* The descriptors the Compute section lists come from Chem, which on an RDKit stand asks its container: a feature
   never waits on that service, so the editor opens while the list cannot be read. */
export const descriptorsUnlisted = Given('the descriptors cannot be listed, put back at feature end', async (page: Page) => {
  await page.evaluate(() => {
    const chem = grok.chem;
    chem.bddDescriptorsTree ??= chem.descriptorsTree;
    chem.descriptorsTree = () => Promise.reject(new Error('the descriptors are not available'));
  });
  atFeatureEnd(page, () => page.evaluate(() => {
    const chem = grok.chem;
    if (chem.bddDescriptorsTree) {
      chem.descriptorsTree = chem.bddDescriptorsTree;
      delete chem.bddDescriptorsTree;
    }
  }));
}, {tier: 'api', description: 'grok.chem.descriptorsTree answers "the descriptors are not available" in the page, as when Chem\'s descriptor service is down; put back at feature end'});

// ---- the per-viewer settings of the landings --------------------------------------------------------------------

const SETTINGS = /^(HD|HT)Campaign|^Accordion:Hit-Design-campaigns-group-/;

export const settingsDefaults = Given('the campaigns tables start from their defaults, put back at feature end',
  async (page: Page) => {
    const saved: [string, string][] = await page.evaluate((re) => {
      const keys = Object.keys(localStorage).filter((k) => new RegExp(re).test(k));
      const out = keys.map((k) => [k, localStorage.getItem(k) ?? ''] as [string, string]);
      for (const k of keys)
        localStorage.removeItem(k);
      return out;
    }, SETTINGS.source);
    atFeatureEnd(page, () => page.evaluate(([re, back]) => {
      for (const k of Object.keys(localStorage).filter((k) => new RegExp(re).test(k)))
        localStorage.removeItem(k);
      for (const [k, v] of back)
        localStorage.setItem(k, v);
    }, [SETTINGS.source, saved] as [string, [string, string][]]));
  }, {tier: 'api', description: 'the columns, grouping and sorting the landings keep in the browser (HDCampaign*, HTCampaign*) and the groups\' expanded state (Accordion:Hit-Design-campaigns-group-*) are cleared now and put back as they were at feature end'});

export const sortedBy = Given('the campaigns table of {string} is sorted by {string} {word}',
  async (page: Page, appName: string, column: string, direction: string) => {
    const app = appOf(appName);
    if (direction !== 'ascending' && direction !== 'descending')
      throw new Error(`a sort is ascending or descending, not "${direction}"`);
    await page.evaluate(([key, value]) => localStorage.setItem(key, value),
      [isDesign(app) ? 'HDCampaignsTableSorting' : 'HTCampaignsTableSorting',
        JSON.stringify({columnName: column, ascending: direction === 'ascending'})] as [string, string]);
  }, {tier: 'api', description: 'the sort a double-click on the header leaves in the browser (HDCampaignsTableSorting for Hit Design and PeptiHit, HTCampaignsTableSorting for Hit Triage and PepTriage), set before the landing opens; with "the campaigns tables start from their defaults" it is put back at feature end'});

// ---- the campaigns table: reading it ------------------------------------------------------------------------------

interface TableRead {headers: string[], arrows: string[], rows: string[][], names: string[]}

/** The campaigns table of the view shown: its header captions and sort arrows (`aria-sort`), the cells of every row
 * it rendered, and the names of all its campaigns in order (the rows of a folded group, which it does not render,
 * included). */
function readTable(page: Page): Promise<TableRead | null> {
  return page.evaluate((selector) => {
    const table = Array.from(document.querySelectorAll(selector)).find((t) => (t as HTMLElement).offsetParent !== null);
    if (!table)
      return null;
    const heads = Array.from(table.querySelectorAll('th[data-u2="column-header"]'));
    const control = (window as any).DG.U2.Control.forElement(table);
    return {
      headers: heads.map((th) => (th.textContent ?? '').trim()),
      arrows: heads.map((th) => ({ascending: 'ascending', descending: 'descending'} as {[k: string]: string})[th.getAttribute('aria-sort') ?? ''] ?? ''),
      rows: Array.from(table.querySelectorAll('tr.u2-table-row'))
        .map((tr) => Array.from(tr.children).map((td) => (td.textContent ?? '').trim())),
      names: control.getItems({all: true}).list.filter((i: any) => !i.key.startsWith('section:')).map((i: any) => i.label),
    };
  }, TABLE);
}

async function pollTable<T>(page: Page, what: string, read: (t: TableRead) => T, want: T): Promise<void> {
  let last: TableRead | null = null;
  let got: unknown = null;
  await expect.poll(async () => {
    last = await readTable(page);
    return got = last ? read(last) : '(no campaigns table is shown)';
  }, {message: what}).toEqual(want).catch(() => {
    const t = last as TableRead | null;
    throw new Error(`${what}: expected ${JSON.stringify(want)}, read ${JSON.stringify(got)}` +
      (t ? `; the table shows the columns "${t.headers.filter((h) => h).join(', ')}" and ${t.rows.length} row(s): ` +
        t.rows.slice(0, 25).map((r) => r[0]).join(' | ') : ''));
  });
}

export const tableColumns = Then('the campaigns table should have the columns {string}', (page: Page, list: string) =>
  pollTable(page, `the columns of the campaigns table`, (t) => t.headers.filter((h) => h !== '').join(', '),
    list.split(',').map((s) => s.trim()).join(', ')),
{description: 'the header captions in order, the sort arrow and the unlabeled row-action columns apart'});

/** The order the table sorts the shown values in (u2 CatalogTable): numbers when every value reads as one, else dates
 * when every value parses as one, else text in natural order (KHL-2 before KHL-12); empty values go last both ways
 * and are left to `sortProblem`. */
function compareValues(values: string[]): (a: string, b: string) => number {
  const filled = values.filter((v) => v !== '');
  const numbers = filled.every((v) => !isNaN(Number(v)));
  const dates = !numbers && filled.every((v) => !isNaN(Date.parse(v)));
  const text = new Intl.Collator(undefined, {numeric: true});
  return (a, b) => numbers ? Number(a) - Number(b) : dates ? Date.parse(a) - Date.parse(b) : text.compare(a, b);
}

function sortProblem(t: TableRead, column: string, direction: string): string {
  const i = t.headers.indexOf(column);
  if (i < 0)
    return `no "${column}" column`;
  if (t.arrows[i] !== direction)
    return `the "${column}" header shows ${t.arrows[i] ? `the ${t.arrows[i]} arrow` : 'no arrow'}`;
  if (t.rows.length < 2)
    return `${t.rows.length} row(s): no order to see`;
  const values = t.rows.map((r) => r[i] ?? '');
  const cmp = compareValues(values);
  for (let k = 1; k < values.length; k++) {
    const [a, b] = [values[k - 1], values[k]];
    if (b === '')
      continue;
    if (a === '' || (direction === 'ascending' ? cmp(a, b) > 0 : cmp(a, b) < 0))
      return `row ${k} ("${t.rows[k - 1][0]}": "${a}") comes before row ${k + 1} ("${t.rows[k][0]}": "${b}")`;
  }
  return '';
}

export const tableSorted = Then('the campaigns table should be sorted by {string} {word}',
  (page: Page, column: string, direction: string) =>
    pollTable(page, `the campaigns table sorted by ${column} ${direction}`, (t) => sortProblem(t, column, direction), ''),
  {description: 'the header shows the arrow of that direction, and the values of the column go that way down every row the table rendered — as numbers when all are numbers, else as dates, else as text in natural order (KHL-2 before KHL-12); empty values last either way'});

export const tableNotSorted = Then('the campaigns table should not be sorted', (page: Page) =>
  pollTable(page, 'the sort arrows of the campaigns table', (t) => t.arrows.filter((a) => a !== '').length, 0),
{description: 'no header shows a sort arrow'});

function cellOf(t: TableRead, column: string, row: string): string | null {
  const i = t.headers.indexOf(column);
  const r = t.rows.find((cells) => cells[0] === row);
  return i < 0 || !r ? null : r[i] ?? '';
}

export const cellIs = Then('the {string} column of the {string} campaign row should be {string}',
  (page: Page, column: string, row: string, value: string) =>
    pollTable(page, `the ${column} of "${row}" in the campaigns table`, (t) => cellOf(t, column, row), value),
  {description: 'the text of that cell, exactly'});

export const cellContains = Then('the {string} column of the {string} campaign row should contain {string}',
  (page: Page, column: string, row: string, text: string) =>
    pollTable(page, `the ${column} of "${row}" in the campaigns table, containing "${text}"`,
      (t) => (cellOf(t, column, row) ?? '(no such cell)').includes(text) ? text : cellOf(t, column, row) ?? '(no such cell)', text),
  {description: 'the text of that cell holds the text (case-sensitive)'});

const namesOf = (list: string): string[] => list.split(',').map((s) => s.trim()).filter((s) => s !== '');

export const listed = Then('the campaigns table should list {string}', (page: Page, list: string) =>
  pollTable(page, `the campaigns ${list} in the table`,
    (t) => namesOf(list).filter((n) => !t.names.includes(n)).join(', ') || 'all listed', 'all listed'),
{description: 'each of the names (a comma-separated list) is a row of the table — in a folded group too'});

export const notListed = Then('the campaigns table should not list {string}', (page: Page, list: string) =>
  pollTable(page, `the campaigns ${list} not in the table`,
    (t) => namesOf(list).filter((n) => t.names.includes(n)).join(', ') || 'none listed', 'none listed'),
{description: 'none of the names is a row of the table'});

export const groupCount = Then('the {string} campaign group should count {int} campaigns', async (page: Page, label: string, count: number) => {
  let text = '';
  await expect.poll(async () => {
    text = await page.evaluate(([selector, l]) => {
      const header = Array.from(document.querySelectorAll(`${selector} tr.u2-table-section`))
        .find((h) => (h as HTMLElement).offsetParent !== null &&
          (h.querySelector('[data-u2-part="item-label"]')?.textContent ?? '') === l);
      return header ? (header.textContent ?? '').trim() : '(no such group)';
    }, [TABLE, label] as [string, string]);
    return (text.replace(label, '').match(/\d+/g) ?? []).map(Number);
  }, {message: `the count in the header of the "${label}" group`}).toContain(count).catch(() => {
    throw new Error(`the header of the "${label}" group reads "${text}": no count ${count} beside its label`);
  });
}, {description: 'the header of the group shows the number of campaigns it holds beside its label (the label\'s own digits apart)'});

export const listedInOrder = Then('the campaigns table should list {string} in this order', (page: Page, list: string) => {
  const names = list.split(',').map((s) => s.trim());
  return pollTable(page, `the campaigns ${names.join(', ')} in this order`, (t) => {
    const order = t.names;
    const missing = names.filter((n) => !order.includes(n));
    if (missing.length > 0)
      return `not listed: ${missing.join(', ')}`;
    return names.slice().sort((a, b) => order.indexOf(a) - order.indexOf(b)).join(', ');
  }, names.join(', '));
}, {description: 'each of them is a row the table rendered, and they come in that order (other rows may lie between)'});

// ---- the campaigns table: gestures ------------------------------------------------------------------------------

export const doubleClickHeader = When('user double-clicks the {string} header of the campaigns table',
  async (page: Page, column: string) => {
    const header = page.locator(`${TABLE} th[data-u2="column-header"]`).filter({visible: true})
      .filter({hasText: exactText(column)});
    await expect(header, {message: `the "${column}" header of the campaigns table`}).toHaveCount(1);
    await header.dblclick();
  }, {tier: 'ui', description: 'a double click on the header cell of that column (a single click does nothing)'});

// ---- menus ---------------------------------------------------------------------------------------------------------

function menuCheck(page: Page, item: string): Promise<string> {
  return page.evaluate((label) => {
    const items = Array.from(document.querySelectorAll('.d4-menu-popup [role="menuitem"], .u2-menu-item'))
      .filter((i) => (i as HTMLElement).offsetParent !== null &&
        (i.querySelector(':scope > .d4-menu-item-label, :scope > .u2-menu-label')?.textContent ?? '').trim() === label);
    if (items.length !== 1)
      return `${items.length} "${label}" items in the open menu`;
    return items[0].getAttribute('aria-checked') === 'true' ? 'checked' : 'unchecked';
  }, item);
}

export const menuItemChecked = Then('the open menu should show {string} as {word}', async (page: Page, item: string, state: string) => {
  if (state !== 'checked' && state !== 'unchecked')
    throw new Error(`a menu item shows as checked or unchecked, not "${state}"`);
  await expect.poll(() => menuCheck(page, item), {message: `"${item}" in the open menu`}).toBe(state);
}, {description: 'the check mark of a popup menu item (aria-checked) — a Dart menu\'s or a u2 one\'s, a submenu\'s items included'});

/** A confirmation's button says what it does ("DELETE"); a platform dialog's is OK. Either is the dialog's confirming
 * button. */
export const confirmDialog = When('user confirms {element}', async (page: Page, target: ElementRef) => {
  const dialog = (await locate(page, target)).filter({visible: true}).first();
  const button = dialog.locator('button').filter({visible: true}).filter({hasText: /^\s*(OK|DELETE)\s*$/i});
  await expect(button, {message: `the OK or DELETE button of ${target.phrase}`}).toHaveCount(1);
  await button.click();
}, {tier: 'ui', description: 'clicks the button that confirms the dialog: OK, or DELETE where the dialog names what it does'});

export const campaignViewsOpen = Then('{int} campaign view(s) should be open', async (page: Page, count: number) => {
  await expect.poll(() => page.evaluate(() => Array.from(grok.shell.views).filter((v: any) => v.type === 'TableView').length),
    {message: 'the table views open in the shell (a campaign is one)'}).toBe(count);
}, {description: 'the open table views, whatever they are named: an app\'s campaigns open as table views'});

// ---- the server -------------------------------------------------------------------------------------------------

function campaignJson(page: Page, app: App, code: string): Promise<any | null> {
  return page.evaluate(async (p) => await grok.dapi.files.exists(p) ? JSON.parse(await grok.dapi.files.readAsText(p)) : null,
    `${share(app)}/campaigns/${code}/campaign.json`);
}

export const campaignOnServerClaim = Then('the {string} campaign {string} should be on the server',
  async (page: Page, appName: string, code: string) => {
    const app = appOf(appName);
    await expect.poll(async () => (await campaignJson(page, app, code)) !== null,
      {message: `${app} campaign ${code} on the server`, timeout: pollMs(20000)}).toBe(true);
  }, {tier: 'api', description: 'its campaign.json is on the share'});

export const campaignNotOnServer = Then('the {string} campaign {string} should not be on the server',
  async (page: Page, appName: string, code: string) => {
    const app = appOf(appName);
    await expect.poll(() => page.evaluate((p) => grok.dapi.files.exists(p), `${share(app)}/campaigns/${code}`),
      {message: `${app} campaign ${code} gone from the server`, timeout: pollMs(20000)}).toBe(false);
  }, {tier: 'api', description: 'its folder is gone from the share'});

const CAMPAIGN_KEYS: {[word: string]: (c: any) => string} = {
  name: (c) => String(c.friendlyName ?? ''), status: (c) => String(c.status ?? ''),
  author: (c) => String(c.authorUserFriendlyName ?? ''), template: (c) => String(c.templateName ?? ''),
};

export const campaignHas = Then('the {string} campaign {string} on the server should have the {word} {string}',
  async (page: Page, appName: string, code: string, word: string, value: string) => {
    const app = appOf(appName);
    const read = CAMPAIGN_KEYS[word];
    if (!read)
      throw new Error(`a campaign's ${word} is not claimed; the claims are ${Object.keys(CAMPAIGN_KEYS).join(', ')}`);
    await expect.poll(async () => {
      const c = await campaignJson(page, app, code);
      return c ? read(c) : '(no such campaign)';
    }, {message: `the ${word} of ${app} campaign ${code} on the server`, timeout: pollMs(20000)}).toBe(value);
  }, {tier: 'api', description: 'read from its campaign.json: name (the friendly name), status, author (friendly name), template (its name)'});

function templateJsonOnServer(page: Page, app: App, name: string): Promise<any | null> {
  return page.evaluate(async (p) => await grok.dapi.files.exists(p) ? JSON.parse(await grok.dapi.files.readAsText(p)) : null,
    `${share(app)}/templates/${name}.json`);
}

export const templateOnServer = Then('the {string} template {string} should be on the server',
  async (page: Page, appName: string, name: string) => {
    const app = appOf(appName);
    await expect.poll(async () => (await templateJsonOnServer(page, app, name)) !== null,
      {message: `${app} template ${name} on the server`, timeout: pollMs(20000)}).toBe(true);
  }, {tier: 'api', description: 'its <name>.json is in the app\'s templates folder'});

export const templateNotOnServer = Then('the {string} template {string} should not be on the server',
  async (page: Page, appName: string, name: string) => {
    const app = appOf(appName);
    await expect.poll(async () => (await templateJsonOnServer(page, app, name)) !== null,
      {message: `${app} template ${name} gone from the server`, timeout: pollMs(20000)}).toBe(false);
  }, {tier: 'api', description: 'no <name>.json in the app\'s templates folder'});

const TEMPLATE_KEYS: {[word: string]: (t: any) => string} = {
  key: (t) => String(t.key ?? ''),
  fields: (t) => (t.campaignFields ?? []).map((f: any) => `${f.name} (${f.type}${f.required ? ', required' : ''})`).join(', '),
  stages: (t) => (t.stages ?? []).join(', '),
  submit: (t) => t.submit ? `${t.submit.package}:${t.submit.fName}` : '',
  layout: (t) => t.layoutViewState ? 'a layout' : 'no layout',
};

export const templateHas = Then('the {string} template {string} on the server should have the {word} {string}',
  async (page: Page, appName: string, name: string, word: string, value: string) => {
    const app = appOf(appName);
    const read = TEMPLATE_KEYS[word];
    if (!read)
      throw new Error(`a template's ${word} is not claimed; the claims are ${Object.keys(TEMPLATE_KEYS).join(', ')}`);
    await expect.poll(async () => {
      const t = await templateJsonOnServer(page, app, name);
      return t ? read(t) : '(no such template)';
    }, {message: `the ${word} of ${app} template ${name} on the server`, timeout: pollMs(20000)}).toBe(value);
  }, {tier: 'api', description: 'read from its JSON: key; fields "Owner (String, required), Due (Date)"; stages "Design, Make"; submit "HitTriage:demoFileSubmit"; layout "a layout" or "no layout"'});

// ---- the template editor ------------------------------------------------------------------------------------------

/** The two RecordsInputs of the Hit Design / PeptiHit editor: the Additional fields (Name, Type, Required per row) and
 * the Stages (Name). "Add field" / "Add stage" under the rows adds an empty row; a row typed into is a row, saved
 * without leaving it. */
const FIELDS = `${TEMPLATE_EDITOR} [data-u2-name="fields"]`;
const STAGES = `${TEMPLATE_EDITOR} [data-u2-name="stages"]`;

async function newRow(page: Page, records: string) {
  await page.locator(`${records} [data-u2-part="add"]`).click();
  const row = page.locator(`${records} .u2-records-body > .u2-records-row`).last();
  await expect(row, {message: `the new row of the ${records === FIELDS ? 'Additional fields' : 'Stages'}`}).toBeVisible();
  return row;
}

async function typeField(page: Page, name: string, type: string, required: boolean, leave: boolean): Promise<void> {
  const row = await newRow(page, FIELDS);
  const nameEditor = row.locator('[data-u2="text-input"] input');
  await gestures.typeVerified(nameEditor, name, 'the name of the new field');
  await row.locator('select').selectOption(type);
  const box = row.locator('input[type="checkbox"]');
  if (await box.isChecked() !== required)
    await box.click();
  await expect(box, {message: 'the Required box of the new field'}).toBeChecked({checked: required});
  if (leave)
    await box.press('Tab');
}

export const addField = When('user adds the field {string} of type {string} in the template editor',
  (page: Page, name: string, type: string) => typeField(page, name, type, false, true),
  {tier: 'ui', description: 'clicks "Add field", fills the new row of Additional fields and leaves it'});
export const addRequiredField = When('user adds the required field {string} of type {string} in the template editor',
  (page: Page, name: string, type: string) => typeField(page, name, type, true, true),
  {tier: 'ui', description: 'clicks "Add field", fills the new row of Additional fields, ticks Required and leaves the row'});
export const typePendingField = When('user types the field {string} of type {string} in the template editor without adding it',
  (page: Page, name: string, type: string) => typeField(page, name, type, false, false),
  {tier: 'ui', description: 'clicks "Add field" and fills the new row of Additional fields, the focus left in it (CREATE takes it)'});

async function typeStage(page: Page, name: string, leave: boolean): Promise<void> {
  const editor = (await newRow(page, STAGES)).locator('[data-u2="text-input"] input');
  await gestures.typeVerified(editor, name, 'the name of the new stage');
  if (leave)
    await editor.press('Tab');
}

export const addStage = When('user adds the stage {string} in the template editor', (page: Page, name: string) => typeStage(page, name, true),
  {tier: 'ui', description: 'clicks "Add stage", types the name into the new row of Stages and leaves it'});
export const typePendingStage = When('user types the stage {string} in the template editor without adding it',
  (page: Page, name: string) => typeStage(page, name, false),
  {tier: 'ui', description: 'clicks "Add stage" and types the name into the new row of Stages, the focus left in it'});

/** The rows of a RecordsInput as its value holds them, the ones with no name apart (CREATE leaves them out). */
function readRecords(page: Page, records: string): Promise<string> {
  return page.evaluate((sel) => {
    const el = Array.from(document.querySelectorAll(sel)).find((e) => (e as HTMLElement).offsetParent !== null);
    if (!el)
      return '(not shown)';
    const rows: {name: string, type?: string, required?: boolean}[] = DG.U2.control(el).peek() ?? [];
    return rows.filter((r) => r.name).map((r) => r.type === undefined ? r.name :
      `${r.name} (${r.type}${r.required ? ', required' : ''})`).join(', ');
  }, records);
}

export const editorFields = Then('the template editor should list the fields {string}', async (page: Page, list: string) => {
  await expect.poll(() => readRecords(page, FIELDS), {message: 'the rows of Additional fields'}).toBe(list);
}, {description: 'the named rows of Additional fields as "Owner (String, required), Due (Date)"'});

export const editorStages = Then('the template editor should list the stages {string}', async (page: Page, list: string) => {
  await expect.poll(() => readRecords(page, STAGES), {message: 'the rows of Stages'}).toBe(list);
}, {description: 'the named rows of Stages as "Design, Make"'});

/** The messages under the editor's inputs: u2's error lines (Hit Design, PeptiHit), and the Hit Triage / PepTriage
 * error divs, always in the DOM and shown by their opacity, which Playwright's visibility ignores. */
function shownErrors(page: Page): Promise<string[]> {
  return page.evaluate(() => {
    const root = grok.shell.v?.root as HTMLElement | undefined;
    const legacy = Array.from(root?.querySelectorAll('.hit-triage-error-div') ?? [])
      .filter((e) => (e as HTMLElement).offsetParent !== null && Number(getComputedStyle(e).opacity) > 0);
    const u2 = Array.from(root?.querySelectorAll('[data-u2-name="template-editor"] [data-u2-part="error"]') ?? [])
      .filter((e) => (e as HTMLElement).offsetParent !== null);
    return [...legacy, ...u2].map((e) => (e.textContent ?? '').trim()).filter((t) => t !== '');
  });
}

export const editorShowsError = Then('the template editor should show the error {string}', async (page: Page, text: string) => {
  await expect.poll(async () => (await shownErrors(page)).includes(text), {message: `the error "${text}" under the template editor's inputs`}).toBe(true);
}, {description: 'an error line of the editor\'s Details shown (not transparent)'});

export const editorShowsNoError = Then('the template editor should not show the error {string}', async (page: Page, text: string) => {
  await expect.poll(async () => (await shownErrors(page)).includes(text), {message: `the error "${text}" under the template editor's inputs`}).toBe(false);
}, {description: 'no error line of the editor\'s Details shows that text'});

export const viewNotOpen = Then('the {string} view should not be open', async (page: Page, name: string) => {
  await expect.poll(() => page.evaluate((n) => Array.from(grok.shell.views).some((v: any) => v.name === n), name),
    {message: `a view named "${name}" open`}).toBe(false);
}, {description: 'no open view has that name (closed, not only behind another)'});
