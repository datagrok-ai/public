/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {waitUntil} from '@datagrok-libraries/test/src/test';
import {_package} from '../package';
import {AppName, HitDesignTemplate, HitTriageTemplate, PepTriageTemplate} from '../app/types';
import {CampaignJsonName, CampaignTableName, HitDesignMolColName, PeptiHitHelmColName} from '../app/consts';
import {HitDesignApp} from '../app/hit-design-app';
import {HitTriageApp} from '../app/hit-triage-app';
import {PeptiHitApp} from '../app/pepti-hit-app';
import {PepTriageApp} from '../app/pep-triage-app';
import {obfuscateSmiles} from '../app/utils/molreg';
import type {INewCampaignResult} from '../app/accordeons/new-campaign-accordeon';
import {defaultCampaignStorage, getDefaultCampaignStorageSettings} from '../packageSettingsEditor';
import * as api from '../package-api';

export const allUsers = DG.Group.defaultGroupsIds['All users'];
export const otherUserId = '00000000-0000-0000-0000-000000000001';
export const unknownGroupId = '00000000-0000-0000-0000-00000000b3b3';

/** JSON with sorted keys, so two objects compare by content. */
function canon(x: unknown): string {
  const sort = (v: any): any => Array.isArray(v) ? v.map(sort) : v && typeof v === 'object' ?
    Object.fromEntries(Object.keys(v).sort().map((k) => [k, sort(v[k])])) : v;
  return JSON.stringify(sort(x));
}

export function expectSame(actual: unknown, expected: unknown, what: string): void {
  if (canon(actual) !== canon(expected))
    throw new Error(`${what}: expected ${canon(expected)}, got ${canon(actual)}`);
}

/** Ad-hoc SQL on the database the package's own queries use. */
export async function hitdesignSql(sql: string): Promise<DG.DataFrame> {
  const connection = (DG.Func.find({package: 'HitTriage', name: 'GetLastModified'})[0] as DG.DataQuery).connection;
  return await connection.query('ZZB3 probe', sql).apply();
}

export const sqlText = (s: string) => `'${s.replace(/'/g, '\'\'')}'`;

/** Records the calls of a package query made while [body] runs (live, in [calls]), each as its [input]'s length. */
export async function countQueryCalls<T>(queryName: string, input: string | null,
  body: (calls: number[]) => Promise<T>): Promise<{result: T, calls: number[]}> {
  const calls: number[] = [];
  const sub = grok.functions.onBeforeRunAction.subscribe((fc) => {
    if (fc?.func?.nqName?.toLowerCase() === `hittriage:${queryName.toLowerCase()}`)
      calls.push(input ? ((fc.inputs as {[name: string]: unknown[] | undefined})[input]?.length ?? -1) : 1);
  });
  try {
    return {result: await body(calls), calls};
  } finally {
    sub.unsubscribe();
  }
}

export function noComputeTemplate(name: string, key: string): HitDesignTemplate {
  return {
    name, key,
    campaignFields: [
      {name: 'Chemist', type: 'String', required: true},
      {name: 'Deadline', type: 'Date'},
      {name: 'Target MW', type: 'Number'},
    ],
    stages: ['Idea', 'Made', 'Tested'],
    compute: {descriptors: {enabled: false, args: []}, functions: [], scripts: [], queries: []},
    submit: {fName: 'demoFileSubmit', package: 'HitTriage'},
  };
}

export function triageTemplate(name: string, key: string): HitTriageTemplate {
  return {
    name, key,
    campaignFields: [
      {name: 'Chemist', type: 'String', required: true},
      {name: 'Deadline', type: 'Date', required: false},
      {name: 'Scaffold', type: 'Molecule', required: false},
    ],
    dataSourceType: 'Query',
    queryFunctionName: 'Demo Molecules variable',
    compute: {descriptors: {enabled: false, args: []}, functions: [], scripts: [], queries: []},
    submit: {fName: 'demoFileSubmit', package: 'HitTriage'},
  };
}

export function pepTriageTemplate(name: string, key: string): PepTriageTemplate {
  return {...triageTemplate(name, key), sequenceColumnName: 'Helm', moleculeColumnName: 'molfile(Helm)',
    queryFunctionName: 'Demo Peptide Sequences'};
}

/** A Date campaign field as the new-campaign form writes it. */
export const formDate = {_value: '2026-03-27T00:00:00.000Z', isUtc: false, date: '2026-03-27T00:00:00.000Z'};
/** The older shape of the demo campaigns (HTD-1). */
export const oldDate = {a: 1759953600000, b: false, date: '2025-10-08T20:00:00.000Z'};

export const designCsv = 'Molecule,Stage,V-iD,Note\nCCO,Idea,,first\nc1ccccc1,Made,,second\n';
export const designColumns = {
  columnSemTypes: {'Molecule': 'Molecule', 'Stage': null, 'V-iD': null, 'Note': null},
  columnTypes: {'Molecule': 'string', 'Stage': 'string', 'V-iD': 'string', 'Note': 'string'},
};

/** A Hit Design campaign.json in its full shape, over {@link designCsv}. */
export function designCampaign(sb: Sandbox, id: string, template: HitDesignTemplate, extra: {[key: string]: unknown} = {}):
  {name: string, [key: string]: any} {
  const user = DG.User.current();
  return {
    name: id, friendlyName: `Friendly ${id}`, templateName: template.name, status: 'In Progress', createDate: '2024/9/4',
    campaignFields: {'Chemist': 'Ann', 'Deadline': formDate, 'Target MW': 350},
    ...designColumns, rowCount: 2, filteredRowCount: 2, savePath: sb.savePath(id), version: 3,
    authorUserId: user.id, authorUserFriendlyName: user.friendlyName, lastModifiedUserName: 'Someone Else',
    permissions: {view: [allUsers], edit: [allUsers]}, template, ...extra,
  };
}

/** The landing's Start (`startNewCampaign`): a new campaign from [template], saved once. The app's first row comes
 * filled, so no sketcher opens. */
export async function startDesignCampaign(app: HitDesignApp, template: HitDesignTemplate, fields: {[key: string]: unknown},
  name?: string, first = 'CCO'): Promise<string> {
  const stubbed = app as unknown as {newCampaignFrame: (t: HitDesignTemplate) => DG.DataFrame};
  const frame = stubbed.newCampaignFrame.bind(app);
  stubbed.newCampaignFrame = (t) => {
    const df = frame(t);
    df.col(app instanceof PeptiHitApp ? PeptiHitHelmColName : HitDesignMolColName)!.set(0, first);
    return df;
  };
  try {
    await app.startNewCampaign(template, fields, name);
  } finally {
    delete (stubbed as {newCampaignFrame?: unknown}).newCampaignFrame;
  }
  return app.campaignId!;
}

/** Hit Triage's Start (the landing's `startCampaign`): a new campaign from [template] over a one-molecule table read
 * from a file, saved once. */
export async function startTriageCampaign(app: HitTriageApp, template: HitTriageTemplate, smiles: string,
  name?: string): Promise<string> {
  const start = app.infoView as unknown as {startCampaign: (t: HitTriageTemplate, camp: INewCampaignResult) => Promise<void>};
  // no trailing newline: fromCsv reads one as an empty row
  await start.startCampaign(template, {df: DG.DataFrame.fromCsv(`smiles\n${smiles}`), type: 'File',
    campaignProps: {Chemist: name ?? smiles}, friendlyName: name});
  return app.campaignId!;
}

/** Whether a balloon showing [text] is on the page. */
export function hasBalloon(text: string): boolean {
  return Array.from(document.querySelectorAll('.d4-balloon')).some((b) => b.textContent?.includes(text));
}

/** Appends rows given by column name; filled molecules keep the row-added handler from opening the sketcher. */
export function addRows(df: DG.DataFrame, rows: {[col: string]: unknown}[]): void {
  for (const r of rows)
    df.rows.addNew(df.columns.names().map((n) => r[n] ?? null), true);
}

let nextIsolatedName: AppName = 'PeptiHit';

/** PeptiHit and PepTriage hardcode their app name; the sandbox gives them its own. */
class IsolatedPeptiHitApp extends PeptiHitApp {
  private isolatedName?: AppName;
  constructor(c: DG.FuncCall, name: AppName) {
    nextIsolatedName = name;
    super(c);
    this.isolatedName = name;
  }
  override get appName(): AppName {return this.isolatedName ?? nextIsolatedName;}
}

class IsolatedPepTriageApp extends PepTriageApp {
  private isolatedName?: AppName;
  constructor(c: DG.FuncCall, name: AppName) {
    nextIsolatedName = name;
    super(c);
    this.isolatedName = name;
  }
  override get appName(): AppName {return this.isolatedName ?? nextIsolatedName;}
}

type AnyApp = HitDesignApp | HitTriageApp;

/** The private parts of the Hit Triage app a test has to reach. */
export const internals = (app: AnyApp) => app as unknown as {
  _pickView?: DG.TableView, _pickViewPromise?: Promise<void> | null,
};

/**
 * A throwaway app name. Every path the apps build (`<app>/campaigns`, `<app>/templates`,
 * `System.AppData/HitTriage/<app>/…`) and every hitdesign row (`app_name`) are private to it, so other
 * sessions on the stand never list its campaigns, and the cleanup is one folder and one `app_name`.
 * The few branches on the app name itself (the AI description, `createHitCampaign`'s PeptiHit table) see
 * neither app, so no test goes through them.
 */
export class Sandbox {
  readonly appName: AppName;
  private readonly href = window.location.href;
  private readonly showHelp = grok.shell.windows.showHelp;
  private readonly apps: AnyApp[] = [];
  private maxVidId = 0;
  /** Obfuscated canonical SMILES only this test registers: their dictionary rows go too. */
  readonly ownMolecules = new Set<string>();

  private constructor(tag: string) {
    this.appName = `ZZB3 ${tag} ${DG.Utils.randomString(6)}` as AppName;
  }

  static async create(tag: string): Promise<Sandbox> {
    const folder = await getDefaultCampaignStorageSettings();
    if (folder.replace(/\/$/, '') !== defaultCampaignStorage)
      throw new Error(`new campaigns would go to ${folder}, outside the sandbox: the tests need the default campaign folder`);
    const sb = new Sandbox(tag);
    sb.maxVidId = (await hitdesignSql('SELECT COALESCE(MAX(id), 0) AS id FROM hitdesign.vid_dictionary')).get('id', 0);
    return sb;
  }

  path(id: string, file: string = CampaignJsonName): string {
    return `${this.appName}/campaigns/${id}/${file}`;
  }

  /** The dot form campaign.json keeps in `savePath`. */
  savePath(id: string, folder = 'campaigns'): string {
    return `System.AppData/HitTriage/${this.appName}/${folder}/${id}/${CampaignTableName}`;
  }

  /** The colon form Hit Triage keeps in `ingest.query`. */
  ingestPath(id: string): string {
    return `System:AppData/HitTriage/${this.appName}/campaigns/${id}/${CampaignTableName}`;
  }

  ownMolecule(smiles: string): string {
    this.ownMolecules.add(obfuscateSmiles(_package.convertToSmiles(smiles)));
    return smiles;
  }

  async writeTemplate(t: {name: string}): Promise<void> {
    await _package.files.writeAsText(`${this.appName}/templates/${t.name}.json`, JSON.stringify(t));
  }

  async writeCampaign(campaign: {name: string}, csv: string, tablePath?: string): Promise<void> {
    if (tablePath)
      await grok.dapi.files.writeAsText(tablePath, csv);
    else
      await _package.files.writeAsText(this.path(campaign.name, CampaignTableName), csv);
    await _package.files.writeAsText(this.path(campaign.name), JSON.stringify(campaign));
  }

  async readCampaignText(id: string): Promise<string> {
    return await _package.files.readAsText(this.path(id));
  }

  async readCampaign(id: string): Promise<{[key: string]: any}> {
    return JSON.parse(await this.readCampaignText(id));
  }

  async lastModified(id: string): Promise<string | null> {
    const df = await api.queries.getLastModified(this.appName, id);
    return df.rowCount ? df.col('updated_at')!.getString(0) : null;
  }

  private call(name: string): DG.FuncCall {
    return DG.Func.find({package: 'HitTriage', name})[0].prepare();
  }

  async hitDesign(): Promise<HitDesignApp> {
    return await this.start(new HitDesignApp(this.call('hitDesignApp'), this.appName));
  }

  async peptiHit(): Promise<PeptiHitApp> {
    await grok.functions.call('Bio:initBio', {});
    return await this.start(new IsolatedPeptiHitApp(this.call('peptiHitApp'), this.appName));
  }

  async hitTriage(): Promise<HitTriageApp> {
    return await this.start(new HitTriageApp(this.call('hitTriageApp'), this.appName));
  }

  async pepTriage(): Promise<PepTriageApp> {
    await grok.functions.call('Bio:initBio', {});
    return await this.start(new IsolatedPepTriageApp(this.call('pepTriageApp'), this.appName));
  }

  private async start<T extends AnyApp>(app: T): Promise<T> {
    this.apps.push(app);
    grok.shell.addView(app.multiView);
    const root = app.infoView.root;
    // the template picker, or its offer to create the first template
    await waitUntil(() => root.querySelector('[data-u2-name="template-picker"] select, [data-u2-name="template-picker"] ' +
      '[data-u2="empty-state"]') != null && root.querySelector('.d4-update-shadow') == null,
    `the ${this.appName} landing did not finish loading`, 30000);
    return app;
  }

  /** Opens a Hit Design / PeptiHit campaign the way its link on the landing does. */
  async openDesign(app: HitDesignApp, id: string): Promise<DG.TableView> {
    await app.infoView.setCampaign(id);
    const view = app.designView;
    if (app.campaignId !== id || !view)
      throw new Error(`campaign ${id} did not open`);
    return view;
  }

  /** Opens a Hit Triage / PepTriage campaign; its pick view finishes setting up on timers. */
  async openTriage(app: HitTriageApp, id: string): Promise<DG.TableView> {
    await app.infoView.setCampaign(id);
    await waitUntil(() => app.campaignId === id && internals(app)._pickViewPromise != null, `campaign ${id} did not open`, 30000);
    await internals(app)._pickViewPromise;
    const preset = app.campaign?.filters;
    await waitUntil(() => app.filterSettings !== preset, `the pick view of ${id} never read its filters`, 10000);
    return internals(app)._pickView!;
  }

  async savesDone(): Promise<void> {
    await waitUntil(() => this.apps.every((a) => !(a instanceof HitDesignApp && a.isSaving)),
      'a campaign save is still running', 60000);
  }

  async dispose(): Promise<void> {
    try {
      for (const d of DG.Dialog.getOpenDialogs())
        d.close();
      await this.savesDone();
      const open = Array.from(grok.shell.views).map((v) => v.dart);
      for (const app of this.apps) {
        for (const v of [app instanceof HitDesignApp ? app.designView : undefined, internals(app)._pickView, app.multiView]) {
          if (v && open.includes(v.dart))
            v.close();
        }
      }
      // a progress tracker left open saves its moves as its view closes
      await this.savesDone();
    } finally {
      history.replaceState(history.state, '', this.href);
      grok.shell.windows.showHelp = this.showHelp;
      delete (_package.campaignsCache as {[key: string]: unknown})[this.appName];
      await this.removeServerState();
    }
  }

  private async removeServerState(): Promise<void> {
    if (await _package.files.exists(this.appName))
      await _package.files.delete(this.appName);
    const app = sqlText(this.appName);
    const vids = (await hitdesignSql(`WITH d AS (DELETE FROM hitdesign.campaign_vids WHERE app_name = ${app} RETURNING vid) SELECT DISTINCT vid FROM d`))
      .col('vid')!.toList() as string[];
    if (vids.length) {
      const own = Array.from(this.ownMolecules).map(sqlText).join(',') || 'NULL';
      await hitdesignSql(`WITH d AS (DELETE FROM hitdesign.vid_dictionary v WHERE v.vid IN (${vids.map(sqlText).join(',')}) ` +
        `AND (v.id > ${this.maxVidId} OR v.mh_string IN (${own})) ` +
        'AND NOT EXISTS (SELECT 1 FROM hitdesign.campaign_vids c WHERE c.vid = v.vid) RETURNING 1) SELECT COUNT(*) AS n FROM d');
    }
    await hitdesignSql(`WITH a AS (DELETE FROM hitdesign.update_logs WHERE app_name = ${app} RETURNING 1), ` +
      `b AS (DELETE FROM hitdesign.campaign_locks WHERE app_name = ${app} RETURNING 1) ` +
      'SELECT (SELECT COUNT(*) FROM a) + (SELECT COUNT(*) FROM b) AS n');
    const left = await hitdesignSql(`SELECT (SELECT COUNT(*) FROM hitdesign.update_logs WHERE app_name = ${app}) + ` +
      `(SELECT COUNT(*) FROM hitdesign.campaign_locks WHERE app_name = ${app}) + ` +
      `(SELECT COUNT(*) FROM hitdesign.campaign_vids WHERE app_name = ${app}) AS n`);
    if (await _package.files.exists(this.appName) || left.get('n', 0) > 0)
      throw new Error(`${this.appName}: the cleanup left files or hitdesign rows behind`);
  }
}
