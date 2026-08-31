/* eslint-disable max-len */
/**
 * AI view functions for the Hit apps (Hit Triage, Hit Design, PeptiHit, PepTriage).
 *
 * The vocabulary is registered ONCE per session through `grok.functions.register`
 * (real platform Funcs, no package.ts entries) and reported through `getFunctions()`
 * overrides: the info (landing) views return the info set, and the campaign table
 * views get an instance `getFunctions` override returning only the campaign set plus
 * a bridge (searchHitTableFunctions / callHitTableFunction) to the table view's own
 * standard functions. Every function takes the generic `view` argument the AI
 * assistant injects and resolves the owning app from it at call time.
 */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {HitAppBase} from '../hit-app-base';
import {HitDesignApp} from '../hit-design-app';
import {_package} from '../../package';
import {AppName, HitDesignCampaign, HitDesignTemplate, HitTriageCampaign, HitTriageTemplate} from '../types';
import {TileCategoriesColName, ViDColName} from '../consts';
import {checkEditPermissions, loadTemplate} from '../utils';
import {registerMol} from '../utils/molreg';
import {createNewHitDesignDataFrame} from '../accordeons/new-hit-design-campaign-accordeon';

type AnyCampaign = HitTriageCampaign | HitDesignCampaign;
type AnyTemplate = HitTriageTemplate | HitDesignTemplate;

/** The members of {@link HitDesignApp} / HitTriageApp the functions below rely on
 * (they live on the concrete apps, not on {@link HitAppBase}). */
interface IHitApp {
  campaign?: AnyCampaign;
  campaignId?: string;
  molColName?: string;
  infoView: IHitInfoView;
  saveCampaign(...args: unknown[]): Promise<AnyCampaign>;
  submitParams?: {fName: string, package: string};
  submitView?: {submit(): Promise<boolean>, submitWithStatus?(status?: string): Promise<boolean>};
}

interface IHitInfoView {
  setCampaign(name: string): Promise<unknown>;
  deleteCampaignAndRefresh(appName: AppName, campaignId: string): Promise<void>;
}

type HitApp = HitAppBase<AnyTemplate> & IHitApp;

const AI_NAMESPACE = 'HitTriageAI';
const DEFAULT_SEARCH_LIMIT = 10;

const _funcs = new Map<string, DG.Func>();

/** Registers one vocabulary Func on first use; later calls return the cached one. */
function aiFunc(name: string, signature: string, description: string, run: (...args: any[]) => any): DG.Func {
  let f = _funcs.get(name);
  if (!f) {
    f = grok.functions.register({signature, run, isAsync: true, namespace: AI_NAMESPACE,
      tags: 'hit-ai-function', options: {description}});
    _funcs.set(name, f);
  }
  return f;
}

/** Campaign table views carry a back reference to their app. */
const appByView = new WeakMap<object, HitAppBase<AnyTemplate>>();

/** Resolves the owning app from the AI-injected view argument: the campaign table
 * view (back reference), an info/submit view (its `.app`), or the app's MultiView
 * (its current view's `.app`). */
function appOf(view: any): HitApp {
  const v = view?.jsView ?? view;
  for (const c of [v, view, v?.currentView]) {
    if (c == null)
      continue;
    const byRef = appByView.get(c);
    if (byRef instanceof HitAppBase)
      return byRef as HitApp;
    if (c.app instanceof HitAppBase)
      return c.app as HitApp;
    if (c.currentView?.app instanceof HitAppBase)
      return c.currentView.app as HitApp;
  }
  throw new Error('The current view does not belong to a Hit Triage / Hit Design app');
}

function designAppOf(view: any): HitDesignApp & IHitApp {
  const app = appOf(view);
  if (!(app instanceof HitDesignApp))
    throw new Error('This action applies only to Hit Design / PeptiHit campaigns');
  return app as HitDesignApp & IHitApp;
}

/** The open campaign table view: the injected view itself, or the app's current one. */
function campaignTableView(view: any): DG.TableView {
  const v = view?.jsView ?? view;
  if (v instanceof DG.TableView)
    return v;
  const app = appOf(view) as any;
  const tv = app._designView ?? app._pickView;
  if (tv instanceof DG.TableView)
    return tv;
  throw new Error('No campaign table view is open');
}

function campaignBrief(c: AnyCampaign) {
  return {
    code: c.name,
    ...(c.friendlyName && c.friendlyName !== c.name ? {name: c.friendlyName} : {}),
    template: c.templateName,
    status: c.status ?? null,
    ...(c.authorUserFriendlyName ? {author: c.authorUserFriendlyName} : {}),
    ...(c.createDate ? {created: c.createDate} : {}),
    ...(c.rowCount != null ? {rows: c.rowCount} : {}),
  };
}

function computeBrief(t?: AnyTemplate) {
  const compute = t?.compute;
  if (!compute)
    return null;
  return {
    descriptors: compute.descriptors?.enabled ? compute.descriptors.args : [],
    functions: (compute.functions ?? []).map((f) => `${f.package}:${f.name}`),
    scripts: (compute.scripts ?? []).map((s) => s.name),
    queries: (compute.queries ?? []).map((q) => q.name),
  };
}

async function canEditCampaign(c: AnyCampaign): Promise<boolean> {
  return !c.authorUserId || !c.permissions || await checkEditPermissions(c.authorUserId, c.permissions);
}

async function findCampaign(app: HitApp, campaignCode: string): Promise<AnyCampaign | null> {
  const campaigns = await _package.loadCampaigns(app.appName, []);
  return campaigns[campaignCode] ??
    Object.values(campaigns).find((c) => (c.friendlyName ?? '').toLowerCase() === campaignCode.toLowerCase()) ?? null;
}

async function listTemplateNames(appName: AppName): Promise<string[]> {
  return (await _package.files.list(`${appName}/templates`))
    .filter((f) => f.name.endsWith('.json'))
    .map((f) => f.name.slice(0, -5));
}

/** Reduce a table-function invocation result to something JSON-serializable. */
function serializeResult(value: any): any {
  if (value == null || typeof value === 'string' || typeof value === 'number' || typeof value === 'boolean')
    return value;
  if (value instanceof DG.DataFrame)
    return {type: 'dataframe', name: value.name, rowCount: value.rowCount, columns: value.columns.names()};
  if (value instanceof DG.Column)
    return {type: 'column', name: value.name, length: value.length};
  if (value instanceof DG.ViewBase)
    return {type: 'view', name: value.name};
  try {
    return JSON.parse(JSON.stringify(value));
  } catch (_) {
    return String(value);
  }
}

// ─────────────────────── the info (landing) view set ───────────────────────

export function hitInfoViewFunctions(): DG.Func[] {
  return [
    aiFunc('getHitAppInfo', 'dynamic getHitAppInfo(view view)',
      'Overview of this Hit app - its name, saved campaign and template counts, the campaign currently open. Call this first',
      async (view: any) => {
        const app = appOf(view);
        const campaigns = await _package.loadCampaigns(app.appName, []);
        return {
          app: app.appName,
          campaigns: Object.keys(campaigns).length,
          templates: await listTemplateNames(app.appName),
          currentCampaign: app.campaignId ?? app.campaign?.name ?? null,
          canCreateCampaignsHere: app instanceof HitDesignApp,
        };
      }),

    aiFunc('listHitCampaigns', 'dynamic listHitCampaigns(view view)',
      'List the saved campaigns of this Hit app with code, name, template, status, author and row count. Campaigns are addressed by code',
      async (view: any) => {
        const app = appOf(view);
        const campaigns = Object.values(await _package.loadCampaigns(app.appName, []));
        return {total: campaigns.length, campaigns: campaigns.map(campaignBrief)};
      }),

    aiFunc('getHitCampaignDetails', 'dynamic getHitCampaignDetails(view view, string campaignCode)',
      'Full detail of one saved campaign - fields, row counts, compute methods, stages, whether the current user can edit it',
      async (view: any, campaignCode: string) => {
        const app = appOf(view);
        const c = await findCampaign(app, campaignCode);
        if (!c)
          return {success: false, error: `No campaign '${campaignCode}' — call listHitCampaigns for codes`};
        return {
          ...campaignBrief(c),
          ...(c.campaignFields && Object.keys(c.campaignFields).length ? {campaignFields: c.campaignFields} : {}),
          ...(c.filteredRowCount != null ? {filteredRows: c.filteredRowCount} : {}),
          ...(c.lastModifiedUserName ? {lastModifiedBy: c.lastModifiedUserName} : {}),
          ...((c.template as AnyTemplate & {stages?: string[]})?.stages?.length ?
            {stages: (c.template as AnyTemplate & {stages?: string[]}).stages} : {}),
          compute: computeBrief(c.template),
          canEdit: await canEditCampaign(c),
        };
      }),

    aiFunc('openHitCampaign', 'dynamic openHitCampaign(view view, string campaignCode)',
      'Open a saved campaign so the user can work on it - loads its table and shows the campaign view',
      async (view: any, campaignCode: string) => {
        const app = appOf(view);
        const c = await findCampaign(app, campaignCode);
        if (!c)
          return {success: false, error: `No campaign '${campaignCode}' — call listHitCampaigns for codes`};
        await app.infoView.setCampaign(c.name);
        if (app.campaignId !== c.name)
          return {success: false, error: 'Failed to open the campaign — you may not have permission to view it'};
        return {success: true, campaign: c.name, note: 'The campaign view is now open'};
      }),

    aiFunc('listHitTemplates', 'dynamic listHitTemplates(view view)',
      'List the campaign template names of this Hit app',
      async (view: any) => ({templates: await listTemplateNames(appOf(view).appName)})),

    aiFunc('getHitTemplateDetails', 'dynamic getHitTemplateDetails(view view, string templateName)',
      'Full detail of one campaign template - its campaign fields, stages, compute methods, data source and submit function',
      async (view: any, templateName: string) => {
        const app = appOf(view);
        const path = `${app.appName}/templates/${templateName}.json`;
        if (!await _package.files.exists(path))
          return {success: false, error: `No template '${templateName}' — call listHitTemplates for names`};
        const t = await loadTemplate<AnyTemplate>(path);
        const asDesign = t as HitDesignTemplate;
        const asTriage = t as HitTriageTemplate;
        return {
          name: t.name,
          key: (t as {key?: string}).key,
          campaignFields: (t.campaignFields ?? []).map((f) => ({name: f.name, type: f.type, required: !!f.required})),
          ...(asDesign.stages?.length ? {stages: asDesign.stages} : {}),
          ...(asTriage.dataSourceType ? {dataSourceType: asTriage.dataSourceType} : {}),
          ...(asTriage.queryFunctionName ? {dataSourceFunction: asTriage.queryFunctionName} : {}),
          compute: computeBrief(t),
          ...(t.submit ? {submitFunction: `${t.submit.package}:${t.submit.fName}`} : {}),
        };
      }),

    aiFunc('createHitCampaign', 'dynamic createHitCampaign(view view, string templateName, string campaignName, map campaignFields)',
      'Create a new campaign from a template and open it. Hit Design and PeptiHit only. campaignName and campaignFields are optional, but fields the template marks required must be provided',
      async (view: any, templateName: string, campaignName?: string, campaignFields?: {[key: string]: unknown}) => {
        const app = designAppOf(view);
        const path = `${app.appName}/templates/${templateName}.json`;
        if (!await _package.files.exists(path))
          return {success: false, error: `No template '${templateName}' — call listHitTemplates for names`};
        const template = await loadTemplate<HitDesignTemplate>(path);
        const missing = (template.campaignFields ?? [])
          .filter((f) => f.required && (campaignFields?.[f.name] == null || campaignFields[f.name] === ''))
          .map((f) => f.name);
        if (missing.length)
          return {success: false, error: `Missing required campaign fields: ${missing.join(', ')} — pass them in campaignFields`};
        app.clearCampaign();
        app.dataFrame = createNewHitDesignDataFrame(template, app.appName === 'PeptiHit');
        await app.setTemplate(template as any);
        app.campaignProps = {...(campaignFields ?? {})};
        const friendly = campaignName?.trim() ? campaignName.trim() : undefined;
        await app.saveCampaign(false, true, {friendlyName: friendly});
        if (template.layoutViewState && app.campaign)
          app.campaign.layout = template.layoutViewState;
        return {success: true, campaign: app.campaignId, note: 'The campaign is created and its design view is open'};
      }),

    aiFunc('deleteHitCampaign', 'dynamic deleteHitCampaign(view view, string campaignCode)',
      'Permanently delete a saved campaign. ALWAYS confirm with the user first',
      async (view: any, campaignCode: string) => {
        const app = appOf(view);
        const c = await findCampaign(app, campaignCode);
        if (!c)
          return {success: false, error: `No campaign '${campaignCode}' — call listHitCampaigns for codes`};
        if (!await canEditCampaign(c))
          return {success: false, error: 'You do not have permission to delete this campaign'};
        await app.infoView.deleteCampaignAndRefresh(app.appName, c.name);
        return {success: true, deleted: c.name};
      }),
  ];
}

// ─────────────────────── the campaign (table view) set ───────────────────────

function campaignSharedFunctions(): DG.Func[] {
  return [
    aiFunc('getHitCampaignState', 'dynamic getHitCampaignState(view view)',
      'State of the open campaign - name, template, status, molecule column, stages, columns, row counts, edit permission, submit function. Call this first. The molecule data itself lives in the attached table',
      async (view: any) => {
        const app = appOf(view);
        const df = app.dataFrame;
        const stages = app instanceof HitDesignApp ? app.stages : [];
        return {
          app: app.appName,
          campaign: app.campaign ? campaignBrief(app.campaign) : null,
          ...(app.campaignId ? {code: app.campaignId} : {}),
          moleculeColumn: app.molColName ?? null,
          ...((app as {helmColName?: string}).helmColName ? {sequenceColumn: (app as {helmColName?: string}).helmColName} : {}),
          ...(stages.length ? {stages} : {}),
          rows: df?.rowCount ?? 0,
          filteredRows: df?.filter.trueCount ?? 0,
          columns: df ? df.columns.toList().filter((c) => !c.name.startsWith('~'))
            .map((c) => ({name: c.name, type: c.type, ...(c.semType ? {semType: c.semType} : {})})) : [],
          canEdit: app.canEdit,
          submitFunction: app.submitParams ? `${app.submitParams.package}:${app.submitParams.fName}` : null,
          compute: computeBrief(app.campaign?.template ?? app.template),
        };
      }),

    aiFunc('saveHitCampaign', 'dynamic saveHitCampaign(view view)',
      'Save the open campaign - its table, layout and metadata',
      async (view: any) => {
        const app = appOf(view);
        if (!app.canEdit)
          return {success: false, error: 'You do not have permission to modify this campaign'};
        await app.saveCampaign();
        return {success: true, campaign: app.campaignId ?? app.campaign?.name};
      }),

    aiFunc('setHitCampaignStatus', 'dynamic setHitCampaignStatus(view view, string status)',
      'Set the status of the open campaign (free text, e.g. In Progress, Approved) and save it',
      async (view: any, status: string) => {
        const app = appOf(view);
        if (!status?.trim())
          return {success: false, error: 'status is required'};
        if (!app.campaign)
          return {success: false, error: 'No campaign is open'};
        if (!app.canEdit)
          return {success: false, error: 'You do not have permission to modify this campaign'};
        app.campaign.status = status.trim();
        await app.saveCampaign();
        return {success: true, status: app.campaign.status};
      }),

    aiFunc('submitHitCampaign', 'dynamic submitHitCampaign(view view, string status)',
      'Submit the campaign molecules to the configured submit function and save. status is optional and applies to Hit Design campaigns (Hit Triage submits always set the status to Submitted). ALWAYS confirm with the user first',
      async (view: any, status?: string) => {
        const app = appOf(view);
        if (!app.canEdit)
          return {success: false, error: 'You do not have permission to modify this campaign'};
        if (!app.submitParams)
          return {success: false, error: 'No submit function is configured for this campaign template'};
        const submitView = app.submitView;
        if (!submitView)
          return {success: false, error: 'The campaign is not fully open yet'};
        const ok = submitView.submitWithStatus ?
          await submitView.submitWithStatus(status ?? undefined) :
          await submitView.submit();
        if (!ok)
          return {success: false, error: 'The configured submit function could not be run'};
        return {success: true, submittedTo: `${app.submitParams.package}:${app.submitParams.fName}`};
      }),

    aiFunc('searchHitTableFunctions', 'dynamic searchHitTableFunctions(view view, string query, int limit)',
      'Search the standard table view functions (hundreds of data commands - add columns, aggregate, select rows, viewers). ALWAYS pass a query. Invoke matches with callHitTableFunction',
      async (view: any, query?: string, limit?: number) => {
        const tv = campaignTableView(view);
        const funcs = DG.View.prototype.getFunctions.call(tv) as DG.Func[];
        const terms = (query ?? '').toLowerCase().split(/\s+/).filter((t) => t.length > 0);
        const score = (f: DG.Func) => {
          if (!terms.length)
            return 1;
          const text = [f.name, f.friendlyName, f.description].filter(Boolean).join(' ').toLowerCase();
          return terms.filter((t) => text.includes(t)).length;
        };
        const matched = funcs.map((f) => ({f, s: score(f)})).filter((x) => x.s > 0)
          .sort((a, b) => b.s - a.s).map((x) => x.f);
        const max = limit && limit > 0 ? limit : DEFAULT_SEARCH_LIMIT;
        return {
          total: matched.length,
          functions: matched.slice(0, max).map((f) => ({
            name: f.name,
            description: f.description || f.friendlyName || f.name,
            inputs: f.inputs.filter((p) => p.name !== 'view' && (p.propertyType as string) !== 'view')
              .map((p) => ({name: p.name, type: p.propertyType})),
          })),
          ...(matched.length > max ? {note: `Showing ${max} of ${matched.length} — refine the query`} : {}),
        };
      }),

    aiFunc('callHitTableFunction', 'dynamic callHitTableFunction(view view, string name, map parameters)',
      'Invoke a standard table view function by name (from searchHitTableFunctions) with parameters keyed by input name',
      async (view: any, name: string, parameters?: {[key: string]: unknown}) => {
        const tv = campaignTableView(view);
        const funcs = DG.View.prototype.getFunctions.call(tv) as DG.Func[];
        const f = funcs.find((x) => x.name === name) ??
          funcs.find((x) => x.name.toLowerCase() === (name ?? '').toLowerCase());
        if (!f)
          return {success: false, error: `No table function '${name}' — call searchHitTableFunctions to find it`};
        const params: {[key: string]: unknown} = {...(parameters ?? {})};
        for (const inp of f.inputs) {
          if (!(inp.name in params) && (inp.name === 'view' || (inp.propertyType as string) === 'view'))
            params[inp.name] = tv;
        }
        const result = await f.apply(params);
        return {success: true, ...(result != null ? {result: serializeResult(result)} : {})};
      }),
  ];
}

function designCampaignFunctions(): DG.Func[] {
  return [
    ...campaignSharedFunctions(),

    aiFunc('setHitMolecule', 'dynamic setHitMolecule(view view, string molecule, int row)',
      'Set the molecule of a row (SMILES or molblock, HELM for PeptiHit). Omit row to add a new row. Registers the V-iD and runs the configured property calculations',
      async (view: any, molecule: string, row?: number) => {
        const app = designAppOf(view);
        if (!app.canEdit)
          return {success: false, error: 'You do not have permission to modify this campaign'};
        if (!molecule?.trim())
          return {success: false, error: 'molecule is required'};
        const df = app.dataFrame;
        if (!df)
          return {success: false, error: 'No campaign table is open'};
        const helmColName = (app as {helmColName?: string}).helmColName;
        const entryCol = helmColName ? df.col(helmColName) : df.col(app.molColName!);
        if (!entryCol)
          return {success: false, error: 'The molecule column was not found'};
        let idx = row ?? null;
        if (idx == null) {
          df.rows.addNew(null, true);
          idx = df.rowCount - 1;
        } else if (idx < 0 || idx >= df.rowCount)
          return {success: false, error: `row must be between 0 and ${df.rowCount - 1}`};
        entryCol.set(idx, molecule, true);
        if (helmColName) {
          // PeptiHit: the calculation converts the sequence and fills the molecule column
          await app.performSingleCellCalculations(idx, molecule);
          const mol = df.col(app.molColName!)?.get(idx);
          if (mol) {
            try {
              const vid = await registerMol(mol, app.campaignId!, app.appName);
              if (vid)
                df.col(ViDColName)?.set(idx, vid, true);
            } catch (e) {
              console.error('Failed to register molecule', e);
            }
          }
        } else {
          try {
            const vid = await registerMol(molecule, app.campaignId!, app.appName);
            if (vid)
              df.col(ViDColName)?.set(idx, vid, true);
          } catch (e) {
            console.error('Failed to register molecule', e);
          }
          await app.performSingleCellCalculations(idx, molecule);
        }
        return {success: true, row: idx, ...(df.col(ViDColName)?.get(idx) ? {vId: df.col(ViDColName)!.get(idx)} : {})};
      }),

    aiFunc('setHitStage', 'dynamic setHitStage(view view, int row, string stage)',
      'Move a row to another stage of the campaign progress tracker and save',
      async (view: any, row: number, stage: string) => {
        const app = designAppOf(view);
        if (!app.canEdit)
          return {success: false, error: 'You do not have permission to modify this campaign'};
        const df = app.dataFrame;
        const stageCol = df?.col(TileCategoriesColName);
        if (!df || !stageCol)
          return {success: false, error: 'This campaign has no stage column'};
        if (row == null || row < 0 || row >= df.rowCount)
          return {success: false, error: `row must be between 0 and ${df.rowCount - 1}`};
        if (!app.stages.includes(stage))
          return {success: false, error: `Unknown stage '${stage}' — the stages are: ${app.stages.join(', ')}`};
        stageCol.set(row, stage, true);
        await app.saveCampaign(false);
        return {success: true, row, stage};
      }),

    aiFunc('setHitStages', 'dynamic setHitStages(view view, list stages)',
      'Replace the stage list of the campaign progress tracker. Rows in a removed stage move to the first stage. Saves the campaign',
      async (view: any, stages: string[]) => {
        const app = designAppOf(view);
        if (!app.canEdit)
          return {success: false, error: 'You do not have permission to modify this campaign'};
        if (!Array.isArray(stages) || !stages.length || stages.some((s) => typeof s !== 'string' || !s.trim()))
          return {success: false, error: 'Pass stages as a non-empty list of names'};
        if (!app.campaign?.template || !app.template)
          return {success: false, error: 'No campaign is open'};
        if (!app.dataFrame?.col(TileCategoriesColName))
          return {success: false, error: 'This campaign has no stage column'};
        await app.setStages(stages.map((s) => s.trim()));
        return {success: true, stages: app.stages};
      }),

    aiFunc('recalculateHitRow', 'dynamic recalculateHitRow(view view, int row)',
      'Re-run the configured property calculations for one row of the campaign',
      async (view: any, row: number) => {
        const app = designAppOf(view);
        const df = app.dataFrame;
        if (!df)
          return {success: false, error: 'No campaign table is open'};
        if (row == null || row < 0 || row >= df.rowCount)
          return {success: false, error: `row must be between 0 and ${df.rowCount - 1}`};
        const helmColName = (app as {helmColName?: string}).helmColName;
        const value = (helmColName ? df.col(helmColName) : df.col(app.molColName!))?.get(row);
        if (!value)
          return {success: false, error: 'The row has no molecule to calculate from'};
        await app.performSingleCellCalculations(row, value);
        return {success: true, row};
      }),
  ];
}

// ─────────────────────── attachment ───────────────────────

export function hitInfoAiDescription(appName: AppName): string {
  const purpose = appName === 'Hit Design' || appName === 'PeptiHit' ?
    'design molecules in a spreadsheet, calculate properties, and track progress through stages' :
    'ingest a molecule dataset, calculate properties, and filter it down to the hits';
  return `${appName} — the landing page of a molecule campaign app where teams ${purpose}. ` +
    'Act on it through the view functions (search list_view_functions with "hit"): getHitAppInfo (call first), ' +
    'listHitCampaigns / getHitCampaignDetails / openHitCampaign to browse and continue campaigns, ' +
    'listHitTemplates / getHitTemplateDetails for the templates, createHitCampaign / deleteHitCampaign to manage them ' +
    '(confirm deletions with the user first).';
}

function campaignAiDescription(appName: AppName, kind: 'design' | 'pick'): string {
  const work = kind === 'design' ?
    'sketching and editing molecules in the grid. Hit-specific functions: getHitCampaignState (call first), ' +
    'setHitMolecule / setHitStage / setHitStages / recalculateHitRow to edit, ' +
    'saveHitCampaign / setHitCampaignStatus / submitHitCampaign to persist' :
    'filtering the enriched molecule table down to the hits. Hit-specific functions: getHitCampaignState (call first), ' +
    'saveHitCampaign / setHitCampaignStatus / submitHitCampaign';
  return `An open ${appName} campaign — a molecule spreadsheet for ${work}. ` +
    'The underlying table view also has hundreds of standard data commands — find them with ' +
    'searchHitTableFunctions (always pass a query) and run them with callHitTableFunction. ' +
    'The data itself lives in the attached table.';
}

/** Called from the app constructors: the MultiView (the shell view hosting the info
 * view) reports the info vocabulary. */
export function attachHitInfoAi(multiView: DG.MultiView, app: HitAppBase<AnyTemplate>): void {
  appByView.set(multiView, app);
  // TODO: Remove before api release
  (multiView as any).aiDescription = hitInfoAiDescription(app.appName);
  multiView.getFunctions = () => hitInfoViewFunctions();
}

/** Called when a campaign table view is created: only the hit-specific vocabulary is
 * reported (the standard table functions stay reachable through the search/call bridge). */
export function attachHitCampaignAi(view: DG.TableView, app: HitAppBase<AnyTemplate>, kind: 'design' | 'pick'): void {
  appByView.set(view, app);
  // TODO: Remove before api release
  (view as any).aiDescription = campaignAiDescription(app.appName, kind);
  view.getFunctions = () => kind === 'design' ? designCampaignFunctions() : campaignSharedFunctions();
}
