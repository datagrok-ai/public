/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Dialog, computed, signal} from '@datagrok-libraries/u2';
import {AppName, HitDesignCampaign, HitDesignTemplate, IFunctionArgs, TemplateCompute, TriagePermissions} from './types';
import {HitDesignInfoView} from './hit-design-views/info-view';
import {CampaignIdKey, CampaignJsonName, CampaignTableName, HTFunctionOrderingLSKey,
  HTQueryPrefix, HTScriptPrefix, HitDesignMolColName, TileCategoriesColName, ViDColName} from './consts';
import {calculateColumns, calculateCellValues} from './utils/calculate-single-cell';
// @ts-ignore
import '../../css/hit-triage.css';
// @ts-ignore
import './hit-design-views/utils.css';
import {_package} from '../package';
import {modifyUrl, toFormatedDateString} from './utils';
import {HitDesignSubmitView} from './hit-design-views/submit-view';
import {CampaignLockHeldError, HitAppBase} from './hit-app-base';
import {HitBaseView} from './base-view';
import {descriptorsTree, rerunOnOpenOption} from './dialogs/functions-dialog';
import {defaultPermissions} from './shared/permissions-editor';
import {computeEditor} from './shared/compute-editor';
import {CampaignSession, CampaignSessionOptions, openCampaignSession, pickMoleculeColumn}
  from './shared/campaign-session';
import {createNewHitDesignDataFrame} from './accordeons/new-hit-design-campaign-accordeon';
import {getDefaultCampaignStorageSettings, getDefaultSharingSettings} from '../packageSettingsEditor';
import {registerMol} from './utils/molreg';
import {PlaceholderTag, autoMergeOnCampaignOpen, isEmptyPlaceholder, mergeIntoCampaign, mergeServerTable}
  from './utils/merge-table';
import {attachHitInfoAi, attachHitCampaignAi} from './ai/view-functions';

interface SaveRequest {
  /** The table of the campaign the save was asked for. */
  df: DG.DataFrame;
  notify: boolean;
  isCreating: boolean;
  customProps?: Partial<HitDesignCampaign>;
}

/** A campaign besides its table, as a save writes it: the app's own while the campaign is open, else what the app
 * kept of it when another one opened. */
interface CampaignState {
  campaign?: HitDesignCampaign;
  template: HitDesignTemplate;
  campaignId: string;
  filePath: string;
  molColName: string;
  /** The table as the app last read it from the server or wrote it there: what a save over another session's save
   * tells the other session's edits by. */
  base?: DG.DataFrame;
  canEdit: boolean;
  /** A new campaign's fields. */
  fields: {[key: string]: any};
  session?: CampaignSession;
}

export class HitDesignApp<T extends HitDesignTemplate = HitDesignTemplate> extends HitAppBase<T> {
  multiView: DG.MultiView;
  protected _infoView: HitDesignInfoView;
  get infoView(): HitDesignInfoView {return this._infoView;}
  public _submitView?: HitDesignSubmitView;
  get submitView(): HitDesignSubmitView | undefined {return this._submitView;}
  protected _designViewName = 'Design';
  protected _filePath = `System.AppData/HitTriage/${this.appName}/campaigns`;
  protected _campaignId?: string;
  protected _molColName: string = HitDesignMolColName;
  protected _campaign?: HitDesignCampaign;
  public campaignProps: {[key: string]: any} = {};
  public mainView: DG.ViewBase;
  protected get version() {return this._campaign?.version ?? 0;};
  public existingStatuses: string[] = [];
  /** Merges running: the rows they add are not the user's. */
  public joining = 0;
  private _session?: CampaignSession;
  private _base?: DG.DataFrame;
  /** Saves asked for and not finished. */
  private _saves = 0;
  /** The last save asked for: the next one starts after it. */
  private _queue: Promise<unknown> = Promise.resolve();
  /** The save waiting for the running one: a save of its campaign asked for meanwhile is folded into it. */
  private _nextSave?: {request: SaveRequest, done: Promise<HitDesignCampaign | undefined>};
  /** What the app kept of each campaign it switched away from, for a save of it still to run. */
  private readonly _closed = new WeakMap<DG.DataFrame, CampaignState>();
  private _opens = 0;
  /** The latest single-cell calculation per row: an earlier one still running drops its results. */
  private readonly _rowRuns = new Map<number, number>();

  constructor(c: DG.FuncCall, an: AppName = 'Hit Design',
    infoViewConstructor: (app: HitDesignApp) => HitDesignInfoView = (app) => new HitDesignInfoView(app)) {
    super(c, an);
    this._infoView = infoViewConstructor(this);
    this.multiView = new DG.MultiView({viewFactories: {[this._infoView.name]: () => this._infoView}});
    this.multiView.helpUrl = this._infoView.helpUrl;
    this.multiView.subs.push(this.multiView.tabs.onTabChanged.subscribe((_) => {
      if (this.multiView.currentView instanceof HitBaseView)
        (this.multiView.currentView as HitBaseView<T, typeof this>).onActivated();
    }));
    this.multiView.subs.push(grok.events.onCurrentViewChanged.subscribe(() => {
      if (grok.shell.v === this.mainView)
        this.setBaseUrl();
    }));

    this.multiView.parentCall = c;
    attachHitInfoAi(this.multiView, this);

    this.mainView = this.multiView;
  }

  /** The open campaign's view. */
  get designView(): DG.TableView | undefined {return this._session?.view;}

  /** A save is running or waiting for the one that runs. */
  get isSaving(): boolean {return this._saves > 0;}

  /** How the campaign view of this app works: Hit Design's. */
  protected sessionOptions(): CampaignSessionOptions {
    const layout = this.campaign?.layout ?? (this.campaign ? undefined : this.template?.layoutViewState);
    return {entryColumn: this.molColName, editor: 'Chem:editMoleculeCell', duplicateLabel: 'Duplicate molecule',
      registerFirst: true, designer: true, layout};
  }

  /** Gets the campaign's table ready before its view shows it. */
  protected prepareTable(_df: DG.DataFrame): void {}

  /**
   * Correctly joins the incoming dataframe with the existing campaign dataframe.
   * Drag-n-drop entry point — uses SMILES-based join with the historical defaults.
   * @param df Dataframe to be joined into existing campaign dataframe
   */
  async handleJoiningDataframe(df: DG.DataFrame) {
    const molCols = HitAppBase.moleculeColumns(df);
    if (molCols.length === 0) {
      grok.shell.error('No molecule column found');
      return;
    }

    let chosenMolColName: string | null = molCols[0]!.name;
    if (molCols.length > 1) {
      const likely = molCols.find((c) => c?.name?.toLowerCase() === 'smiles' || c?.name?.toLowerCase() === 'molecule') ?? molCols[0];
      const session = this._session;
      const table = this.dataFrame;
      chosenMolColName = await pickMoleculeColumn(molCols.map((c) => c.name), likely.name, session?.host);
      // the table goes into the campaign it was dropped on, or nowhere
      if (session?.host.scope.isDisposed || this.dataFrame !== table)
        return;
      if (!chosenMolColName) {
        grok.shell.info('Operation cancelled by user');
        return;
      }
    }

    if (!df.col(chosenMolColName)) {
      grok.shell.error('No molecule column found');
      return;
    }

    await mergeIntoCampaign(this, df, {
      mode: 'smiles',
      molColName: chosenMolColName,
      addNewRows: true,
      runComputeOnNewRows: true,
      clashStrategy: 'overwrite',
    });
    await this.saveCampaign(false);
  }

  public get submitParams() {
    return this.campaign?.template?.submit ?? this.template?.submit;
  }

  public get stages() {
    return this.campaign?.template?.stages ?? this.template?.stages ?? [];
  }

  /** The campaign's stages from now on, duplicates merged; rows in a stage that goes move to the first one. Saved;
   * answers whether the save happened (a balloon says why not). */
  public async setStages(st: string[]): Promise<boolean> {
    if (!this.campaign || !this.campaign.template || !this.template) {
      grok.shell.error('Campaign or template is not set');
      return false;
    }
    if (!st?.length) {
      grok.shell.error('Removing all stages is not allowed');
      return false;
    }
    const stageCol = this.dataFrame?.col(TileCategoriesColName);
    if (!stageCol) {
      grok.shell.error('No stage column found');
      return false;
    }
    const stageSet = new Set(st);
    const uniqueStages = [...stageSet];
    for (let i = 0; i < this.dataFrame!.rowCount; i++) {
      if (!stageSet.has(stageCol.get(i)))
        stageCol.set(i, st[0], false);
    }
    this.campaign.template.stages = uniqueStages;
    this.template.stages = uniqueStages;
    stageCol.setTag(DG.TAGS.CHOICES, JSON.stringify(uniqueStages));
    return await this.saveCampaign(true) !== undefined;
  }

  /** Reads the campaign's table from the server again into a new view; the rest of the campaign stays as it is. */
  public async refreshCampaign() {
    if (!this._campaignId || !this.template) {
      grok.shell.warning('No campaign to refresh');
      return; // should not happen
    }
    if (this.isSaving) {
      grok.shell.info('Campaign is being saved. Please wait a moment and try again.');
      return;
    }
    await this.setTemplate(this.template, this._campaignId);
  }

  /** The ticket of an open asked for now: of the opens asked for, only the latest one takes its campaign. */
  openTicket(): number {
    return ++this._opens;
  }

  /** Opens the campaign of `campaignId` (`campaign` is its campaign.json, by default the open one: a refresh) with
   * `template`, or a new campaign of the template without one. Everything is read first; then the app takes the
   * campaign and its view at once, so no save sees half of one campaign and half of another. An open asked for later
   * wins: `ticket` ({@link openTicket}) is taken when the open is asked for, before what the caller reads first.
   * Answers whether this one took the campaign. */
  public async setTemplate(template: T, campaignId?: string,
    campaign: HitDesignCampaign | undefined = campaignId ? this._campaign : undefined,
    ticket = this.openTicket()): Promise<boolean> {
    let df: DG.DataFrame;
    let filePath: string;
    if (!campaignId) {
      campaignId = await this.getNewCampaignName(`${this.appName}/campaigns`, template.key);
      let defaultFolderPath = await getDefaultCampaignStorageSettings();
      if (!defaultFolderPath.endsWith('/'))
        defaultFolderPath += '/';
      filePath = `${defaultFolderPath}${this.appName}/campaigns/${campaignId}/${CampaignTableName}`;
      df = this.newCampaignFrame(template);
    } else {
      const fileLoc = `System.AppData/HitTriage/${this.appName}/campaigns`;
      filePath = campaign?.savePath ?? `${fileLoc}/${campaignId}/${CampaignTableName}`;
      df = await grok.dapi.files.readCsv(filePath);
      for (const name of df.columns.names().filter((name) => name.startsWith('~')))
        df.columns.remove(name);
      for (const [colName, semtype] of Object.entries(campaign?.columnSemTypes ?? {})) {
        const col = df.columns.byName(colName);
        if (col && semtype)
          col.semType = semtype;
      }
      for (const [colName, type] of Object.entries(campaign?.columnTypes ?? {})) {
        const col = df.columns.byName(colName);
        if (!col || col.type === type)
          continue;
        try {
          df.changeColumnType(colName, type as DG.COLUMN_TYPE);
        }
        catch (e) {
          console.error(e);
        }
      }
      // restored before the view binds the table, which subscribes formula columns to recalculation only if it has any
      for (const [colName, tags] of Object.entries(campaign?.columnTags ?? {})) {
        const col = df.columns.byName(colName);
        if (!col)
          continue;
        for (const [tag, value] of Object.entries(tags))
          !!tag && !!value && col.setTag(tag, value);
        // the tag alone does not reorder the categories
        if (tags['.category-order'] && col.type === DG.COLUMN_TYPE.STRING)
          col.setCategoryOrder(col.getCategoryOrder());
      }
      // the CSV reader formats a number column to every digit it read ('0.0000000000000' for a calculated MW): unless
      // the campaign chose another format, the column shows as the grid shows numbers, as it did before the reload
      for (const col of df.columns.toList()) {
        const digits = col.getTag(DG.TAGS.SOURCE_PRECISION);
        if (digits !== null && col.getTag(DG.TAGS.FORMAT) === `0.${'0'.repeat(parseInt(digits, 10))}`) {
          col.tags.delete(DG.TAGS.FORMAT);
          col.tags.delete(DG.TAGS.SOURCE_PRECISION);
        }
      }
    }

    if (!df) {
      grok.shell.error('DataFrame is empty');
      return false;
    }
    const base = df.clone();
    await df.meta.detectSemanticTypes();
    const molColName = HitAppBase.moleculeColumn(df)?.name ?? HitDesignMolColName;
    const canEdit = await this.canEditCampaign(campaign);

    // Implicit merge pass before the view is shown: only runs when the campaign carries a
    // mergeConfig with autoMergeOnOpen and a remembered fileshare path. Compute is forced off
    // and we deliberately do not save — the user can save manually after the merge.
    if (campaign) {
      await autoMergeOnCampaignOpen(this, campaign.mergeConfig, {df, molColName, campaignId,
        stages: campaign.template?.stages ?? template.stages ?? []});
    }

    // Recompute any compute functions flagged "re-run when campaign opens" over the whole molecule
    // column. Like the auto-merge pass, this runs only for existing campaigns and is not auto-saved.
    if (campaign)
      await this.reComputeOnCampaignOpen(df, template, molColName);
    if (ticket !== this._opens)
      return false;

    if (this.dataFrame && this.dataFrame !== df)
      this._closed.set(this.dataFrame, {...this.state(), session: undefined});
    this._campaign = campaign;
    if (campaign)
      campaign.template = template;
    this.template = template;
    this._campaignId = campaignId;
    this._filePath = filePath;
    this.dataFrame = df;
    this._base = base;
    this._molColName = molColName;
    this.hasEditPermission = canEdit;
    this._submitView ??= new HitDesignSubmitView(this);
    grok.shell.windows.showHelp = false;
    this.openView();
    this.setBaseUrl();
    modifyUrl(CampaignIdKey, campaignId);
    return true;
  }

  /** The campaign's view, the one open before closed first. */
  private openView(): void {
    this._session?.view.close();
    const df = this.dataFrame!;
    this.prepareTable(df);
    this._designViewName = this.campaign?.friendlyName ?? this.campaign?.name ?? this._designViewName;
    const view = grok.shell.addTableView(df);
    view.name = this._designViewName;
    view.parentCall = this.parentCall;
    attachHitCampaignAi(view, this, 'design');
    const session = this._session = openCampaignSession(this, view, this.sessionOptions());
    session.host.own(() => {
      if (this._session === session)
        this._session = undefined;
    });
  }

  get campaignId(): string | undefined {return this._campaignId;}

  get designViewName(): string {return this._designViewName;}

  set designViewName(name: string) {
    this._designViewName = name;
    if (this.designView)
      this.designView.name = name;
  }

  get molColName() {
    return this._molColName ??= (this.dataFrame && HitAppBase.moleculeColumn(this.dataFrame)?.name) ?? HitDesignMolColName;
  }

  /** Where the campaign's table is kept. */
  get filePath(): string {return this._filePath;}
  set filePath(path: string) {this._filePath = path;}

  /** The folder the campaign's own folder is in. */
  get campaignFolder(): string {
    const at = this._filePath.indexOf((this.campaignId ?? this._campaign?.name)!);
    return at === -1 ? this._filePath : this._filePath.substring(0, at - 1);
  }

  get campaign(): HitDesignCampaign | undefined {return this._campaign;}

  set campaign(campaign: HitDesignCampaign | undefined) {this._campaign = campaign;}

  /** A new campaign from `template`, opened and saved at once: nothing of the campaign open before carries over, the
   * name typed names it, and its first save takes the code that is free then. */
  async startNewCampaign(template: T, fields: {[key: string]: unknown}, name?: string): Promise<void> {
    if (!await this.setTemplate(template))
      return;
    this.campaignProps = fields;
    await this.saveCampaign(false, true, {friendlyName: name});
  }

  /** The one-row table a new campaign starts from. */
  protected newCampaignFrame(template: T): DG.DataFrame {
    return createNewHitDesignDataFrame(template);
  }

  private _duplicateVidCache?: {
    colVersion: number,
    valueCounts: Int32Array | Uint32Array,
    indexes: Int32Array | Uint32Array,
  };

  private cacheDuplicateVIDs() {
    const col = this.dataFrame?.col(ViDColName);
    if (!col) {
      this._duplicateVidCache = undefined;
      return;
    }
    const colVersion = col.version;
    const valueCounts = new Uint32Array(col.categories.length).fill(0);
    const indexes = (col.getRawData() as Uint32Array).subarray(0, col.length); // because for small data DG returns list of length min 16
    indexes.forEach((v) => valueCounts[v]++);
    this._duplicateVidCache = {colVersion, valueCounts, indexes};
  }

  get duplicateVidCache() {
    const col = this.dataFrame?.col(ViDColName);
    if (!col)
      return undefined;
    if (!this._duplicateVidCache || col.version !== this._duplicateVidCache.colVersion)
      this.cacheDuplicateVIDs();
    return this._duplicateVidCache;
  }

  /** Starts the row's calculation; answers its number, which {@link isLatestRun} checks after each await. */
  protected nextRun(row: number): number {
    const run = (this._rowRuns.get(row) ?? 0) + 1;
    this._rowRuns.set(row, run);
    return run;
  }

  protected isLatestRun(row: number, run: number): boolean {
    return this._rowRuns.get(row) === run;
  }

  /** The template's compute on the row's molecule, the results into the row (new columns added); saved. A later
   * calculation of the row started meanwhile wins. */
  async performSingleCellCalculations(newValueIdx: number, newValue?: string) {
    const computeObj = this.template!.compute;
    // OK on an empty sketcher writes an empty molfile into a molblock cell: no molecule to calculate
    if (!newValue || !_package.convertToSmiles(newValue))
      return;
    const df = this.dataFrame!;
    const run = this.nextRun(newValueIdx);
    const calcDf =
              await calculateCellValues(
                [newValue], computeObj.descriptors.args, computeObj.functions, computeObj.scripts, computeObj.queries);
    if (!this.isLatestRun(newValueIdx, run) || this.dataFrame !== df)
      return;
    for (const col of calcDf.columns.toList()) {
      if (col.name === HitDesignMolColName) continue;
      const target = df.col(col.name);
      if (!target) {
        const newCol = df.columns.addNew(col.name, col.type);
        newCol.semType = col.semType;
      }
      else if (target.type !== col.type && isEmptyPlaceholder(target) &&
        // a number placeholder keeps counts as they come: a retype would take the layout's filter on it away
        !(target.type === DG.COLUMN_TYPE.FLOAT && col.type === DG.COLUMN_TYPE.INT))
        df.changeColumnType(col.name, col.type as DG.COLUMN_TYPE);
      df.col(col.name)!.set(newValueIdx, col.get(0), false);
    }
    df.fireValuesChanged();
    void this.saveCampaign(false);
  }

  /** The V-iD of the row's molecule, registered with the campaign; written only while the row still holds that
   * molecule, so the registration of a later edit wins. */
  async registerRow(row: number): Promise<void> {
    const df = this.dataFrame;
    const molecules = df?.col(this.molColName);
    const molecule = molecules?.get(row);
    if (!df || !molecule)
      return;
    const vid = await registerMol(molecule, this._campaignId!, this.appName);
    if (vid && molecules!.get(row) === molecule)
      df.col(ViDColName)?.set(row, vid, true);
  }

  /**
   * Recomputes the compute functions/scripts/queries flagged with `rerunOnOpen` over the whole
   * molecule column of `df`, the table of a campaign being opened, before the design view is
   * shown. Mirrors the auto-merge-on-open pass: values are updated in memory but the campaign is
   * not auto-saved, and failures never bubble up so the campaign still opens.
   */
  private async reComputeOnCampaignOpen(df: DG.DataFrame, template: T, molColName: string): Promise<void> {
    const compute = template.compute;
    if (!compute)
      return;

    // Descriptors are deterministic from structure, so they are intentionally never re-run on open.
    const functions = (compute.functions ?? []).filter((f) => f.rerunOnOpen);
    const scripts = (compute.scripts ?? []).filter((s) => s.rerunOnOpen);
    const queries = (compute.queries ?? []).filter((q) => q.rerunOnOpen);
    if (!functions.length && !scripts.length && !queries.length)
      return;

    const molCol = df.col(molColName);
    if (!molCol)
      return;

    try {
      // Compute against a temporary molecule-only table, then merge the results back in. Merging by
      // column name lets us overwrite an existing column in place instead of appending a duplicate
      // "<name> (2)" column — e.g. a query that re-resolves a ChEMBL molregno already in the table.
      const values = molCol.toList();
      const calcDf = await calculateCellValues(values, [], functions, scripts, queries);
      for (const resCol of calcDf.columns.toList()) {
        if (resCol.name === molColName || resCol.name === HitDesignMolColName)
          continue;
        let targetCol = df.col(resCol.name);
        if (!targetCol || targetCol.type !== resCol.type) {
          if (targetCol) // existing column of a different type: replace it wholesale
            df.columns.remove(resCol.name);
          targetCol = df.columns.addNew(resCol.name, resCol.type);
        }
        targetCol.semType = resCol.semType;
        const rows = Math.min(resCol.length, targetCol.length);
        for (let i = 0; i < rows; i++)
          targetCol.set(i, resCol.isNone(i) ? null : resCol.get(i), false);
      }
      df.fireValuesChanged();
    } catch (e) {
      _package.logger.error(e);
      grok.shell.warning('Re-run on open failed for some functions; opening the campaign anyway.');
    }
  }

  /** "Calculate additional properties": the campaign's compute in the compute editor; on OK what is new or changed
   * is calculated over the whole table, and the campaign saved. */
  editCompute(): void {
    const dialog = new Dialog('Compute', {name: 'compute'});
    dialog.build(() => {
      const editor = computeEditor({functions: this.computeFunctions,
        value: this.campaign?.template?.compute ?? this.template!.compute, descriptors: descriptorsTree,
        orderKey: HTFunctionOrderingLSKey, rerunOnOpen: rerunOnOpenOption});
      // the descriptors picked are written in the tree's order, as the compute dialog wrote them
      const listed = signal(false);
      void editor.ready.then(() => listed.value = true);
      dialog.add(editor);
      dialog.okDisabled(computed(() => listed.value ? null : 'Listing the descriptors…'), {footer: true});
      dialog.onOK(() => void this.applyCompute(editor.read()));
    });
    this._session!.host.show(dialog, {width: 800, resizable: true});
  }

  /** The campaign's compute from now on: what it calculates anew or with other settings is calculated over the
   * whole table, and the campaign saved. */
  private async applyCompute(compute: TemplateCompute): Promise<void> {
    const old = this.template!.compute;
    const df = this.dataFrame!;
    this.template!.compute = compute;
    if (this.campaign)
      this.campaign.template = this.template;
    const changed = <F extends {args: IFunctionArgs}>(items: F[], before: F[], key: (f: F) => string) =>
      Object.fromEntries(items.filter((f) => {
        const was = before.find((b) => key(b) === key(f));
        return !was || !Object.entries(f.args).every(([name, value]) => was.args[name] === value);
      }).map((f) => [key(f), f.args]));
    const grid = this.designView?.grid;
    grid && ui.setUpdateIndicator(grid.root, true);
    try {
      await calculateColumns({
        descriptors: compute.descriptors.args.filter((d) => !old.descriptors.args.includes(d)),
        externals: changed(compute.functions, old.functions, (f) => `${f.package}:${f.name}`),
        scripts: changed(compute.scripts ?? [], old.scripts ?? [], (s) => `${HTScriptPrefix}:${s.name}:${s.id}`),
        queries: changed(compute.queries ?? [], old.queries ?? [], (q) => `${HTQueryPrefix}:${q.name}:${q.id}`),
      }, df, this.molColName!);
    } finally {
      grid && ui.setUpdateIndicator(grid.root, false);
      if (this.dataFrame === df)
        void this.saveCampaign(false);
    }
  }

  /** Saves the open campaign: one save at a time; a save of a campaign asked for while one runs is folded into the
   * next one of that campaign. Answers the saved campaign, or nothing when it was not saved (a balloon says why). */
  saveCampaign(notify = true, isCreating = false, customProps?: Partial<HitDesignCampaign>): Promise<HitDesignCampaign | undefined> {
    const df = this.dataFrame;
    if (!df)
      return Promise.resolve(undefined);
    const next = this._nextSave;
    if (next?.request.df === df) {
      next.request.notify ||= notify;
      next.request.isCreating ||= isCreating;
      next.request.customProps ??= customProps;
      return next.done;
    }
    const request: SaveRequest = {df, notify, isCreating, customProps};
    const queued = this._saves > 0;
    this._saves++;
    const start = () => {
      if (this._nextSave?.request === request)
        this._nextSave = undefined;
      return this._write(request);
    };
    const done = (queued ? this._queue.then(start, start) : start()).finally(() => this._saves--);
    this._queue = done;
    if (queued)
      this._nextSave = {request, done};
    return done;
  }

  /** The open campaign besides its table. */
  private state(): CampaignState {
    return {campaign: this._campaign, template: this.template!, campaignId: this._campaignId!, filePath: this._filePath,
      molColName: this.molColName, base: this._base, canEdit: this.hasEditPermission, fields: this.campaignProps,
      session: this._session};
  }

  /** What a save changed of the campaign of `df`: the app's own if it still has it open, else what it kept of it. */
  private recordSave(df: DG.DataFrame, changes: Partial<CampaignState>): void {
    if (this.dataFrame !== df) {
      const kept = this._closed.get(df);
      if (kept)
        Object.assign(kept, changes);
      return;
    }
    if (changes.campaign)
      this._campaign = changes.campaign;
    if (changes.base)
      this._base = changes.base;
    if (changes.filePath !== undefined)
      this._filePath = changes.filePath;
    if (changes.campaignId !== undefined) {
      this._campaignId = changes.campaignId;
      modifyUrl(CampaignIdKey, changes.campaignId);
    }
  }

  /** One save of the campaign of `df` under the campaign's lock: refused before the lock when the user may not edit
   * the campaign. The app shows the saved campaign only if it still has that one open. */
  private async _write({df, notify, isCreating, customProps}: SaveRequest): Promise<HitDesignCampaign | undefined> {
    const open = () => this.dataFrame === df;
    const state = open() ? this.state() : this._closed.get(df);
    if (!state)
      return undefined;
    if (!state.canEdit) {
      grok.shell.error('You do not have permission to modify this campaign');
      return undefined;
    }
    // a closed view's session shows no feedback: its timers would start on a disposed host
    const session = state.session?.host.scope.isDisposed ? undefined : state.session;
    let {campaignId, filePath} = state;
    session?.waitingForLock(true);
    try {
      if (state.campaign) {
        const lock = await this.aquireAndWaitCampaignLock(campaignId);
        if (!lock?.aquired)
          throw new CampaignLockHeldError('Could not acquire lock for campaign during 30 seconds');
      } else {
        const locked = await this.lockNewCampaignCode(campaignId, state.template.key);
        if (locked !== campaignId) {
          filePath = filePath.replace(`/${campaignId}/${CampaignTableName}`, `/${locked}/${CampaignTableName}`);
          campaignId = locked;
          this.recordSave(df, {campaignId, filePath});
        }
      }
    }
    catch (e) {
      _package.logger.error(e);
      // only a lock another session kept is "someone else"; a query that failed says why
      grok.shell.error(e instanceof CampaignLockHeldError ?
        'Failed to acquire campaign lock. Someone else might be editing the campaign. Please try again later.' :
        `Failed to save the campaign: ${e instanceof Error ? e.message : e}`);
      return undefined;
    } finally {
      session?.waitingForLock(false);
    }
    session?.saving(true);
    const save: CampaignState = {...state, campaignId, filePath};
    const version = state.campaign?.version ?? 0;
    let written = false;
    try {
      const serverVersion = await this.mergeStaleSave(df, save);
      const campaign = await this.campaignToSave(df, save, customProps);
      if (serverVersion !== undefined)
        campaign.version = Math.max(version, serverVersion) + 1;
      const saved = df.columns.toList().filter((col) => !col.name.startsWith('~') && !isEmptyPlaceholder(col));
      const csvDf = DG.DataFrame.fromColumns(saved).toCsv();
      const base = df.clone();
      await grok.dapi.files.writeAsText(filePath, csvDf);
      written = true;
      // a view closed meanwhile keeps the layout saved before
      const view = session && !session.host.scope.isDisposed ? session.view : undefined;
      const newLayout = view?.saveLayout();
      if (view && !newLayout)
        grok.shell.warning('Layout cound not be saved');
      else if (newLayout)
        campaign.layout = newLayout.viewState;
      campaign.template = save.template;

      await _package.saveCampaignJson(this.appName, campaign);
      notify && grok.shell.info('Campaign saved successfully.');
      !notify && isCreating && grok.shell.info('Campaign created successfully.');
      const updatedAt = await this.releaseCampaignLock(campaignId);
      this.recordSave(df, {campaign, base});
      if (open()) {
        this.designViewName = campaign.friendlyName ?? campaign.name;
        session?.saved(updatedAt);
      }
      return campaign;
    } catch (e) {
      _package.logger.error(e);
      // nothing written yet: the lock goes without logging a save other sessions would be told of
      await (written ? this.releaseCampaignLock(campaignId) : this.dropCampaignLock(campaignId))
        .catch((err) => _package.logger.error(err));
      grok.shell.error(`Failed to save the campaign: ${e instanceof Error ? e.message : e}`);
      return undefined;
    } finally {
      session?.saving(false);
    }
  }

  /** campaign.json as the save writes it: the keys this app knows, and every key of the stored one it does not. A
   * layout's column nothing filled yet is left out. */
  private async campaignToSave(df: DG.DataFrame, save: CampaignState,
    customProps?: Partial<HitDesignCampaign>): Promise<HitDesignCampaign> {
    const columns = df.columns.toList().filter((col) => !col.name.startsWith('~') && !isEmptyPlaceholder(col));
    const columnSemTypes: {[_: string]: string} = {};
    const colTypeMap: {[_: string]: string} = {};
    const columnTags: {[_: string]: {[_: string]: string}} = {};
    for (const col of columns) {
      columnSemTypes[col.name] = col.semType;
      colTypeMap[col.name] = col.type;
      // the Stage choices are the template's stages, set on every open: a stored copy would go stale
      columnTags[col.name] = Object.fromEntries(Array.from(col.tags).filter(([tag]) => tag !== PlaceholderTag &&
        !(col.name === TileCategoriesColName && tag === DG.TAGS.CHOICES)));
    }

    const getDefaultPerms = async () => {
      try {
        const perms = await getDefaultSharingSettings();
        const idPerms: TriagePermissions = {view: perms.view.map((p) => p.id), edit: perms.edit.map((p) => p.id)};
        if (idPerms.edit.length === 0)
          idPerms.edit = [defaultPermissions.edit[0]];
        if (idPerms.view.length === 0)
          idPerms.view = [defaultPermissions.view[0]];
        return idPerms;
      } catch (e) {
        grok.shell.error('Failed to get default permissions');
        console.error(e);
      }
      return defaultPermissions;
    };
    const previous = save.campaign;
    // if its first time save author as current user, else keep the same
    const authorUserId = previous?.authorUserId ?? grok.shell.user.id;
    const permissions = previous?.permissions ?? await getDefaultPerms();
    const authorName = authorUserId ? previous?.authorUserFriendlyName ?? (await grok.dapi.users.find(authorUserId))?.friendlyName : undefined;
    return {
      ...previous,
      name: save.campaignId,
      friendlyName: previous?.friendlyName ?? customProps?.friendlyName ?? save.campaignId,
      templateName: save.template.name,
      status: previous?.status ?? 'In Progress',
      createDate: previous?.createDate ?? toFormatedDateString(new Date()),
      campaignFields: previous?.campaignFields ?? save.fields,
      columnSemTypes,
      rowCount: df.rowCount,
      filteredRowCount: df.filter.trueCount,
      savePath: save.filePath,
      columnTypes: colTypeMap,
      columnTags,
      version: (previous?.version ?? 0) + 1,
      authorUserId,
      authorUserFriendlyName: authorName,
      lastModifiedUserName: grok.shell.user.friendlyName,
      permissions,
    };
  }

  /** Another session saved the campaign since this one read it: its table, where its campaign.json says, is merged
   * into `df`. Answers the version it saved, nothing when there was no such save. */
  private async mergeStaleSave(df: DG.DataFrame, save: CampaignState): Promise<number | undefined> {
    const campaignPath = `${this.appName}/campaigns/${save.campaignId}/${CampaignJsonName}`;
    if (!await _package.files.exists(campaignPath))
      return undefined;
    const server: HitDesignCampaign = JSON.parse(await _package.files.readAsText(campaignPath));
    const serverVersion = server.version ?? 0;
    if (serverVersion <= (save.campaign?.version ?? 0))
      return undefined;
    const table = server.savePath ?? save.filePath;
    if (await grok.dapi.files.exists(table))
      mergeServerTable(this, df, save.molColName, await grok.dapi.files.readCsv(table), save.base);
    return serverVersion;
  }
}
