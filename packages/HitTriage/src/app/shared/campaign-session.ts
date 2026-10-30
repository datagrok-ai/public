/* A campaign open in Hit Design or PeptiHit: the platform TableView over its table, and everything the app adds to
   it — the ribbon group, the crumb back to the landing, the address, what the grid's edits start, its menu items, a
   dropped file, the check for other sessions' saves, the save feedback and the progress tracker. All of it is owned
   by the view's host and ends when the view closes. The two apps differ only by {@link CampaignSessionOptions}. */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ChoiceInput, Dialog, computed, signal, span, subscribe} from '@datagrok-libraries/u2';
import type {ReadonlySignal} from '@datagrok-libraries/u2';
import {ViewHost, funcForm, isView} from '@datagrok-libraries/u2/dg';
import type {RibbonAction, RibbonButton} from '@datagrok-libraries/u2/dg';
import type {HitDesignApp} from '../hit-design-app';
import {HitAppBase} from '../hit-app-base';
import {CampaignIdKey, HitDesignerFunctionTag, TileCategoriesColName, ViDColName, ViDSemType} from '../consts';
import {modifyUrl} from '../utils';
import {_package} from '../../package';
import * as api from '../../package-api';
import {openMergeTableDialog} from '../dialogs/merge-table-dialog';
import {MERGE_ALLOWED_EXTENSIONS, PlaceholderTag, readFileToDataFrame} from '../utils/merge-table';
import {showProgressTracker} from '../hit-design-views/tiles-view';
import {PermissionsEditor} from './permissions-editor';

export interface CampaignSessionOptions {
  /** The column a design is entered in: Molecule (Hit Design) or Helm (PeptiHit). */
  entryColumn: string;
  /** Opens the editor of an entry cell: Chem's sketcher or the HELM editor. */
  editor: string;
  /** The grid menu's item that copies an entry into a new row. */
  duplicateLabel: string;
  /** An entry's V-iD is registered before its calculation (Hit Design), or after the calculation converted the
   * entry into the molecule the V-iD is of (PeptiHit). */
  registerFirst: boolean;
  /** The designer function's button (Hit Design). */
  designer: boolean;
  /** The layout the view opens with: the campaign's, else the template's for a new campaign (both apps) or for any
   * campaign without one (PeptiHit). The columns its viewers are bound to that the table lacks (a new campaign's
   * computed ones) are added empty, so those viewers start empty instead of failing; a save leaves such a column out
   * until it holds values. */
  layout?: string;
}

/** The session of a campaign view, as the app's saves see it. */
export interface CampaignSession {
  readonly view: DG.TableView;
  /** The view's host, which owns the session: it ends when the view closes. */
  readonly host: ViewHost;
  /** Another session saved the campaign since this one read it: the Refresh icon says so. */
  readonly needsRefresh: ReadonlySignal<boolean>;
  /** The save's feedback: the grid covered while it writes. */
  saving(on: boolean): void;
  /** The save waits for another session's lock: said in the status bar once the wait is noticeable. */
  waitingForLock(on: boolean): void;
  /** A save of this session ended at `updatedAt`: that is no other session's change. */
  saved(updatedAt: string | null): void;
  /** The progress tracker over the view, unless it is open already; answers its dialog. */
  openTracker(): Dialog;
}

const NEW_CHANGES = 'New changes available. Click to refresh campaign.';
const MULTIPLE_MOLECULES = 'Multiple molecule columns found. ' +
  'Please select the primary molecule column to be used for calculations.';

/** The session of `view`, built under the view's {@link ViewHost}: it ends when the view closes. */
export function openCampaignSession(app: HitDesignApp, view: DG.TableView,
  options: CampaignSessionOptions): CampaignSession {
  const host = ViewHost.of(view);
  return host.build(() => {
    const needsRefresh = signal(false);
    const waiting = signal(false);
    let cancelWaitNote: (() => void) | undefined;
    let lastModified: string | null = null;
    let checking = false;
    let tracker: Dialog | undefined;

    // before the signals: a layout without a grid is reset to one, which tells no one
    if (options.layout !== undefined)
      loadLayout(options.layout);
    else
      checkGrid();
    const signals = host.signals;
    prepareColumns();
    follow();
    chrome();
    openEditorOnNew();
    watchForChanges();

    return {
      view, host, needsRefresh,
      saving: (on) => {
        const grid = view.grid;
        if (grid)
          ui.setUpdateIndicator(grid.root, on, 'Saving Campaign...');
      },
      waitingForLock: (on) => {
        cancelWaitNote?.();
        waiting.value = false;
        if (on)
          cancelWaitNote = host.scope.setTimeout(() => waiting.value = true, 500);
      },
      saved: (updatedAt) => {
        lastModified = updatedAt;
        needsRefresh.value = false;
      },
      openTracker,
    };

    function openTracker(): Dialog {
      if (!tracker?.isOpen.peek())
        tracker = showProgressTracker(app, host);
      return tracker;
    }

    /** The template's layout file again (read anew), else its stored layout; not saved by itself. */
    async function applyTemplateLayout(): Promise<void> {
      const template = app.template;
      if (template?.localLayoutPath) {
        try {
          view.loadLayout(DG.ViewLayout.fromJson(await grok.dapi.files.readAsText(template.localLayoutPath)));
          return;
        }
        catch (e) {
          grok.shell.error('Failed to apply layout from mounted file storage. Falling back to static layout');
          _package.logger.error(e);
        }
      }
      if (!template?.layoutViewState) {
        grok.shell.error('No layout found');
        return;
      }
      try {
        view.loadLayout(DG.ViewLayout.fromViewState(template.layoutViewState));
      }
      catch (e) {
        grok.shell.error('Failed to apply template layout. Check console for more details.');
        _package.logger.error(e);
      }
    }

    function loadLayout(viewState: string): void {
      try {
        const layout = DG.ViewLayout.fromViewState(viewState);
        addBoundColumns(layout, view.dataFrame);
        view.loadLayout(layout);
      }
      catch (e) {
        grok.shell.error('Failed to apply layout. Falling back to default layout.');
        console.error(e);
      }
      checkGrid();
    }

    function checkGrid(): void {
      if (!view.grid) {
        grok.shell.error('Applied layout created view without grid. Resetting layout.');
        view.resetLayout();
      }
    }

    /** A Stage cell is a choice of the template's stages (a row in another stage keeps it until changed). V-iDs are
     * the app's, never typed; they sort by their number (V000999 before V001000). */
    function prepareColumns(): void {
      const stages = app.stages;
      if (stages.length > 0)
        view.dataFrame.col(TileCategoriesColName)?.setTag(DG.TAGS.CHOICES, JSON.stringify(stages));
      const vid = view.dataFrame.col(ViDColName);
      if (!vid)
        return;
      vid.semType = ViDSemType;
      vid.valueComparer = (a: string | null, b: string | null) =>
        a === b ? 0 : a === null ? -1 : b === null ? 1 : parseInt(a.substring(1)) - parseInt(b.substring(1));
    }

    function openEditorOnNew(): void {
      const entry = view.dataFrame.col(options.entryColumn);
      if (entry && entry.toList().every((v) => !v))
        grok.functions.call(options.editor, {cell: view.grid.cell(options.entryColumn, 0)});
    }

    /** What the app does on the view: the address while it is current, the grid's edits and paint, the rows the user
     * adds, its menu items, a dropped file. */
    function follow(): void {
      const df = view.dataFrame;
      host.effect(() => {
        if (!signals.isCurrent.value)
          return;
        grok.shell.windows.showHelp = false;
        app.setBaseUrl();
        modifyUrl(CampaignIdKey, app.campaignId ?? '');
      });
      // a layout replaces the grid: each grid's subscriptions go with it
      host.effect(() => {
        const grid = signals.grid.value;
        if (grid === null)
          return;
        const vid = grid.col(ViDColName);
        if (vid)
          vid.editable = false;
        subscribe(grid.onCellValueEdited, (cell) => void edited(cell));
        subscribe(grid.onCellRender, (args) => paintDuplicateVid(args));
      });
      subscribe(df.onRowsAdded, () => {
        if (app.joining === 0)
          rowsAdded();
      });
      subscribe(grok.events.onContextMenu, (args) => contextMenu(args));
      subscribe(grok.events.onFileImportRequest, (args) => fileDropped(args as DG.EventData<DG.FileImportArgs>));
    }

    async function edited(cell: DG.GridCell): Promise<void> {
      const column = cell.tableColumn?.name;
      if (column === TileCategoriesColName) {
        await app.saveCampaign(false);
        return;
      }
      const value = cell.cell.value;
      const row = cell.tableRowIndex;
      if (column !== options.entryColumn || row == null)
        return;
      // the app may have opened another campaign meanwhile: this view's work ends with it
      const register = () => host.scope.isDisposed ? undefined :
        app.registerRow(row).catch((e) => console.error('Failed to register molecule', e));
      const calculate = () => host.scope.isDisposed ? undefined : app.performSingleCellCalculations(row, value);
      try {
        if (options.registerFirst) {
          if (value)
            await register();
          await calculate();
        } else {
          await calculate();
          await register();
        }
      }
      catch (e) {
        grok.shell.error('Failed to process molecule edit, check console for details');
        console.error(e);
      }
    }

    /** Every row without a stage gets the first one, and the editor opens on the last empty entry. */
    function rowsAdded(): void {
      const df = view.dataFrame;
      const stages = app.stages;
      const stageCol = df.col(TileCategoriesColName);
      if (stages.length > 0 && stageCol) {
        for (let i = 0; i < df.rowCount; i++) {
          if (!stageCol.get(i) || stageCol.isNone(i))
            stageCol.set(i, stages[0]);
        }
      }
      let empty: DG.GridCell | null = null;
      for (let i = 0; i < df.rowCount; i++) {
        const cell = view.grid.cell(options.entryColumn, i);
        if (cell && (cell.cell.value === '' || cell.cell.value === null))
          empty = cell;
      }
      if (empty)
        grok.functions.call(options.editor, {cell: empty});
    }

    function paintDuplicateVid(args: DG.GridCellRenderArgs): void {
      const cell = args.cell;
      const cache = app.duplicateVidCache;
      if (!cache || !cell?.isTableCell || cell.tableColumn?.name !== ViDColName || (cell.tableRowIndex ?? -1) < 0)
        return;
      const value = cache.indexes[cell.tableRowIndex!];
      if (cache.valueCounts[value] > 1)
        args.cell.style.backColor = DG.Color.setAlpha(DG.Color.getCategoricalColor(value), 150);
    }

    function contextMenu(args: DG.EventData): void {
      const viewer: DG.Viewer | undefined = args?.args?.context;
      const item: DG.GridCell | undefined = args?.args?.item;
      const menu: DG.Menu | undefined = args?.args?.menu;
      if (viewer?.type !== DG.VIEWER.GRID || viewer.tableView?.id !== view.id || !menu ||
        item?.tableColumn?.name !== options.entryColumn || !item.isTableCell)
        return;
      const df = view.dataFrame;
      menu.item('Add new row', () => df.rows.addNew(null, true));
      menu.item(options.duplicateLabel, () => {
        const row = df.rows.addNew(null, true);
        if (row.idx > -1)
          df.cell(row.idx, options.entryColumn).value = item.cell?.value ?? '';
      });
      const value = item.cell?.value;
      const row = item.tableRowIndex ?? -1;
      if (value && row > -1) {
        menu.item('Re-Run Calculations',
          () => void app.performSingleCellCalculations(row, value).catch((e) => console.error(e)));
      }
    }

    /** A molecule file dropped on the view: into this campaign, or a table of its own. */
    function fileDropped(args: DG.EventData<DG.FileImportArgs>): void {
      const file = args?.args?.file;
      const ext = file?.name?.split('.').pop()?.toLowerCase() ?? '';
      if (!file || !isView(grok.shell.tv, view) || !MERGE_ALLOWED_EXTENSIONS.includes(ext))
        return;
      args.preventDefault();
      void upload(file);
    }

    async function upload(file: File): Promise<void> {
      let df: DG.DataFrame | null;
      try {
        df = await readFileToDataFrame(file);
        if (!df)
          throw new Error('Failed to read file');
        df.name = file.name.slice(0, file.name.lastIndexOf('.'));
        await df.meta.detectSemanticTypes();
        await grok.data.detectSemanticTypes(df);
      }
      catch (e) {
        grok.shell.error('Hit Design Failed to read file');
        console.error(e);
        return;
      }
      const table = df;
      if (!HitAppBase.moleculeColumn(table)) {
        grok.shell.addTableView(table);
        return;
      }
      const dialog = new Dialog('Uploading table with molecules', {name: 'upload-molecules'});
      dialog.add(span('Add new table to workspace or join with existing Hit Design campaign'));
      dialog.addButton('Add to workspace', () => {
        dialog.close();
        grok.shell.addTableView(table);
      });
      dialog.onOK(() => app.handleJoiningDataframe(table), {text: 'Add to existing campaign', primary: true});
      host.show(dialog);
    }

    /** The ribbon: the package's group after the platform's icons, the crumb back to the landing, Save left off, a
     * status line for a view-only campaign and for a save waiting for another session's. */
    function chrome(): void {
      const template = app.template;
      const actions: RibbonAction[] = [
        {icon: 'sync', tooltip: 'Refresh Campaign', name: 'Refresh', run: () => app.refreshCampaign(),
          hint: computed(() => needsRefresh.value ? NEW_CHANGES : null), hintDetail: 'New changes'},
        {icon: 'plus', tooltip: 'Add new row', run: () => {view.dataFrame.rows.addNew(null, true);}},
        {icon: 'code-merge', tooltip: 'Merge another table into the campaign',
          text: 'Merge another table into the campaign...', run: () => openMergeTableDialog(app, host)},
        {icon: 'wrench', tooltip: 'Calculate additional properties', text: 'Calculate additional properties...',
          run: () => app.editCompute()},
      ];
      const designer = options.designer ? DG.Func.find({meta: {role: HitDesignerFunctionTag}})
        .find((f) => f.outputs.length === 1 && f.outputs[0].propertyType === DG.TYPE.DATA_FRAME) : undefined;
      if (designer) {
        actions.push({icon: 'pen-nib', tooltip: designer.description || designer.friendlyName,
          text: `${designer.friendlyName}...`, run: () => runDesigner(designer)});
      }
      if (app.canEdit) {
        actions.push({icon: 'users', tooltip: 'Edit campaign permissions', text: 'Edit campaign permissions...',
          run: () => editPermissions()});
      }
      if (template?.layoutViewState || template?.localLayoutPath) {
        const source = template.localLayoutPath ? 'Loaded from mounted file storage' : 'Static';
        actions.push({icon: 'svg-view-layout', tooltip: `Apply template layout (${source})`,
          run: () => applyTemplateLayout()});
      }
      const buttons: RibbonButton[] = [
        ...app.stages.length > 0 ? [{text: 'Progress tracker', run: () => {openTracker();}}] : [],
        {text: 'Submit...', run: () => app.submitView?.render(), cls: 'hit-design-submit-button'},
      ];
      host.ribbon({
        hideSave: true, name: 'hit-design-actions', actions, buttons,
        breadcrumbs: {parent: app.mainView, onClick: () => {
          view.close();
          void app.infoView.init();
        }},
        status: computed(() => !app.canEdit ? 'View only — changes are not saved' :
          waiting.value ? 'Waiting for another session to finish saving…' : ''),
      });
    }

    // pinned to the view: opening another campaign closes the view, and the dialog with it
    function editPermissions(): Promise<void> {
      const opened = app.campaign;
      const code = app.campaignId!;
      const name = opened?.friendlyName;
      return PermissionsEditor.show(opened?.permissions, {
        subject: name && name !== code ? `${code} (${name})` : code,
        author: opened?.authorUserFriendlyName, host,
        save: async (permissions) => {
          const campaign = app.campaign;
          if (!campaign)
            return false;
          const before = campaign.permissions;
          campaign.permissions = permissions;
          if (await app.saveCampaign(true) !== undefined)
            return true;
          campaign.permissions = before;
          return false;
        },
      });
    }

    /** The designer function's form, the current row's molecule given to it; what it returns joins the campaign as a
     * dropped table does. */
    function runDesigner(designer: DG.Func): void {
      const df = view.dataFrame;
      const inputs: {[name: string]: unknown} = {};
      const first = designer.inputs[0];
      if (first?.semType === DG.SEMTYPE.MOLECULE && df.currentRowIdx >= 0) {
        const molecule = _package.convertToSmiles(df.col(app.molColName)?.get(df.currentRowIdx));
        if (molecule)
          inputs[first.name] = molecule;
      }
      const call = designer.prepare(inputs);
      const dialog = new Dialog(designer.friendlyName, {name: 'designer'});
      dialog.add(() => funcForm(call));
      dialog.onOK(async () => {
        const progress = DG.TaskBarProgressIndicator.create(`Running ${designer.friendlyName}`);
        try {
          await call.call();
          const result: DG.DataFrame | undefined = call.getOutputParamValue();
          if (!result?.rowCount) {
            grok.shell.warning(`${designer.friendlyName} returned an empty result`);
            return;
          }
          await result.meta.detectSemanticTypes();
          await grok.data.detectSemanticTypes(result);
          await app.handleJoiningDataframe(result);
        }
        catch (e) {
          grok.shell.error(`Failed to run ${designer.friendlyName}`);
          _package.logger.error(e);
          return false;
        } finally {
          progress.close();
        }
      });
      host.show(dialog);
    }

    /** Every 3 s while the view is current and this session neither saves nor merges: has another session saved the
     * campaign since? The check starts once the last save is read, also when that read failed (any save it then
     * reads is news); a session closed meanwhile starts none. */
    function watchForChanges(): void {
      const campaignId = app.campaignId;
      if (!campaignId)
        return;
      host.task(async (t) => {
        try {
          const known = lastModified;
          const time = (await api.queries.getLastModified(app.appName, campaignId))
            ?.col('updated_at')?.getString(0) ?? null;
          // a save of this session that ended meanwhile told a later time
          if (lastModified === known)
            lastModified = time;
        }
        catch (e) {
          console.error(e);
        }
        t.build(() => t.owner.setInterval(() => void checkForChanges(), 3000));
      });
    }

    async function checkForChanges(): Promise<void> {
      if (checking || needsRefresh.peek() || !signals.isCurrent.peek() || app.isSaving || app.joining > 0 ||
        !app.campaignId)
        return;
      checking = true;
      try {
        const known = lastModified;
        const time = (await api.queries.getLastModified(app.appName, app.campaignId))?.col('updated_at')?.getString(0);
        // a save of this session that started or ended meanwhile logs its own time
        if (time && time !== known && lastModified === known && !app.isSaving && !host.scope.isDisposed) {
          lastModified = time;
          needsRefresh.value = true;
        }
      }
      catch (e) {
        console.error(e);
      } finally {
        checking = false;
      }
    }
  });
}

/** The Primary Molecule Column of a joined table with several: answers null when the user cancels, or when the
 * campaign view of `host`, which the dialog is pinned to, closes. */
export async function pickMoleculeColumn(names: string[], likely: string, host?: ViewHost): Promise<string | null> {
  if (host?.scope.isDisposed)
    return null;
  let chosen: string | null = null;
  const dialog = new Dialog('Select Primary Molecule Column', {name: 'primary-molecule-column'});
  dialog.build(() => {
    const input = new ChoiceInput({label: 'Primary Molecule Column', items: names, value: likely, nullable: false});
    dialog.add(span(MULTIPLE_MOLECULES));
    dialog.add(input);
    dialog.onOK(() => {
      chosen = input.peek();
    });
  });
  await new Promise<void>((resolve) => {
    dialog.onClose(() => resolve());
    host ? host.show(dialog) : dialog.show();
  });
  return chosen;
}

/** The columns the viewers of `layout` are bound to and `df` lacks yet (a new campaign's computed ones), added
 * empty and marked: those viewers start empty instead of failing on a column that does not exist, and fill in as
 * the calculations fill the columns. Which columns are bound is the platform's answer
 * (`DG.ViewLayout.boundColumnNames`: a viewer's column settings, the Filters viewer's filters); a platform without
 * it gets none added. The type is the app's guess, as the platform does not say what an absent column should hold:
 * a column the layout's grid shows as text, or a filter other than a histogram is on, is a string one
 * (`textColumns`), any other a number; a merge gives it the incoming column's type. A hidden `~` column is never
 * added. */
function addBoundColumns(layout: DG.ViewLayout, df: DG.DataFrame): void {
  const text = textColumns(layout.viewState);
  for (const name of layout.boundColumnNames ?? []) {
    if (!name.startsWith('~') && !df.columns.contains(name)) {
      df.columns.addNew(name, text.has(name) ? DG.COLUMN_TYPE.STRING : DG.COLUMN_TYPE.FLOAT)
        .setTag(PlaceholderTag, 'true');
    }
  }
}

function textColumns(viewState: string): Set<string> {
  const text = new Set<string>();
  for (const o of layoutObjects(viewState)) {
    if (typeof o.columnName === 'string' && o.cellType !== undefined && o.cellType !== 'number')
      text.add(o.columnName);
    if (o['#type'] === 'FiltersLook' && Array.isArray(o.filters)) {
      for (const filter of o.filters as {[key: string]: unknown}[]) {
        const name = filter?.column ?? filter?.columnName;
        if (typeof name === 'string' && filter.type !== DG.FILTER_TYPE.HISTOGRAM)
          text.add(name);
      }
    }
  }
  return text;
}

/** Every object of a layout's view state, each before the ones it holds. */
export function layoutObjects(viewState: string): {[key: string]: unknown}[] {
  const objects: {[key: string]: unknown}[] = [];
  const visit = (node: unknown): void => {
    if (Array.isArray(node)) {
      for (const item of node)
        visit(item);
    }
    else if (node !== null && typeof node === 'object') {
      objects.push(node as {[key: string]: unknown});
      for (const value of Object.values(node))
        visit(value);
    }
  };
  visit(JSON.parse(viewState));
  return objects;
}
