import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {requiredFeaturesOf} from '../apply/apply-model';
import {isSuggested} from '../apply/column-matching';
import {COMPARE_COLUMNS, CompareModelRow} from '../catalog/compare-models';
import {GRID_DATE_FORMAT, MODEL_TYPE, VIEW_NAME, VIEW_PATH, VIEW_TYPE} from '../constants';
import {forgeDb, ModelRow, ModelStorageMode} from '../generated/db';
import {isRecord} from '../preparation/preparation-options';
import {modelsChanged, storedModels, syncStatus} from '../storage/model-store';
import {isOpen} from './apply-model-dialog';
import {Gallery, SortField} from './gallery';
import {addModelItems} from './model-actions';
import {ModelComparison, openComparisonView} from './model-comparison';
import {ForgeModelHandler, forgeModelHandler} from './model-handler';
import {reportError} from './report-error';
import {STORAGE_CAPTIONS} from './save-model-dialog';
import {launchTraining} from './train-launcher';
import '../../css/forge.css';

const SHOWN = ['name', 'engine_name', 'task', 'target_name', 'storage_mode', 'row_count', 'tags'] as const;
// features, options and blob are not shown: the card, Applicable to and Download read them.
const CATALOG_COLUMNS = [...SHOWN, 'features', 'options', 'blob'] as const;
const VISIBLE_COLUMNS = [...SHOWN, 'created_on'];
const COMPARISON_DELAY_MS = 300;
const RELOAD_DELAY_MS = 300;
const MODEL_HEADERS: {[column: string]: string} = {
  name: 'Model name.',
  engine_name: 'The method that trained the model.',
  task: 'classification (categories) or regression (numbers).',
  target_name: 'The column the model predicts.',
  storage_mode: 'What the model keeps of its training data.',
  row_count: 'Rows the model learned from.',
  tags: 'The model\'s tags.',
  created_on: 'When the model was saved.',
};
const MODEL_SORT: SortField[] = [
  {caption: 'Created', column: 'created_on'},
  {caption: 'Name', column: 'name'},
  {caption: 'Method', column: 'engine_name'},
  {caption: 'Target', column: 'target_name'},
  {caption: 'Task', column: 'task'},
  {caption: 'Training rows', column: 'row_count'},
];

interface ContextMenuArgs { item?: unknown; menu?: DG.Menu }

export class ForgeApp extends DG.ViewBase {
  readonly applicableToInput: DG.InputBase<DG.DataFrame | null>;
  /** The view's content: the catalog as cards, a list or the grid. */
  readonly modelsGallery: Gallery;
  /** What the gallery shows while the catalog is empty. */
  readonly emptyHint: HTMLElement;
  private readonly grid: DG.Grid;
  // The current frame's subscriptions: a reload replaces the frame and drops them.
  private catalogSubs: rxjs.Subscription[] = [];
  // Numbers the comparison reads, so a slow read never replaces a newer selection's comparison.
  private comparisonRead = 0;
  // Numbers the catalog's due loads: a change or a Refresh takes the next one, and a load is shown only while its
  // number is the newest, so a slow reload never replaces a newer one's frame.
  private catalogLoad = 0;
  // The open table the catalog is filtered by.
  private filteredBy: DG.DataFrame | null = null;
  // Whether a catalog row fits the table of Applicable to; null while none is chosen.
  private fit: ((i: number) => boolean) | null = null;

  private constructor(models: DG.DataFrame) {
    super();
    this.name = VIEW_NAME;
    this.box = true;
    this.aiDescription = 'Models (address /models): the catalog of trained predictive models (the Forge package), ' +
      'with a search box (text and #tags), list, cards and table modes, a sort menu and a model count; it remembers ' +
      'its mode and sort, and shows cards by default (name, target, features, method, creation date). Applicable to ' +
      'keeps only the models whose features the chosen open table has; a clicked model is shown in the context panel ' +
      '(details, performance, activity, sharing, history), Ctrl+click selects more and two or more selected models ' +
      'are compared there; the right-click menu has Apply..., Edit model..., Share..., Download, Delete model and ' +
      'Compare (opens the selected models side by side); the toolbar\'s icon refreshes the catalog; with no ' +
      'models, Train a model... opens the Train view. From code or the AI ' +
      'assistant, the Forge functions listMethods, listModels, getModel, trainModel, applyModel, evaluateModel, ' +
      'updateModel, shareModel and deleteModel do the same without the UI (see the package\'s agents/skills.md).';

    this.grid = DG.Viewer.grid(models);
    this.grid.onCellTooltip((cell, x, y) => {
      const tooltip = this.modelTooltip(cell) ?? ForgeApp.headerTooltip(cell);
      if (tooltip !== null)
        ui.tooltip.show(tooltip, x, y);
      return tooltip !== null;
    });
    this.grid.onCellPrepare((cell) => {
      if (cell.isTableCell && cell.tableColumn?.name === 'storage_mode')
        cell.customText = ForgeApp.storageText(cell.cell.value);
    });
    this.applicableToInput = ui.input.table('Applicable to', {nullable: true,
      tooltipText: 'Show only the models whose features all have a close column in the table; empty shows every model.',
      onValueChanged: () => this.filterModels()});
    const refreshIcon = ui.icons.sync(() => {
      this.markStale();
      void this.refresh();
    }, 'Refresh');
    this.emptyHint = ui.divV([ui.divText('No models yet.'),
      ui.link('Train a model...', () => void launchTraining(null), 'Open the Train view on a table')],
    'forge-gallery-empty');
    // The background sync: the icon spins while it runs and warns of what was not refreshed, in its tooltip.
    ui.tooltip.bind(refreshIcon, () => syncStatus.value.message || 'Refresh');
    const rowOf = (i: number) => ForgeModelHandler.rowOf(this.modelValues(i));
    this.modelsGallery = new Gallery({key: 'models', noun: 'models',
      placeholder: 'Search models by name, target, method or #tag', grid: this.grid,
      searchColumns: ['name', 'target_name', 'engine_name', 'task'], tagsColumn: 'tags', sortFields: MODEL_SORT,
      defaultSort: {column: 'created_on', desc: true},
      card: (i) => ForgeModelHandler.cardContent(rowOf(i)),
      listItem: (i) => ForgeModelHandler.listLine(rowOf(i)),
      tooltip: (i) => forgeModelHandler.renderTooltip(rowOf(i)),
      menu: (i, menu) => this.fillMenu(i, menu),
      dragObject: rowOf,
      actions: ui.divH([refreshIcon]),
      below: this.applicableToInput.root,
      emptyHint: this.emptyHint,
      extraFilter: (i) => this.fit?.(i) ?? true});
    this.modelsGallery.root.classList.add('forge-catalog');
    // The gallery rendered the cards already; Applicable to starts empty.
    this.subscribeCatalog();
    this.publishLoad(0);
    this.subs.push(
      // A sync's writes reload the catalog as any change does (modelsChanged).
      syncStatus.subscribe(({checking, message}) => {
        refreshIcon.classList.toggle('fa-spin', checking);
        refreshIcon.classList.toggle('forge-card-warn', !checking && message !== '');
        // The toolbar's icon colour outweighs the class.
        refreshIcon.style.color = !checking && message !== '' ? 'var(--orange-3)' : '';
      }),
      // Stale from the change, not from the debounced reload.
      modelsChanged.subscribe(() => this.markStale()),
      // A burst of writes (tag chips) reloads once.
      DG.debounce(modelsChanged, RELOAD_DELAY_MS).subscribe(() => this.refresh()),
      // The table input drops a closed table without reporting a change.
      grok.events.onTableRemoved.subscribe(() => {
        if (this.filteredBy !== null && !isOpen(this.filteredBy))
          this.filterModels();
      }),
      grok.events.onContextMenu.subscribe((event: {args: ContextMenuArgs}) => this.extendMenu(event.args)));

    this.path = VIEW_PATH;
    this.root.classList.add('forge-app');
    this.root.appendChild(this.modelsGallery.root);
  }

  get models(): DG.DataFrame {
    return this.grid.dataFrame;
  }

  get type(): string {
    return VIEW_TYPE;
  }

  /** `/models`, in any case and with or without a trailing `/`; the base `handlePath` keeps the view as it is. */
  acceptsPath(path: string): boolean {
    return path.toLowerCase().replace(/\/$/, '') === VIEW_PATH;
  }

  static async create(): Promise<ForgeApp> {
    return new ForgeApp(await ForgeApp.loadModels());
  }

  /** Focuses an open view, or opens one. */
  static async open(): Promise<void> {
    try {
      const open = ForgeApp.openView();
      if (open === undefined) {
        grok.shell.addView(await ForgeApp.create());
        return;
      }
      grok.shell.v = open;
    } catch (e) {
      reportError(e);
    }
  }

  private static openView(): ForgeApp | undefined {
    for (const view of grok.shell.views) {
      if (view instanceof ForgeApp)
        return view;
    }
    return undefined;
  }

  static async loadModels(): Promise<DG.DataFrame> {
    const [models, properties] = await Promise.all([
      storedModels()
        .select(...CATALOG_COLUMNS)
        .orderBy('created_on', true)
        .df(),
      grok.dapi.domains.registry.rowProperties(MODEL_TYPE),
    ]);
    for (const p of properties) {
      const col = models.col(p.name);
      if (col !== null)
        col.meta.friendlyName = p.friendlyName;
    }
    const createdOn = models.col('created_on');
    if (createdOn !== null) {
      createdOn.meta.friendlyName = 'Created';
      // The cards' date in the table mode; the grid still sorts by the full datetime.
      createdOn.meta.format = GRID_DATE_FORMAT;
    }
    return models;
  }

  /** Opens the selected models side by side in the table view **Compare models**. */
  async compareSelected(): Promise<void> {
    try {
      openComparisonView(await this.selectedModels());
    } catch (e) {
      reportError(e);
    }
  }

  detach(): void {
    this.unsubscribeCatalog();
    this.modelsGallery.detach();
    // A comparison read still running finds itself outdated and sets nothing.
    this.comparisonRead++;
    super.detach();
  }

  /** A newer load is due: a load still running is superseded and shows nothing, and the pane is busy until the due load
   * is shown. */
  private markStale(): void {
    this.catalogLoad++;
    this.modelsGallery.root.setAttribute('aria-busy', 'true');
  }

  /** Loads the catalog under the due load's number ({@link markStale} first); the current row and the selected rows
   * that still exist are found again by id. A context panel left on a comparison with a model the catalog no longer has
   * (a delete) shows the selection again (see `showSelection`). */
  private async refresh(): Promise<void> {
    const load = this.catalogLoad;
    try {
      const models = await ForgeApp.loadModels();
      if (load !== this.catalogLoad)
        return;
      const ids = this.models.getCol('id');
      const current = this.models.currentRowIdx;
      const currentId: unknown = current < 0 ? null : ids.get(current);
      const selectedIds = new Set(Array.from(this.models.selection.getSelectedIndexes(), (i) => ids.get(i)));
      this.grid.dataFrame = models;
      this.bindCatalog();
      const reloaded = this.models.getCol('id');
      if (selectedIds.size > 0)
        this.models.selection.init((i) => selectedIds.has(reloaded.get(i)));
      const restored = currentId === null ? -1 : reloaded.toList().indexOf(currentId);
      if (restored >= 0)
        this.models.currentRowIdx = restored;
      this.publishLoad(load);
      // A restored selection reports its change, which shows the selection again; without one it is checked here.
      if (this.models.selection.trueCount === 0 && this.isStaleComparison(grok.shell.o))
        await this.showSelection();
    } catch (e) {
      if (load === this.catalogLoad)
        this.modelsGallery.root.setAttribute('aria-busy', 'false');
      reportError(e);
    }
  }

  /** The catalog's readiness for automation (ARCHITECTURE.md, "Automation surface"): load [load] is shown. */
  private publishLoad(load: number): void {
    this.modelsGallery.root.setAttribute('data-forge-catalog-load', String(load));
    this.modelsGallery.root.setAttribute('aria-busy', 'false');
  }

  /** Whether [x] is a comparison with a model that the catalog no longer has; a single model that is gone leaves the
   * panel through `ForgeModelHandler.followChanges`. */
  private isStaleComparison(x: unknown): boolean {
    if (!(x instanceof ModelComparison))
      return false;
    const ids = new Set<unknown>(this.models.getCol('id').toList());
    return x.rows.some((row) => !ids.has(row.id));
  }

  /** Follows a reloaded frame: the subscriptions, Applicable to and the gallery. */
  private bindCatalog(): void {
    this.subscribeCatalog();
    this.updateFit();
    this.modelsGallery.rebind();
  }

  private subscribeCatalog(): void {
    this.grid.columns.setOrder(VISIBLE_COLUMNS);
    this.grid.columns.setVisible(VISIBLE_COLUMNS);
    this.unsubscribeCatalog();
    this.catalogSubs = [
      this.models.onCurrentRowChanged.subscribe(() => this.showCurrentModel()),
      // A selection is often made in several steps: the comparison is read once it settles.
      DG.debounce(this.models.onSelectionChanged, COMPARISON_DELAY_MS).subscribe(() => void this.showSelection()),
    ];
  }

  private unsubscribeCatalog(): void {
    for (const sub of this.catalogSubs)
      sub.unsubscribe();
    this.catalogSubs = [];
  }

  /** The chosen model becomes the current object, shown in the context panel; no current row, or two or more selected
   * rows (the panel compares them), leave the panel as it is. */
  private showCurrentModel(): void {
    if (this.models.selection.trueCount < 2)
      this.showModel(this.models.currentRowIdx);
  }

  /** Shows the model of catalog row [i] (none for -1) in the context panel; the model the panel already shows (a row
   * found again after a reload) is not set again. */
  private showModel(i: number): void {
    const shown: unknown = grok.shell.o;
    if (i >= 0 && !(shown instanceof DG.DomainRow && shown.id === this.models.get('id', i)))
      // Forced: a plain set within a second of the previous one (a click after a click) is ignored.
      grok.shell.setCurrentObject(ForgeModelHandler.rowOf(this.modelValues(i)), true, true);
  }

  /** Two or more selected rows: the context panel compares them. Fewer: it shows the current row's model, without a
   * current row the selected row's, and with neither it is emptied rather than left on a comparison. */
  private async showSelection(): Promise<void> {
    const read = ++this.comparisonRead;
    const selection = this.models.selection;
    if (selection.trueCount < 2) {
      const current = this.models.currentRowIdx;
      const i = current >= 0 ? current : selection.trueCount === 1 ? selection.getSelectedIndexes()[0] : -1;
      if (i >= 0)
        this.showModel(i);
      else if (grok.shell.o instanceof ModelComparison)
        grok.shell.setCurrentObject(null, true, true);
      return;
    }
    try {
      const rows = await this.selectedModels();
      if (read === this.comparisonRead && rows.length >= 2)
        grok.shell.setCurrentObject(new ModelComparison(rows), true, true);
    } catch (e) {
      reportError(e);
    }
  }

  /** The selected models with the columns a comparison reads, in the catalog's order. */
  private async selectedModels(): Promise<CompareModelRow[]> {
    const ids: string[] = Array.from(this.models.selection.getSelectedIndexes(), (i) => this.models.get('id', i));
    const rows: CompareModelRow[] = await forgeDb.models.query().where('id', '=', ids).select(...COMPARE_COLUMNS)
      .top(ids.length);
    const byId = new Map(rows.map((r) => [r.id, r]));
    return ids.map((id) => byId.get(id)).filter((r): r is CompareModelRow => r !== undefined);
  }

  /** With a table chosen in **Applicable to**, only the models it fits stay in the catalog; nothing leaves the
   * browser. */
  private filterModels(): void {
    this.updateFit();
    this.modelsGallery.refilter();
  }

  private updateFit(): void {
    const chosen = this.applicableToInput.value;
    const table = chosen !== null && isOpen(chosen) ? chosen : null;
    this.filteredBy = table;
    const features = this.models.col('features');
    const options = this.models.col('options');
    // Retrained models share their feature lists: each list is matched against the table once.
    const fitting = new Map<string, boolean>();
    this.fit = table === null ? null : (i: number): boolean => {
      const required = requiredFeaturesOf({features: jsonOf(features?.get(i)), options: jsonOf(options?.get(i))});
      const key = JSON.stringify(required);
      let fit = fitting.get(key);
      if (fit === undefined) {
        fit = required !== null && isSuggested(required, table);
        fitting.set(key, fit);
      }
      return fit;
    };
  }

  /** The commands of catalog row [i], and Compare for a selection of two or more: the gallery's and the grid's menu. */
  private fillMenu(i: number, menu: DG.Menu): void {
    addModelItems(menu, ForgeModelHandler.rowOf(this.modelValues(i)), this.applicableToInput.value);
    if (this.models.selection.trueCount >= 2)
      menu.separator().item('Compare', () => void this.compareSelected());
  }

  /** Right-clicking a row of the table mode: the grid's menu gets the model's commands. */
  private extendMenu({item, menu}: ContextMenuArgs): void {
    if (!(item instanceof DG.GridCell) || menu === undefined || item.grid.dart !== this.grid.dart)
      return;
    const i = item.tableRowIndex;
    if (i !== null && i >= 0)
      this.fillMenu(i, menu);
  }

  /** The model's tooltip for a **Name** [cell]; null for any other cell, which keeps the grid's own tooltip. */
  modelTooltip(cell: DG.GridCell): HTMLElement | null {
    if (!cell.isTableCell || cell.tableColumn?.name !== 'name')
      return null;
    const i = cell.tableRowIndex;
    return i === null || i < 0 ? null : forgeModelHandler.renderTooltip(ForgeModelHandler.rowOf(this.modelValues(i)));
  }

  /** A **Data storage** value as the tooltips and the model panel show it: `Reference`, `None`, `Copy`. */
  static storageText(value: unknown): string {
    return STORAGE_CAPTIONS[value as ModelStorageMode] ?? String(value ?? '');
  }

  /** The column's meaning for a header [cell]; null for any other cell. */
  private static headerTooltip(cell: DG.GridCell): HTMLElement | null {
    const text = cell.isColHeader ? MODEL_HEADERS[cell.tableColumn?.name ?? ''] : undefined;
    return text === undefined ? null : ui.divText(text);
  }

  /** The catalog row [i] as model values, its features and options as objects. */
  private modelValues(i: number): Pick<ModelRow, 'id' | 'name'> & Partial<ModelRow> {
    const values: {[column: string]: unknown} = {};
    for (const col of this.models.columns.toList())
      values[col.name] = col.isNone(i) ? undefined : col.get(i);
    return {...values, id: this.models.get('id', i), name: this.models.get('name', i),
      features: jsonOf(values.features), options: jsonOf(values.options)};
  }
}

/** A json column's value as an object: the typed frame may carry it as text. */
function jsonOf(value: unknown): {[key: string]: unknown} | undefined {
  const parsed: unknown = typeof value === 'string' && value.startsWith('{') ? JSON.parse(value) : value;
  return isRecord(parsed) ? parsed : undefined;
}
