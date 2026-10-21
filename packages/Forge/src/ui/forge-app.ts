import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {requiredFeaturesOf} from '../apply/apply-model';
import {isSuggested} from '../apply/column-matching';
import {COMPARE_COLUMNS, CompareModelRow} from '../catalog/compare-models';
import {APP_NAME, MODEL_TYPE} from '../constants';
import {Engine, EngineKind, EngineRole, hyperparametersOf, rolesOf} from '../engines/engine';
import {EngineRegistry} from '../engines/engine-registry';
import {methodsFrame} from '../engines/method-list';
import {forgeDb, ModelRow} from '../generated/db';
import {isRecord} from '../preparation/preparation-options';
import {modelsChanged} from '../storage/model-store';
import {isOpen} from './apply-model-dialog';
import {readOnlyGrid} from './data-grid';
import {addModelItems, confirmDeleteModel, openModelApply} from './model-actions';
import {ModelComparison, openComparisonView} from './model-comparison';
import {ForgeModelHandler, forgeModelHandler} from './model-handler';
import {reportError} from './report-error';
import '../../css/forge.css';

const SHOWN = ['name', 'engine_name', 'task', 'target_name', 'storage_mode', 'row_count', 'tags'] as const;
// features, options and blob are not shown: the card, Applicable to and Download read them.
const CATALOG_COLUMNS = [...SHOWN, 'features', 'options', 'blob'] as const;
const VISIBLE_COLUMNS = [...SHOWN, 'created_on'];
const COMPARISON_DELAY_MS = 300;
const RELOAD_DELAY_MS = 300;
const METHOD_HEADERS: {[column: string]: string} = {
  'Method': 'Machine learning method, as its package names it.',
  'Package': 'Package or script namespace that provides the method.',
  'Method type': 'How the method is provided.',
  'Roles': 'What the method can do.',
  'Hyperparameters': 'Settings of the method\'s train function.',
};
const KIND_MEANINGS: Record<EngineKind, string> = {
  function: 'A package function.',
  script: 'A script with the method\'s roles.',
};
const ROLE_MEANINGS: Record<EngineRole, string> = {
  train: 'Trains a model on a table and a target column, returns the model',
  apply: 'Applies a trained model to a table, returns the predictions',
  isApplicable: 'Tells whether the method can learn from the given data',
  isInteractive: 'Tells whether training is fast enough to rerun on every change',
  visualize: 'Shows method-specific views of a trained model',
};

interface ContextMenuArgs { item?: unknown; menu?: DG.Menu }

export class ForgeApp extends DG.ViewBase {
  readonly applyIcon: HTMLElement;
  readonly deleteIcon: HTMLElement;
  readonly compareIcon: HTMLElement;
  readonly applicableToInput: DG.InputBase<DG.DataFrame | null>;
  readonly methodsGrid: DG.Grid;
  private readonly grid: DG.Grid;
  private readonly catalogPane: HTMLElement;
  private readonly emptyHint: HTMLElement;
  // The current frame's subscriptions: a reload replaces the frame and drops them.
  private catalogSubs: rxjs.Subscription[] = [];
  // Numbers the comparison reads, so a slow read never replaces a newer selection's comparison.
  private comparisonRead = 0;
  // Numbers the catalog's due loads: a change or a Refresh takes the next one, and a load is shown only while its
  // number is the newest, so a slow reload never replaces a newer one's frame.
  private catalogLoad = 0;
  // The open table the grid is filtered by.
  private filteredBy: DG.DataFrame | null = null;

  private constructor(engines: Engine[], models: DG.DataFrame) {
    super();
    this.name = APP_NAME;
    this.box = true;
    this.aiDescription = 'Forge: predictive modeling. Methods lists the installed machine learning methods with ' +
      'their roles and settings. Models is the catalog of trained models: name, method, task, target, data storage, ' +
      'training rows, tags and creation date; Applicable to keeps only the models whose features the chosen open ' +
      'table has; a chosen model is shown in the context panel (details, performance, activity, sharing, history), ' +
      'the play icon applies it to a table, the trash icon deletes it, Compare opens two or more selected models ' +
      'side by side. From code or the AI assistant, the Forge functions listMethods, listModels, getModel, ' +
      'trainModel, applyModel, evaluateModel, updateModel, shareModel and deleteModel do the same without the UI ' +
      '(see the package\'s agents/skills.md).';

    this.methodsGrid = readOnlyGrid(methodsFrame(engines), {headers: METHOD_HEADERS, cells: {
      'Method type': (row) => KIND_MEANINGS[engines[row].kind],
      'Roles': (row) => rolesOf(engines[row]).map((r) => `${r}: ${ROLE_MEANINGS[r]}`).join('\n'),
      'Hyperparameters': (row) => hyperparametersOf(engines[row]).map((p) => `${p.name}: ${p.description}`).join('\n'),
    }});
    const methodsPane = ui.panel([ui.h2('Methods'), this.methodsGrid.root]);

    this.grid = DG.Viewer.grid(models);
    // Inline: the platform's 400px width for a box in a panel (ui.css) gives way only to an element style.
    this.grid.root.style.width = '100%';
    this.grid.onCellTooltip((cell, x, y) => {
      const tooltip = this.modelTooltip(cell);
      if (tooltip !== null)
        ui.tooltip.show(tooltip, x, y);
      return tooltip !== null;
    });
    this.emptyHint = ui.divText('No models yet.');
    this.applyIcon = ui.iconFA('play', () => this.applyCurrentModel(), 'Apply model');
    this.deleteIcon = ui.icons.delete(() => this.deleteCurrentModel(), 'Delete model');
    this.compareIcon = ui.iconFA('columns', () => this.compareSelected(), 'Compare in a new view');
    this.applicableToInput = ui.input.table('Applicable to', {nullable: true,
      tooltipText: 'Show only the models whose features all have a close column in the table; empty shows every model.',
      onValueChanged: () => this.filterModels()});
    const refreshIcon = ui.icons.sync(() => {
      this.markStale();
      void this.refresh();
    }, 'Refresh');
    const header = ui.divH([ui.h2('Models'), refreshIcon, this.applyIcon, this.deleteIcon, this.compareIcon],
      'forge-pane-header');
    this.catalogPane = ui.panel([header, this.applicableToInput.root, this.emptyHint, this.grid.root], 'forge-catalog');
    this.bindCatalog();
    this.publishLoad(0);
    this.subs.push(
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

    this.root.classList.add('forge-app');
    this.root.appendChild(ui.splitV([methodsPane, this.catalogPane]));
  }

  get models(): DG.DataFrame {
    return this.grid.dataFrame;
  }

  static async create(): Promise<ForgeApp> {
    return new ForgeApp(EngineRegistry.discover(), await ForgeApp.loadModels());
  }

  /** Focuses an open catalog, or opens one. */
  static async open(): Promise<void> {
    try {
      const open = Array.from(grok.shell.views).find((v) => v instanceof ForgeApp);
      if (open !== undefined)
        grok.shell.v = open;
      else
        grok.shell.addView(await ForgeApp.create());
    } catch (e) {
      reportError(e);
    }
  }

  static async loadModels(): Promise<DG.DataFrame> {
    const [models, properties] = await Promise.all([
      forgeDb.models.query()
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
    if (createdOn !== null)
      createdOn.meta.friendlyName = 'Created';
    return models;
  }

  /** Opens the models selected in the grid side by side in the table view **Compare models**. */
  async compareSelected(): Promise<void> {
    try {
      openComparisonView(await this.selectedModels());
    } catch (e) {
      reportError(e);
    }
  }

  detach(): void {
    this.unsubscribeCatalog();
    // A comparison read still running finds itself outdated and sets nothing.
    this.comparisonRead++;
    super.detach();
  }

  /** A newer load is due: a load still running is superseded and shows nothing, and the pane is busy until the due load
   * is shown. */
  private markStale(): void {
    this.catalogLoad++;
    this.catalogPane.setAttribute('aria-busy', 'true');
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
        this.catalogPane.setAttribute('aria-busy', 'false');
      reportError(e);
    }
  }

  /** The catalog's readiness for automation (ARCHITECTURE.md, "Automation surface"): load [load] is shown. */
  private publishLoad(load: number): void {
    this.catalogPane.setAttribute('data-forge-catalog-load', String(load));
    this.catalogPane.setAttribute('aria-busy', 'false');
  }

  /** Whether [x] is a comparison with a model that the catalog no longer has; a single model that is gone leaves the
   * panel through `ForgeModelHandler.followChanges`. */
  private isStaleComparison(x: unknown): boolean {
    if (!(x instanceof ModelComparison))
      return false;
    const ids = new Set<unknown>(this.models.getCol('id').toList());
    return x.rows.some((row) => !ids.has(row.id));
  }

  private bindCatalog(): void {
    this.grid.columns.setOrder(VISIBLE_COLUMNS);
    this.grid.columns.setVisible(VISIBLE_COLUMNS);
    ui.setDisplay(this.emptyHint, this.models.rowCount === 0);
    this.unsubscribeCatalog();
    this.catalogSubs = [
      this.models.onCurrentRowChanged.subscribe(() => this.showCurrentModel()),
      this.models.onSelectionChanged.subscribe(() => this.updateCompareIcon()),
      // A selection is often made in several steps: the comparison is read once it settles.
      DG.debounce(this.models.onSelectionChanged, COMPARISON_DELAY_MS).subscribe(() => void this.showSelection()),
    ];
    this.updateRowIcons();
    this.updateCompareIcon();
    this.filterModels();
  }

  private unsubscribeCatalog(): void {
    for (const sub of this.catalogSubs)
      sub.unsubscribe();
    this.catalogSubs = [];
  }

  /** The chosen model becomes the current object, shown in the context panel; no current row, or two or more selected
   * rows (the panel compares them), leave the panel as it is. */
  private showCurrentModel(): void {
    this.updateRowIcons();
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

  /** The selected models with the columns a comparison reads, in the grid's order. */
  private async selectedModels(): Promise<CompareModelRow[]> {
    const ids: string[] = Array.from(this.models.selection.getSelectedIndexes(), (i) => this.models.get('id', i));
    const rows: CompareModelRow[] = await forgeDb.models.query().where('id', '=', ids).select(...COMPARE_COLUMNS)
      .top(ids.length);
    const byId = new Map(rows.map((r) => [r.id, r]));
    return ids.map((id) => byId.get(id)).filter((r): r is CompareModelRow => r !== undefined);
  }

  private updateRowIcons(): void {
    const hasRow = this.models.currentRowIdx >= 0;
    ui.setDisabled(this.applyIcon, !hasRow);
    ui.setDisabled(this.deleteIcon, !hasRow);
  }

  private updateCompareIcon(): void {
    ui.setDisabled(this.compareIcon, this.models.selection.trueCount < 2);
  }

  /** With a table chosen in **Applicable to**, only the models it fits stay in the grid; nothing leaves the browser. */
  private filterModels(): void {
    const chosen = this.applicableToInput.value;
    const table = chosen !== null && isOpen(chosen) ? chosen : null;
    this.filteredBy = table;
    const features = this.models.col('features');
    const options = this.models.col('options');
    // Retrained models share their feature lists: each list is matched against the table once.
    const fitting = new Map<string, boolean>();
    const fits = (i: number, open: DG.DataFrame): boolean => {
      const required = requiredFeaturesOf({features: jsonOf(features?.get(i)), options: jsonOf(options?.get(i))});
      const key = JSON.stringify(required);
      let fit = fitting.get(key);
      if (fit === undefined) {
        fit = required !== null && isSuggested(required, open);
        fitting.set(key, fit);
      }
      return fit;
    };
    this.models.filter.init((i) => table === null || fits(i, table));
  }

  /** Right-clicking a catalog row: the grid's menu gets the model's commands, and Compare for a selection of two or
   * more. */
  private extendMenu({item, menu}: ContextMenuArgs): void {
    if (!(item instanceof DG.GridCell) || menu === undefined || item.grid.dart !== this.grid.dart)
      return;
    const i = item.tableRowIndex;
    if (i !== null && i >= 0)
      addModelItems(menu, ForgeModelHandler.rowOf(this.modelValues(i)), this.applicableToInput.value);
    if (this.models.selection.trueCount >= 2)
      menu.item('Compare', () => void this.compareSelected());
  }

  /** The model's tooltip for a **Name** [cell]; null for any other cell, which keeps the grid's own tooltip. */
  modelTooltip(cell: DG.GridCell): HTMLElement | null {
    if (!cell.isTableCell || cell.tableColumn?.name !== 'name')
      return null;
    const i = cell.tableRowIndex;
    return i === null || i < 0 ? null : forgeModelHandler.renderTooltip(ForgeModelHandler.rowOf(this.modelValues(i)));
  }

  /** The catalog row [i] as model values, its features and options as objects. */
  private modelValues(i: number): Pick<ModelRow, 'id' | 'name'> & Partial<ModelRow> {
    const values: {[column: string]: unknown} = {};
    for (const col of this.models.columns.toList())
      values[col.name] = col.isNone(i) ? undefined : col.get(i);
    return {...values, id: this.models.get('id', i), name: this.models.get('name', i),
      features: jsonOf(values.features), options: jsonOf(values.options)};
  }

  private async applyCurrentModel(): Promise<void> {
    const row = this.models.currentRowIdx;
    if (row >= 0)
      await openModelApply(this.models.getCol('id').get(row), this.applicableToInput.value);
  }

  private deleteCurrentModel(): void {
    const row = this.models.currentRowIdx;
    if (row >= 0)
      confirmDeleteModel(this.models.getCol('id').get(row), this.models.getCol('name').get(row));
  }
}

/** A json column's value as an object: the typed frame may carry it as text. */
function jsonOf(value: unknown): {[key: string]: unknown} | undefined {
  const parsed: unknown = typeof value === 'string' && value.startsWith('{') ? JSON.parse(value) : value;
  return isRecord(parsed) ? parsed : undefined;
}
