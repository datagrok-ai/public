import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {providedModelOf} from '../apply/apply-model';
import {providedModels, providedRow} from '../catalog/provided-models';
import {Engine} from '../engines/engine';
import {EngineRegistry} from '../engines/engine-registry';
import {ForgeError} from '../forge-error';
import {openApplyDialog} from './apply-model-dialog';
import {Gallery} from './gallery';
import {ForgeModelHandler, forgeModelHandler} from './model-handler';
import {reportError} from './report-error';
import '../../css/forge.css';

const VISIBLE_COLUMNS = ['name', 'target_name', 'created_on'];
const RELOAD_MS = 60000;

interface ContextMenuArgs { item?: unknown; menu?: DG.Menu }

/** The models a method provides (`meta.mlrole: models`), such as the MLflow models the platform syncs: the catalog's
 * gallery over them, with **Apply...** only. The list is read again every minute and on Refresh. */
export class ProvidedModelsView extends DG.ViewBase {
  readonly gallery: Gallery;
  private readonly engine: Engine;
  private readonly grid: DG.Grid;

  private constructor(engine: Engine, models: DG.DataFrame) {
    super();
    this.engine = engine;
    this.name = models.name;
    this.box = true;
    this.grid = DG.Viewer.grid(models);
    const rowOf = (i: number) => ForgeModelHandler.rowOf(providedRow(this.models, engine, i));
    this.gallery = new Gallery({key: 'provided models', noun: 'models', placeholder: 'Search models by name or target',
      grid: this.grid, searchColumns: ['name', 'target_name'],
      sortFields: [{caption: 'Created', column: 'created_on'}, {caption: 'Name', column: 'name'},
        {caption: 'Target', column: 'target_name'}],
      defaultSort: {column: 'created_on', desc: true},
      card: (i) => ForgeModelHandler.cardContent(rowOf(i)),
      listItem: (i) => ForgeModelHandler.listLine(rowOf(i)),
      tooltip: (i) => forgeModelHandler.renderTooltip(rowOf(i)),
      menu: (i, menu) => this.fillMenu(i, menu),
      actions: ui.divH([ui.icons.sync(() => void this.refresh(), 'Refresh')]),
      emptyHint: ui.divText(`No ${engine.name} models yet.`, 'forge-gallery-empty')});
    this.showColumns();
    this.subs.push(
      rxjs.interval(RELOAD_MS).subscribe(() => void this.refresh()),
      grok.events.onContextMenu.subscribe((event: {args: ContextMenuArgs}) => this.extendMenu(event.args)));
    this.root.classList.add('forge-app');
    this.root.appendChild(this.gallery.root);
  }

  get models(): DG.DataFrame {
    return this.grid.dataFrame;
  }

  /** The view of the models the method named [method] provides. */
  static async create(method: string): Promise<ProvidedModelsView> {
    const engine = EngineRegistry.discover().find((e) => e.name === method && e.functions.models !== undefined);
    if (engine === undefined)
      throw new ForgeError(`The method '${method}' is not installed or provides no models.`);
    return new ProvidedModelsView(engine, await providedModels(engine));
  }

  detach(): void {
    this.gallery.detach();
    super.detach();
  }

  /** Reads the list again; the gallery is rebuilt only when the models changed. */
  async refresh(): Promise<void> {
    try {
      const models = await providedModels(this.engine);
      if (models.toCsv() === this.models.toCsv())
        return;
      this.grid.dataFrame = models;
      this.showColumns();
      this.gallery.rebind();
    } catch (e) {
      reportError(e);
    }
  }

  private showColumns(): void {
    this.grid.columns.setOrder(VISIBLE_COLUMNS);
    this.grid.columns.setVisible(VISIBLE_COLUMNS);
  }

  private fillMenu(i: number, menu: DG.Menu): void {
    menu.item('Apply...', () => void this.applyModel(i));
  }

  /** The Apply dialog for the model of row [i], on the current table or, outside a table view, the first open one. */
  async applyModel(i: number): Promise<void> {
    try {
      const provided = providedModelOf(providedRow(this.models, this.engine, i), this.engine);
      await openApplyDialog(grok.shell.currentTable ?? grok.shell.tables[0] ?? null, {provided, switchToTable: true});
    } catch (e) {
      reportError(e);
    }
  }

  /** Right-clicking a row of the table mode: the grid's menu gets the model's command. */
  private extendMenu({item, menu}: ContextMenuArgs): void {
    if (!(item instanceof DG.GridCell) || menu === undefined || item.grid.dart !== this.grid.dart)
      return;
    const i = item.tableRowIndex;
    if (i !== null && i >= 0)
      this.fillMenu(i, menu);
  }
}
