import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {APP_NAME} from '../constants';
import {Engine, hyperparametersOf, rolesOf} from '../engines/engine';
import {EngineRegistry} from '../engines/engine-registry';
import {forgeDb} from '../generated/db';
import {deleteModel, modelsChanged} from '../storage/model-store';
import {reportError} from './report-error';
import '../../css/forge.css';

const CATALOG_COLUMNS = ['name', 'engine_name', 'task', 'target_name', 'storage_mode', 'row_count'] as const;

export class ForgeApp extends DG.ViewBase {
  readonly deleteIcon: HTMLElement;
  private readonly grid: DG.Grid;
  private readonly emptyHint: HTMLElement;
  private currentRowSub: rxjs.Subscription | undefined;

  private constructor(engines: Engine[], models: DG.DataFrame) {
    super();
    this.name = APP_NAME;
    this.box = true;

    const methodsPane = ui.panel([
      ui.h2('Methods'),
      ui.table(engines, (e) => [e.name, e.namespace, e.kind, rolesOf(e).join(', '),
        hyperparametersOf(e).map((p) => p.name).join(', ')],
      ['Method', 'Package', 'Method type', 'Roles', 'Hyperparameters']),
    ]);

    this.grid = DG.Viewer.grid(models);
    this.emptyHint = ui.divText('No models yet.');
    this.deleteIcon = ui.icons.delete(() => this.deleteCurrentModel(), 'Delete model');
    const header = ui.divH([ui.h2('Models'), ui.icons.sync(() => this.refresh(), 'Refresh'), this.deleteIcon],
      'forge-pane-header');
    const modelsPane = ui.panel([header, this.emptyHint, this.grid.root]);
    this.bindCatalog();
    this.subs.push(modelsChanged.subscribe(() => this.refresh()));

    this.root.appendChild(ui.splitV([methodsPane, modelsPane]));
  }

  get models(): DG.DataFrame {
    return this.grid.dataFrame;
  }

  static async create(): Promise<ForgeApp> {
    return new ForgeApp(EngineRegistry.discover(), await ForgeApp.loadModels());
  }

  static async open(): Promise<void> {
    try {
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
      grok.dapi.domains.registry.rowProperties('forge.model'),
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

  detach(): void {
    this.currentRowSub?.unsubscribe();
    super.detach();
  }

  private async refresh(): Promise<void> {
    try {
      this.grid.dataFrame = await ForgeApp.loadModels();
      this.bindCatalog();
    } catch (e) {
      reportError(e);
    }
  }

  private bindCatalog(): void {
    const visibleColumns = [...CATALOG_COLUMNS, 'created_on'];
    this.grid.columns.setOrder(visibleColumns);
    this.grid.columns.setVisible(visibleColumns);
    ui.setDisplay(this.emptyHint, this.models.rowCount === 0);
    this.currentRowSub?.unsubscribe();
    this.currentRowSub = this.models.onCurrentRowChanged.subscribe(() => this.updateDeleteIcon());
    this.updateDeleteIcon();
  }

  private updateDeleteIcon(): void {
    ui.setDisabled(this.deleteIcon, this.models.currentRowIdx < 0);
  }

  private deleteCurrentModel(): void {
    const row = this.models.currentRowIdx;
    if (row < 0)
      return;
    const id = this.models.getCol('id').get(row);
    const name = this.models.getCol('name').get(row);
    ui.dialog('Delete model')
      .add(ui.divText(`Delete the model "${name}"?`))
      .onOK(async () => {
        try {
          await deleteModel(id);
        } catch (e) {
          reportError(e);
        }
      })
      .show();
  }
}
