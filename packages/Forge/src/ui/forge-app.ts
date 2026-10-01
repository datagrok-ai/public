import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {APP_NAME} from '../constants';
import {Engine, hyperparametersOf, rolesOf} from '../engines/engine';
import {EngineRegistry} from '../engines/engine-registry';
import {ForgeError} from '../forge-error';
import {forgeDb} from '../generated/db';
import {_package} from '../package';

const CATALOG_COLUMNS = ['name', 'engine_name', 'task', 'target_name', 'storage_mode', 'row_count'] as const;

export class ForgeApp extends DG.ViewBase {
  readonly models: DG.DataFrame;

  private constructor(engines: Engine[], models: DG.DataFrame) {
    super();
    this.models = models;
    this.name = APP_NAME;
    this.box = true;

    const enginesPane = ui.panel([
      ui.h2('Engines'),
      ui.table(engines, (e) => [e.name, e.namespace, e.kind, rolesOf(e).join(', '),
        hyperparametersOf(e).map((p) => p.name).join(', ')], ['Engine', 'Package', 'Kind', 'Roles', 'Hyperparameters']),
    ]);

    const grid = DG.Viewer.grid(models);
    const visibleColumns = [...CATALOG_COLUMNS, 'created_on'];
    grid.columns.setOrder(visibleColumns);
    grid.columns.setVisible(visibleColumns);
    const hint = models.rowCount === 0 ? [ui.divText('No models yet.')] : [];
    const modelsPane = ui.panel([ui.h2('Models'), ...hint, grid.root]);

    this.root.appendChild(ui.splitV([enginesPane, modelsPane]));
  }

  static async create(): Promise<ForgeApp> {
    const engines = EngineRegistry.discover();
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
    return new ForgeApp(engines, models);
  }

  static async open(): Promise<void> {
    try {
      grok.shell.addView(await ForgeApp.create());
    } catch (e) {
      if (e instanceof ForgeError)
        grok.shell.warning(e.message);
      else {
        grok.shell.error(e instanceof Error ? e.message : String(e));
        _package.logger.error(e);
      }
    }
  }
}
