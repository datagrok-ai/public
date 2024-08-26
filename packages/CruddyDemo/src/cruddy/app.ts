import {DbEntity, DbEntityType} from "./entity";
import * as DG from "datagrok-api/dg";
import {DbQueryEntityCrud, IQueryOptions} from "./crud";
import * as ui from "datagrok-api/ui";
import {CruddyFilterHost} from "./filters";
import {debounceTime} from "rxjs/operators";
import {DbSchema, DbTable} from "./table";
import * as grok from "datagrok-api/grok";
import {Db} from "datagrok-api/dg";
import {Subject} from "rxjs";

export class CruddyViewConfig {
  entityType: DbEntityType;

  constructor(entityType: DbEntityType) {
    this.entityType = entityType;
  }
}

/** Query, browse, navigate, edit entities */
export class CruddyEntityView<TEntity extends DbEntityType = DbEntityType> extends DG.ViewBase {
  app: CruddyApp;
  entityType: DbEntityType;
  crud: DbQueryEntityCrud;
  dataFrame?: DG.DataFrame;
  grid?: DG.Grid;
  mainViewer?: DG.Viewer;
  host: HTMLDivElement = ui.divH([]);
  filters: CruddyFilterHost = new CruddyFilterHost();
  ribbonPanel: HTMLDivElement = ui.ribbonPanel([]);
  queryOptions: IQueryOptions = { offset: 0, limit: 100 };
  totalRowCount?: number;

  onInitialized: Subject<any> = new Subject<any>();
  onQueryStarted: Subject<any> = new Subject<any>();
  onQueryCompleted: Subject<any> = new Subject<any>();
  onFilterRowsCounted: Subject<any> = new Subject<any>();
  onFilterChanged: Subject<any> = new Subject<any>();

  constructor(app: CruddyApp, entityType: DbEntityType) {
    super();
    this.app = app;
    this.entityType = entityType;
    this.name = entityType.type;
    this.crud = new DbQueryEntityCrud(app.config.connection, entityType);
    this.root.style.display = 'flex';
    this.root.classList.add('cruddy-view');
    this.setRibbonPanels([[this.ribbonPanel]]);

    this.append(this.host);
    this.host.appendChild(this.filters.root);

    this.initBehaviors();
    this.initFilters();

    this.refresh().then((_) => this.refreshCount());
  }

  async refresh(): Promise<DG.DataFrame> {
    this.onQueryStarted.next(null);
    this.dataFrame = await this.crud.read(this.filters.getCondition(), this.queryOptions);
    await grok.data.detectSemanticTypes(this.dataFrame);

    if (!this.grid) {
      this.grid = this.dataFrame.plot.grid();
      this.mainViewer = this.entityType.defaultView == 'cards' ? DG.Viewer.tile(this.dataFrame) : this.grid;
      this.mainViewer.root.style.flexGrow = '1';
      this.mainViewer.root.style.height = 'inherit';

      this.host.appendChild(this.mainViewer.root);
      this.onInitialized.next(null);
    }

    (this.mainViewer ?? this.grid)!.dataFrame = this.dataFrame;
    this.onQueryCompleted.next(null);
    return this.dataFrame;
  }

  async refreshCount(): Promise<any> {
    this.totalRowCount = await this.crud.count(this.filters.getCondition(), this.queryOptions);
    this.onFilterRowsCounted.next(null);
  }

  initFilters() {
    this.filters.init(this.entityType);
    this.filters.onChanged.pipe(debounceTime(100)).subscribe(async (_) => {
      await this.refresh();
      await this.refreshCount();
    });
  }

  initBehaviors() {
    const features = [
      CruddyViewFeature.contextDetails(),
      CruddyViewFeature.editable(),
      CruddyViewFeature.insertable(),
      CruddyViewFeature.navigationBar()
    ]
    this.onInitialized.subscribe((_) => features.forEach((f) => f.attach(this)));
  }
}

/** Features can be dynamically added to the view */
export class CruddyViewFeature {
  name: string;
  attach: (view: CruddyEntityView) => void;

  constructor(name: string, attach: (view: CruddyEntityView) => void) {
    this.name = name;
    this.attach = attach;
  }

  static navigationBar(): CruddyViewFeature {
    return new CruddyViewFeature('navigationBar', (v) => {

      let divRange = ui.divText(``);

      function q(offset: number, limit?: number) {
        v.queryOptions.offset = offset;
        v.queryOptions.limit = limit ?? v.queryOptions.limit;
        v.refresh().then((df) => {
          divRange.innerHTML = `${v.queryOptions.offset! + 1} - ${v.queryOptions.offset! + df.rowCount}${v.totalRowCount != null ? ` of ${v.totalRowCount}` : ''}`;
        });
      }

      const left = ui.iconFA('chevron-left', (_) => q(Math.max(0, v.queryOptions.offset! - v.queryOptions.limit!)));
      const right = ui.iconFA('chevron-right', (_) => q(Math.max(0, v.queryOptions.offset! + v.queryOptions.limit!)));
      v.setRibbonPanels([...v.getRibbonPanels(), ...[[divRange, left, right]]]);

      v.onQueryCompleted.subscribe()
    });
  }

  /** Shows information for the referenced rows in the context panel.
   * Referenced row: shows all details
   * Rows that reference this row: in a grid */
  static contextDetails(): CruddyViewFeature {
    return new CruddyViewFeature('context details', (v) => {
      v.grid!.onCurrentCellChanged.pipe(debounceTime(500)).subscribe((gridCell: DG.GridCell) => {
        let row = gridCell.tableRow!;

        // tables that this row references
        let acc = ui.accordion(`entity-details-${this.name}`);
        for (let col of v.entityType.columns.filter((c) => c.references)) {
          const detailsEntity = v.app.config.getEntityType(col.references!.table);
          detailsEntity.crud
            .readSingle({[col.references!.name]: row.get(col.name)})
            .then((values) => {
              if (values)
                acc.addPane(detailsEntity.type, () => ui.tableFromMap(values));
            });
        }

        // references to this row
        let referenced = v.entityType.table.columns.filter((c) => c.referencedBy.length > 0);
        if (referenced.length != 0) {
          for (let col of referenced)
            for (let c2 of col.referencedBy) {
              let host = ui.div();
              acc.addPane(c2.table.name, () => host);

              let detailsEntityType = v.app.config.getEntityType(c2.table);
              detailsEntityType.crud
                .read({[col.name]: row.get(col.name)})
                .then((df) => host.appendChild(df.plot.grid().root));
            }
        }

        acc.end();
        grok.shell.o = acc.root;
      });
    });
  }

  /** Saves edited cells to the database immediately */
  static editable(): CruddyViewFeature {
    return new CruddyViewFeature('editable', (v) => {
      v.grid!.onCellValueEdited.subscribe((gridCell: DG.GridCell) => {
        //let crud = new DbQueryEntityCrud(v.app.config.connection, v.entityType);
        v.entityType.crud
          .update(v.entityType.rowToEntity(gridCell.tableRow!))
          .then((_) => grok.shell.info('Saved'));
      });
    });
  }

  /** Add the '+' icon on top that lets you insert new row in the corresponding table */
  static insertable(): CruddyViewFeature {
    return new CruddyViewFeature('insertable', (v) => {
      v.ribbonPanel.appendChild(ui.iconFA('plus', (_) => {
        const entity = DbEntity.createNew(v.entityType);
        ui.dialog({title: `Add ${v.entityType.type}`})
          .add(ui.input.form(entity, entity.entityType.props))
          .onOK(() => {
            v.entityType.crud.create(entity).then((_) => grok.shell.info('Created'));
          })
          .show();
      }));
    });
  }
}

export class CruddyConfig {
  connection: string = '';
  schema: DbSchema = new DbSchema('', []);
  entityTypes: DbEntityType[] = [];

  getTable(name: string): DbTable {
    return this.schema.tables.find((t) => t.name == name)!;
  }

  getEntityType(table: DbTable): DbEntityType {
    return this.entityTypes.find((et) => et.table == table)!;
  }

  constructor(init?: Partial<CruddyConfig>) {
    Object.assign(this, init);
    for (let e of this.entityTypes) {
      e.crud.connectionId = this.connection;
    }
  }
}

/** Represents a running app. */
export class CruddyApp {
  config: CruddyConfig;
  views: CruddyEntityView[] = [];

  constructor(config: CruddyConfig, views?: CruddyEntityView[]) {
    this.config = config;
    this.views = views ?? this.views;
  }

  run(): void {
    for (let entityType of this.config.entityTypes) {
      let view = new CruddyEntityView(this, entityType);
      grok.shell.addView(view);
      this.views.push(view);
    }
  }
}