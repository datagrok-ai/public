import {DbEntity, DbEntityType} from "./entity";
import * as DG from "datagrok-api/dg";
import {DbQueryEntityCrud} from "./crud";
import * as ui from "datagrok-api/ui";
import {CruddyFilterHost} from "./filters";
import {debounceTime} from "rxjs/operators";
import {DbSchema, DbTable} from "./table";
import * as grok from "datagrok-api/grok";
import {Db} from "datagrok-api/dg";

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
  grid?: DG.Grid;
  mainViewer?: DG.Viewer;
  host: HTMLDivElement = ui.divH([]);
  filters: CruddyFilterHost = new CruddyFilterHost();
  ribbonPanel: HTMLDivElement = ui.ribbonPanel([]);

  constructor(app: CruddyApp, entityType: DbEntityType) {
    super();
    this.app = app;
    this.entityType = entityType;
    this.name = entityType.type;
    this.crud = new DbQueryEntityCrud(app.config.connection, entityType);
    this.root.style.display = 'flex';
    this.root.classList.add('cruddy-view');
    this.setRibbonPanels([[this.ribbonPanel]]);

    this.crud.read(undefined, {limit: 100}).then(async (df) => {
      this.grid = df.plot.grid();
      await grok.data.detectSemanticTypes(df);
      this.append(this.host);
      this.host.appendChild(this.filters.root);

      this.mainViewer = this.entityType.defaultView == 'cards' ? DG.Viewer.tile(df) : this.grid;
      this.host.appendChild(this.mainViewer.root);

      this.mainViewer.root.style.flexGrow = '1';
      this.mainViewer.root.style.height = 'inherit';
      this.initBehaviors();
      this.initFilters();
    });
  }

  initFilters() {
    this.filters.init(this.entityType);
    this.filters.onChanged.pipe(debounceTime(100)).subscribe((_) => {
      this.crud.read(this.filters.getCondition()).then((df) => {
        (this.mainViewer ?? this.grid)!.dataFrame = df;
      });
    });
  }

  initBehaviors() {
    CruddyViewFeature.contextDetails().attach(this);
    CruddyViewFeature.editable().attach(this);
    CruddyViewFeature.insertable().attach(this);
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