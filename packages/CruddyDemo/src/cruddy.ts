import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {debounceTime} from 'rxjs/operators';


export class DbTable {
  schema: string = '';
  name: string = '';
  columns: DbColumn[] = [];

  getColumn(name: string): DbColumn {
    return this.columns.find((c) => c.name == name)!;
  }

  constructor(init?: Partial<DbTable>) {
    Object.assign(this, init);
    for (let col of this.columns)
      col.table = this;
  }
}


/** Column belongs to a table. */
export class DbColumn {
  table: DbTable = new DbTable;
  name: string = '';
  type: string = '';
  isKey: boolean = false;
  ref: string = '';
  references?: DbColumn;
  referencedBy: DbColumn[] = [];
  prop?: DG.Property;  // additional properties such as semantic type, editor, formatting, etc.

  constructor(init?: Partial<DbColumn>) {
    Object.assign(this, init);
  }
}


/** DbEntity is defined by a list of columns */
export class DbEntityType {
  type: string = '';
  table: DbTable = new DbTable();
  behaviors: string[] = [];

  constructor(init?: Partial<DbEntityType>) {
    Object.assign(this, init);
  }

  rowToEntity(row: DG.Row): DbEntity {
    return new DbEntity(this, Object.fromEntries(this.columns.map(c => [c.name, row.get(c.name)])));
  }

  get columns(): DbColumn[] { return this.table.columns; };
  getColumn(name: string): DbColumn { return this.columns.find(c => c.name === name)!; };
}


export class DbEntity<T extends DbEntityType = DbEntityType> {
  entityType: T;
  values: {[key: string]: any };

  constructor(entityType: T, values: { [p: string]: any }) {
    this.entityType = entityType;
    this.values = values;
  }

  columns(options?: {isKey: boolean}): DbColumn[] {
    return this.entityType.columns
      .filter(c => options?.isKey == undefined || options.isKey == c.isKey);
  }
}


export abstract class EntityCrud<TEntity extends DbEntity> {

  /** Creates the entity in the database.
   * Returns the entity on success, or throws an exception otherwise. */
  abstract create(entity: TEntity): Promise<TEntity>;

  /** Reads the entities, according to the specified filter. */
  abstract read(filter?: {[name: string]: string}): Promise<DG.DataFrame>;

  /** Updates the entity in the database.
   * Returns the entity on success, or throws an exception otherwise. */
  abstract update(entity: TEntity): Promise<TEntity>;

  /** Deletes the entity in the database.
   * Returns the entity on success, or throws an exception otherwise. */
  abstract delete(entity: TEntity): Promise<TEntity>;
}


// currently works with columns from one table only
export class DbQueryEntityCrud<TEntity extends DbEntity = DbEntity> extends EntityCrud<TEntity> {
  connectionId: string;
  entityType: DbEntityType;

  constructor(connectionId: string, entityType: DbEntityType) {
    super();
    this.connectionId = connectionId;
    this.entityType = entityType;
  }

  sqlValue(x: any): string {
    if (x == null)
      return 'null';
    if (typeof x === 'string')
      return "'" + x + "'";
    return `${x}`;
  }

  sql(entity: TEntity, column: DbColumn): string {
    const value = entity.values[column.name]
    if (!value)
      return 'null';
    return column.type == 'string' ? "'" + value + "'" : `${value}`;
  }

  getWhereSql(filter: {[name: string]: string}): string {
    return Object
      .keys(filter)
      .map((columnName) => columnName + ' = ' + this.sqlValue(filter[columnName]))
      .join(' AND ');
  }

  getWhereKeySql(entity: TEntity): string {
    return `${entity.columns({isKey: true}).map(c => c.name + ' = ' + this.sql(entity, c)).join(' AND ')}`;
  }

  getCreateSql(entity: TEntity): string {
    return `insert into ${entity.entityType.table.name} 
(${entity.entityType.columns.map(c => c.name).join(', ')}) 
values (${entity.entityType.columns.map(c => this.sql(entity, c)).join(', ')})`;
  }

  getUpdateSql(entity: TEntity): string {
    return `update ${entity.entityType.table.name} ` +
      `set ${entity.columns({isKey: false}).map(c => c.name + ' = ' + this.sql(entity, c)).join(', ')} ` +
      `where ${this.getWhereKeySql(entity)}`;
  }

  getReadSql(filter?: {[name: string]: string}): string {
    return `select ${this.entityType.columns.map(c => c.name).join(', ')} ` +
           `from ${this.entityType.table.name}` +
            (filter ? ` where ${this.getWhereSql(filter)}` : '');
  }

  getDeleteSql(entity: TEntity) {
    return `delete from ${entity.entityType.table.name} where ${this.getWhereKeySql(entity)}`;
  }

  async create(entity: TEntity):  Promise<TEntity> {
    await grok.data.db.query(this.connectionId, this.getCreateSql(entity));
    return entity;
  }

  async read(filter?: {[name: string]: string}): Promise<DG.DataFrame> {
    return await grok.data.db.query(this.connectionId, this.getReadSql(filter));
  }

  async update(entity: TEntity): Promise<TEntity> {
    await grok.data.db.query(this.connectionId, this.getUpdateSql(entity));
    return entity;
  }

  async delete(entity: TEntity): Promise<TEntity> {
    await grok.data.db.query(this.connectionId, this.getDeleteSql(entity));
    return entity;
  }
}


export class CruddyViewConfig {
  entityType: DbEntityType;

  constructor(entityType: DbEntityType) {
    this.entityType = entityType;
  }
}


export class CruddyEntityView<TEntity extends DbEntityType = DbEntityType> extends DG.ViewBase {
  app: CruddyApp;
  entityType: DbEntityType;
  crud: DbQueryEntityCrud;
  grid?: DG.Grid;

  constructor(app: CruddyApp, entityType: DbEntityType) {
    super();
    this.app = app;
    this.entityType = entityType;
    this.name = entityType.type;
    this.crud = new DbQueryEntityCrud(app.config.connection, entityType);

    this.crud.read().then((df) => {
      this.grid = df.plot.grid();
      this.append(this.grid.root);
      this.grid.root.style.width = '100%';
      this.grid.root.style.height = '100%';
      this.initBehaviors();
    });
  }

  initBehaviors() {
    CruddyViewFeature.contextDetails().attach(this);
    CruddyViewFeature.editable().attach(this);
  }
}


export class CruddyViewFeature {
  name: string;
  attach: (view: CruddyEntityView) => void;

  constructor(name: string, attach: (view: CruddyEntityView) => void) {
    this.name = name;
    this.attach = attach;
  }

  static contextDetails(): CruddyViewFeature {
    return new CruddyViewFeature('context details', (v) => {
      v.grid!.onCurrentCellChanged.pipe(debounceTime(500)).subscribe((gridCell: DG.GridCell) => {
        let referenced = v.entityType.table.columns.filter((c) => c.referencedBy.length > 0);
        if (referenced.length == 0)
          return;

        let acc = ui.accordion(`entity-details-${this.name}`);
        for (let col of referenced)
          for (let c2 of col.referencedBy) {
            let host = ui.div();
            acc.addPane(c2.table.name, () => host);

            let detailsEntityType = v.app.config.getEntityType(c2.table);
            new DbQueryEntityCrud(v.app.config.connection, detailsEntityType)
              .read({[col.name]: gridCell.tableRow!.get(col.name)})
              .then((df) => host.appendChild(df.plot.grid().root));
          }
        acc.end();
        grok.shell.o = acc.root;
      });
    });
  }

  static editable(): CruddyViewFeature {
    return new CruddyViewFeature('editable', (v) => {
      v.grid!.onCellValueEdited.subscribe((gridCell: DG.GridCell) => {
        let crud = new DbQueryEntityCrud(v.app.config.connection, v.entityType);
        crud.update(v.entityType.rowToEntity(gridCell.tableRow!))
          .then((_) => grok.shell.info('Saved'));
      });
    });
  }
}


export class CruddyConfig {
  connection: string = '';
  tables: DbTable[] = [];
  entityTypes: DbEntityType[] = [];

  getTable(name: string): DbTable {
    return this.tables.find((t) => t.name == name)!;
  }

  getEntityType(table: DbTable): DbEntityType {
    return this.entityTypes.find((et) => et.table == table)!;
  }

  constructor(init?: Partial<CruddyConfig>) {
    Object.assign(this, init);
    for (let e of this.entityTypes)
      for (let c of e.table.columns)
        if (c.ref) {
          const [table, column] = c.ref.split('.');
          c.references = this.getTable(table).getColumn(column);
          c.references.referencedBy.push(c);
        }
  }
}


export class CruddyApp {
  config: CruddyConfig;
  views: CruddyEntityView[] = [];

  constructor(config: CruddyConfig, views?: CruddyEntityView[]) {
    this.config = config;
    this.views = views ?? this.views;
  }

  run(): void {
    // let v = new DG.MultiView({viewFactories: })

    for (let entityType of this.config.entityTypes) {
      let view = new CruddyEntityView(this, entityType);
      grok.shell.addView(view);
      this.views.push(view);
    }
  }
}