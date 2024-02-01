import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {debounceTime} from 'rxjs/operators';
import {Subject} from "rxjs";

type TFilter = {[name: string]: string | null};

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


export interface IFilterDescription {
  type: 'distinct' | 'search' | 'range';
  column: string;
}


/** DbEntity is defined by a list of columns */
export class DbEntityType {
  type: string = '';
  table: DbTable = new DbTable();
  behaviors: string[] = [];
  filters: IFilterDescription[] = [];
  crud: DbQueryEntityCrud = new DbQueryEntityCrud('', this);

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


export abstract class EntityCrud<TEntity extends DbEntity = DbEntity> {

  /** Creates the entity in the database.
   * Returns the entity on success, or throws an exception otherwise. */
  abstract create(entity: TEntity): Promise<TEntity>;

  /** Reads the entities, according to the specified filter. */
  abstract read(filter?: TFilter): Promise<DG.DataFrame>;

  /** Updates the entity in the database.
   * Returns the entity on success, or throws an exception otherwise. */
  abstract update(entity: TEntity): Promise<TEntity>;

  /** Deletes the entity in the database.
   * Returns the entity on success, or throws an exception otherwise. */
  abstract delete(entity: TEntity): Promise<TEntity>;
}


export interface IQueryOptions {
  distinct?: boolean;
  columnNames?: string[];
  limit?: number;
  offset?: number;
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

  getWhereSql(filter: TFilter): string {
    return Object
      .keys(filter)
      .map((columnName) => filter[columnName] != null ? columnName + ' = ' + this.sqlValue(filter[columnName]) : columnName)
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

  getReadSql(filter?: TFilter, options?: IQueryOptions): string {
    const columnNames = options?.columnNames ?? this.entityType.columns.map(c => c.name);
    return `select ${options?.distinct ? 'distinct ' : ''} ${columnNames.join(', ')} ` +
           `from ${this.entityType.table.name}` +
            (filter && Object.keys(filter).length > 0 ? ` where ${this.getWhereSql(filter)}` : '') +
            (options?.limit ? ` limit ${options.limit}` : '') +
            (options?.offset ? ` limit ${options.offset}` : '');
  }

  getDeleteSql(entity: TEntity) {
    return `delete from ${entity.entityType.table.name} where ${this.getWhereKeySql(entity)}`;
  }

  async create(entity: TEntity):  Promise<TEntity> {
    await grok.data.db.query(this.connectionId, this.getCreateSql(entity));
    return entity;
  }

  async read(filter?: TFilter, options?: IQueryOptions): Promise<DG.DataFrame> {
    return await grok.data.db.query(this.connectionId, this.getReadSql(filter,  options));
  }

  async readSingle(keyFilter?: {[name: string]: string}): Promise<{[name: string]: any} | null> {
    const df = await grok.data.db.query(this.connectionId, this.getReadSql(keyFilter));
    if (df.rowCount == 0)
      return null;
    const resultRow: {[name: string]: any} = {};
    for (const column of df.columns)
      resultRow[column.name] = column.get(0);
    return resultRow;
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


export class CruddyFilterHost {
  root: HTMLDivElement = ui.divV([]);
  onChanged: Subject<TFilter> = new Subject<TFilter>();
  filters: CruddyFilter[] = [];

  constructor() {
    this.onChanged
      .pipe(debounceTime(100))
      .subscribe((_) => grok.shell.info('Filter changed: ' + JSON.stringify(this.getCondition())));
  }

  getCondition(): TFilter {
    return CruddyFilter.mergeConditions(this.filters.map((f) => f.getCondition()));
  }

  init(entityType: DbEntityType) {
    for (const f of entityType.filters) {
      const filter = CruddyFilter.create(entityType, f);
      filter.onChanged.subscribe((x) => this.onChanged.next(x));
      this.filters.push(filter);
      this.root.appendChild(filter.root);
    }
  }
}


export abstract class CruddyFilter {
  root: HTMLDivElement = ui.divV([]);
  entityType: DbEntityType;
  filter: IFilterDescription;
  onChanged: Subject<any> = new Subject<any>();

  abstract getCondition(): TFilter;

  static mergeConditions(conditions: TFilter[]): TFilter {
    const query = {};
    for (const f of conditions)
      Object.assign(query, f);
    return query;
  }

  protected constructor(entityType: DbEntityType, filter: IFilterDescription) {
    this.entityType = entityType;
    this.filter = filter;
  }

  static create(entityType: DbEntityType, filter: IFilterDescription): CruddyFilter {
    switch (filter.type) {
      case 'distinct': return new CruddyFilterCategorical(entityType, filter);
      case 'range': return new CruddyFilterRange(entityType, filter);
    }

    throw `Unknown filter type: ${filter.type}`;
  }
}


export class CruddyFilterCategorical extends CruddyFilter {
  choices?: DG.InputBase;

  constructor(entityType: DbEntityType, filter: IFilterDescription) {
    super(entityType, filter);
    entityType.crud
      .read({}, { distinct: true, columnNames: [filter.column]})
      .then((df) => {
        this.choices = ui.multiChoiceInput('values', [], df.columns.byIndex(0).toList());
        this.choices.onChanged(() => this.onChanged.next(this));
        this.root.appendChild(ui.h2(filter.column));
        this.root.appendChild(this.choices.input);
      });
  }

  getCondition(): TFilter {
    return {
      [`${this.filter.column} in (${(this.choices!.value as Array<any>).map((x) => `'${x}'`).join(', ')})`]: null
    };
  }
}


export class CruddyFilterRange extends CruddyFilter {
  slider = ui.rangeSlider(0, 1, 0, 1, false, 'thin_barbell');

  constructor(entityType: DbEntityType, filter: IFilterDescription) {
    super(entityType, filter);
    const minCol = `min(${filter.column}) as min`;
    const maxCol = `max(${filter.column}) as max`;
    entityType.crud
      .read({}, { columnNames: [minCol, maxCol]})
      .then((df) => {
        this.slider.setValues(df.get('min', 0), df.get('max', 0), df.get('min', 0), df.get('max', 0));
        this.slider.onValuesChanged.subscribe((_) => this.onChanged.next(this));
        this.root.appendChild(ui.h2(filter.column));
        this.root.appendChild(this.slider.root);
      });
  }

  getCondition(): TFilter {
    return {
      [`(${this.filter.column} >= ${this.slider.min} and ${this.filter.column} <= ${this.slider.max})`]: null
    };
  }
}


export class CruddyEntityView<TEntity extends DbEntityType = DbEntityType> extends DG.ViewBase {
  app: CruddyApp;
  entityType: DbEntityType;
  crud: DbQueryEntityCrud;
  grid?: DG.Grid;
  host: HTMLDivElement = ui.divH([]);
  filters: CruddyFilterHost = new CruddyFilterHost();

  constructor(app: CruddyApp, entityType: DbEntityType) {
    super();
    this.app = app;
    this.entityType = entityType;
    this.name = entityType.type;
    this.crud = new DbQueryEntityCrud(app.config.connection, entityType);

    this.crud.read().then((df) => {
      this.grid = df.plot.grid();
      this.append(this.host);
      this.host.appendChild(this.filters.root);
      this.host.appendChild(this.grid.root);
      this.grid.root.style.flexGrow = '1';
      //this.grid.root.style.height = '100%';
      this.initBehaviors();
      this.initFilters();
    });
  }

  initFilters() {
    this.filters.init(this.entityType);
    this.filters.onChanged.subscribe((q) => {
      this.crud.read(this.filters.getCondition()).then((df) => {
        this.grid!.dataFrame = df;
      });
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
    for (let e of this.entityTypes) {
      e.crud.connectionId = this.connection;
      for (let c of e.table.columns)
        if (c.ref) {
          const [table, column] = c.ref.split('.');
          c.references = this.getTable(table).getColumn(column);
          c.references.referencedBy.push(c);
        }
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
    for (let entityType of this.config.entityTypes) {
      let view = new CruddyEntityView(this, entityType);
      grok.shell.addView(view);
      this.views.push(view);
    }
  }
}