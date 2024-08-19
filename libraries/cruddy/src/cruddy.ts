import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';


export class DbTable {
  schema: string = '';
  name: string = '';

  constructor(schema: string, name: string) {
    this.schema = schema;
    this.name = name;
  }
}


/** Column belongs to a table. */
export class DbColumn {
  table: DbTable;
  name: string;
  type: string;
  isKey: boolean = false;
  prop?: DG.Property;  // additional properties such as semantic type, editor, formatting, etc.

  constructor(table: DbTable, name: string, type: string) {
    this.table = table;
    this.name = name;
    this.type = type;
  }
}


/** DbEntity is defined by a list of columns */
export class DbEntityType {
  type: string;
  columns: DbColumn[];

  constructor(type: string, columns: DbColumn[]) {
    this.type = type;
    this.columns = columns;
  }

  get table(): DbTable { return this.columns[0].table };
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
export class DbQueryEntityCrud<TEntity extends DbEntity> extends EntityCrud<TEntity> {
  connectionId: string;
  entityType: DbEntityType;

  constructor(connectionId: string, entityType: DbEntityType) {
    super();
    this.connectionId = connectionId;
    this.entityType = entityType;
  }

  getWhereKeySql(entity: TEntity): string {
    return `${entity.columns({isKey: true}).map(c => c.name == entity.values[c.name]).join(' AND ')}`;
  }

  getCreateSql(entity: TEntity): string {
    return `insert into ${entity.entityType.table.name} 
(${entity.entityType.columns.map(c => c.name).join(', ')}) 
values (${entity.entityType.columns.map(c => entity.values[c.name] ?? 'null')})`;
  }

  getUpdateSql(entity: TEntity): string {
    return `update ${entity.entityType.table.name} ` +
      `set ${entity.columns({isKey: false}).map(c => c.name + ' = ' + entity.values[c.name]).join(', ')} ` +
      `where ${this.getWhereKeySql(entity)}`;
  }

  getReadSql(filter?: {[name: string]: string}): string {
    return `select ${this.entityType.columns.join(', ')} from ${this.entityType.table}`;
  }

  getDeleteSql(entity: TEntity) {
    return `delete from ${entity.entityType.table.name} where `;
  }

  async create(entity: TEntity):  Promise<TEntity> {
    await grok.data.db.query(this.connectionId, this.getCreateSql(entity));
    return entity;
  }

  async read(filter?: {[name: string]: string}): Promise<DG.DataFrame> {
    return await grok.data.db.query(this.connectionId, this.getReadSql());
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


export abstract class EntityView<TEntity extends DbEntityType> extends DG.View {

}
