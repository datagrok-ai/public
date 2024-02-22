import * as DG from "datagrok-api/dg";
import * as grok from "datagrok-api/grok";
import {DbColumn} from "./table";
import {DbEntity, DbEntityType} from "./entity";
import {TFilter} from "./cruddy";

export abstract class Crud<TEntity extends DbEntity = DbEntity> {

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
export class DbQueryEntityCrud<TEntity extends DbEntity = DbEntity> extends Crud<TEntity> {
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
    return `select ${options?.distinct ? 'distinct ' : ''}${columnNames.join(', ')} ` +
      `from ${this.entityType.table.name}` +
      (filter && Object.keys(filter).length > 0 ? ` where ${this.getWhereSql(filter)}` : '') +
      (options?.limit ? ` limit ${options.limit}` : '') +
      (options?.offset ? ` limit ${options.offset}` : '');
  }

  getDeleteSql(entity: TEntity) {
    return `delete from ${entity.entityType.table.name} where ${this.getWhereKeySql(entity)}`;
  }

  async query(sql: string): Promise<DG.DataFrame> {
    console.log(sql);
    return await grok.data.db.query(this.connectionId, sql);
  }

  async create(entity: TEntity): Promise<TEntity> {
    await this.query(this.getCreateSql(entity));
    return entity;
  }

  async read(filter?: TFilter, options?: IQueryOptions): Promise<DG.DataFrame> {
    return await this.query(this.getReadSql(filter, options));
  }

  async readSingle(keyFilter?: { [name: string]: string }): Promise<{ [name: string]: any } | null> {
    const df = await this.query(this.getReadSql(keyFilter));
    if (df.rowCount == 0)
      return null;
    const resultRow: { [name: string]: any } = {};
    for (const column of df.columns)
      resultRow[column.name] = column.get(0);
    return resultRow;
  }

  async update(entity: TEntity): Promise<TEntity> {
    await this.query(this.getUpdateSql(entity));
    return entity;
  }

  async delete(entity: TEntity): Promise<TEntity> {
    await this.query(this.getDeleteSql(entity));
    return entity;
  }
}