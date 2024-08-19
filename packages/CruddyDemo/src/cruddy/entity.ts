import {DbColumn, DbTable} from "./table";
import * as DG from "datagrok-api/dg";
import {IFilterDescription} from "./cruddy";
import {DbQueryEntityCrud} from "./crud";

/** DbEntity is defined by a list of columns */
export class DbEntityType {
  type: string = '';
  table: DbTable = new DbTable();
  behaviors: string[] = [];
  filters: IFilterDescription[] = [];
  crud: DbQueryEntityCrud = new DbQueryEntityCrud('', this);
  defaultView: string = 'grid';
  gridColumns?: DbColumn[];
  gridColumnsNames?: string[];

  constructor(init?: Omit<Partial<DbEntityType>, 'gridColumns'>) {
    Object.assign(this, init);
  }

  rowToEntity(row: DG.Row): DbEntity {
    return new DbEntity(this, Object.fromEntries(this.columns.map(c => [c.name, row.get(c.name)])));
  }

  get columns(): DbColumn[] {
    return this.table.columns;
  };

  getColumn(name: string): DbColumn {
    return this.columns.find(c => c.name === name)!;
  };

  get props(): DG.Property[] {
    return this.table.columns.map((c) => c.prop);
  };
}

export class DbEntity<T extends DbEntityType = DbEntityType> {
  entityType: T;
  values: { [key: string]: any };
  row?: DG.Row;

  constructor(entityType: T, values: { [p: string]: any }) {
    this.entityType = entityType;
    this.values = values;
  }

  /** Initializes an entity from the table row */
  static fromRow(entityType: DbEntityType, row: DG.Row) {
    const values = Object.fromEntries(entityType.table.columns.map((c) => [c.name, row.get(c.name)]));
    const entity = new DbEntity(entityType, values);
    entity.row = row;
    return entity;
  }

  /** Creates an empty instance of the entity, with all values assigned to null. */
  static createNew(entityType: DbEntityType): DbEntity {
    const values = Object.fromEntries(entityType.table.columns.map((c) => [c.name, null]));
    return new DbEntity(entityType, values);
  }

  columns(options?: { isKey: boolean }): DbColumn[] {
    return this.entityType.columns
      .filter(c => options?.isKey == undefined || options.isKey == c.isKey);
  }
}