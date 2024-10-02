import * as DG from "datagrok-api/dg";

import {DbEntity} from "./entity";
import {Db} from "datagrok-api/dg";


/** Represents a database table. */
export class DbTable {
  schema?: DbSchema;
  name: string = '';
  columns: DbColumn[] = [];

  getColumn(name: string): DbColumn {
    return this.columns.find((c) => c.name == name)!;
  }

  // /** Finds a column that references the specified column in another table, written in the 'table.column' format */
  // getRefColumn(ref: string): DbColumn {
  //   return this.columns.find((c) => c.ref == ref)!;
  // }

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
  type: DG.Type = DG.TYPE.STRING;
  isKey: boolean = false;
  nullable: boolean = false;
  ref: string = '';
  references?: DbColumn;
  referencedBy: DbColumn[] = [];
  prop: DG.Property;  // works with the DbEntity

  /** Table-qualified name, such as "employee.first_name". */
  get tqName(): string { return `${this.table.name}.${this.name}`; }

  constructor(init?: Partial<DbColumn>) {
    Object.assign(this, init);

    const col = this;
    this.prop = DG.Property.create(this.name, this.type,
      (x: DbEntity) => x.values[this.name],
      function (x: DbEntity, v: any) {
        x.values[col.name] = v;
      });
  }
}


export class DbSchema {
  name: string;
  tables: DbTable[];

  constructor(name: string, tables: DbTable[]) {
    this.name = name;
    this.tables = tables;

    for (let table of this.tables)
      table.schema = this;

    for (const table of tables)
      for (const c of table.columns)
        if (c.ref) {
          const [table, column] = c.ref.split('.');
          c.references = this.getTable(table).getColumn(column);
          c.references.referencedBy.push(c);
        }
  }

  getTable(name: string): DbTable {
    return this.tables.find((c) => c.name == name)!;
  }

  /** [name] format: <tableName>.<columnName> */
  getColumn(name: string): DbColumn {
    const cols: DbColumn[] = this.tables.map((t) => t.columns).flat();
    return cols.find((c) => `${c.table.name}.${c.name}` === name)!;
  }
}