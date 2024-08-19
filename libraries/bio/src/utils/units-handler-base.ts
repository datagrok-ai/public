import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import wu from 'wu';
import {MoleculeUnits} from '../molecule/molecule-units-handler';
import {MoleculeBase, MoleculeBuildDataFunc} from '../molecule/types';

export abstract class DataBase {
  protected constructor(
    public readonly name?: string
  ) {}
}

export type BuildDataFunc<TValue, TData extends DataBase> = (src: TValue) => TData;


export abstract class UnitsHandlerBase<TValue, TData extends DataBase = DataBase> {
  public readonly units: string;

  protected constructor(
    protected readonly column: DG.Column<TValue>,
    protected readonly requiredSemType: string,
    protected readonly nameColNames: string[]
  ) {
    if (this.column.semType !== requiredSemType) {
      throw new Error(`Invalid semantic type '${this.column.semType} of column '${this.column.name}', ` +
        `expected '${requiredSemType}'.`);
    }
    this.units = this.column.meta.units as MoleculeUnits;
    if (!this.units)
      throw new Error(`Units not specified in column '${this.column.name}'.`);
  }

  protected _data: TData[] | undefined;

  public get data(): TData[] {
    if (this._data === undefined) {
      const buildDataFunc: BuildDataFunc<TValue, TData> = this.getBuildDataFunc();
      this._data = wu.count(0).take(this.column.length)
        .map((rowI) => {
          const v = this.column.get(rowI);
          return buildDataFunc(v!);
        }).toArray();
    }
    return this._data;
  }

  protected _nameCol: DG.Column<string> | null | undefined;

  public get nameCol(): DG.Column<string> | null {
    if (this._nameCol === undefined) {
      this._nameCol = wu(this.column.dataFrame.columns)
        .find((c) => {
          return c.type === DG.COLUMN_TYPE.STRING && this.nameColNames.includes(c.name.toLowerCase());
        }) as DG.Column<string>;
    }
    if (this._nameCol === undefined) {
      this._nameCol = wu(this.column.dataFrame.columns)
        .find((c) => {
          return c.type === DG.COLUMN_TYPE.STRING && this.nameColNames.some((cn) => c.name.includes(cn));
        });
    }
    if (this._nameCol === undefined)
      this._nameCol = null;
    return this._nameCol;
  }

  protected _fileExt: string | undefined;

  public fileExt(): string {
    if (this._fileExt === undefined)
      this._fileExt = this.getFileExt();
    return this._fileExt;
  }

  protected abstract getBuildDataFunc(): BuildDataFunc<TValue, TData>;

  protected abstract getFileExt(): string;

  public getFileName(cell: DG.Cell): string {
    let name: string | undefined = this.data[cell.rowIndex].name;
    const nameCol = this.nameCol;
    if (!name && nameCol != null)
      name = nameCol.get(cell.rowIndex) ?? undefined;
    if (!name)
      name = 'file';
    return `${name}.${this.getFileExt()}`;
  }
}
