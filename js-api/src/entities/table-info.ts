/**
 * TableInfo, ColumnInfo, and FileInfo classes.
 * @module entities/table-info
 */

import {toJs} from "../wrappers";
import {MapProxy} from "../proxies";
import {DataFrame} from "../dataframe";
import {IDartApi} from "../api/grok_api.g";
import {Tags} from "../api/ddt.api.g";
import {Entity} from "./entity";
import {DataConnection} from "./data-connection";
import dayjs from "dayjs";

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;


/** @extends Entity
 * Represents a Table metadata
 * */
export class TableInfo extends Entity {
  /** @constructs TableInfo */
  public tags: {[key: string]: any};

  constructor(dart: any) {
    super(dart);
    this.tags = new MapProxy(api.grok_TableInfo_Get_Tags(this.dart), 'tags');
  }

  static fromDataFrame(t: DataFrame): TableInfo {return toJs(api.grok_DataFrame_Get_TableInfo(t.dart)); }

  get dataFrame(): DataFrame { return toJs(api.grok_TableInfo_Get_DataFrame(this.dart)); }

  get columns(): ColumnInfo[] { return toJs(api.grok_TableInfo_Get_Columns(this.dart)); }
}


/** @extends Entity
 * Represents Column metadata */
export class ColumnInfo extends Entity {
  public tags: {[key: string]: any};
  /** @constructs ColumnInfo */
  constructor(dart: any) {
    super(dart);
    this.tags = new MapProxy(api.grok_ColumnInfo_Get_Tags(this.dart), 'tags');
  }

  get type(): string { return toJs(api.grok_ColumnInfo_Get_Type(this.dart)); }

  /** Returns reference information if the column is referencing another column in another table */
  get referenceInfo(): {table: string, column: string} | null {
    return this.tags[Tags.ReferencesTable] && this.tags[Tags.ReferencesColumn] ? {
      table: this.tags[Tags.ReferencesTable],
      column: this.tags[Tags.ReferencesColumn]
    } : null;
  }
}

/** @extends Entity
 * Allows for files handling in JS-based info panels
 * {@link https://datagrok.ai/help/discover/info-panels}
 * */
export class FileInfo extends Entity {
  /** @constructs FileInfo */
  constructor(dart: any) {
    super(dart);
  }

  get connection(): DataConnection { return api.grok_FileInfo_Get_Connection(toJs(this.dart)); }

  /** Returns path, i.e. `geo/dmv_offices.csv` */
  get path(): string { return api.grok_FileInfo_Get_Path(this.dart); }

  /** Returns full path, i.e. `Demo:TestJobs:Files:DemoFiles/geo/dmv_offices.csv` */
  get fullPath(): string { return api.grok_FileInfo_Get_FullPath(this.dart); }

  /** Returns URL path, i.e. `Demo.TestJobs.Files.DemoFiles/geo/dmv_offices.csv` */
  get viewPath(): string { return api.grok_FileInfo_Get_ViewPath(this.dart); }

  /** Returns file extension, i.e. `csv` */
  get extension(): string { return api.grok_FileInfo_Get_Extension(this.dart); }

  /** Returns file name, i.e. `dmv_offices.csv` */
  get fileName(): string { return api.grok_FileInfo_Get_FileName(this.dart); }

  /** Returns file URL */
  get url(): string { return api.grok_FileInfo_Get_Url(this.dart); }

  /** Checks if file */
  get isFile(): boolean { return api.grok_FileInfo_Get_IsFile(this.dart); }

  /** Checks if directory */
  get isDirectory(): boolean { return api.grok_FileInfo_Get_IsDirectory(this.dart); }

  get updatedOn(): dayjs.Dayjs | null {
    const d = api.grok_FileInfo_Get_UpdatedOn(this.dart);
    return d ? dayjs(d) : null;
  }

  /** @returns {Promise<string>} */
  // readAsString(): Promise<string> {
  //   return new Promise((resolve, reject) => api.grok_FileInfo_ReadAsString(this.dart, (x: any) => resolve(x), (x: any) => reject(x)));
  // }

  get data(): Uint8Array {return api.grok_FileInfo_Get_Data(this.dart);}

  readAsString(): Promise<string> {
    return api.grok_FileInfo_ReadAsString(this.dart);
  }

  /** @returns {Promise<Uint8Array>} */
  readAsBytes(): Promise<Uint8Array> {
    return api.grok_FileInfo_ReadAsBytes(this.dart);
  }

  static fromBytes(path: string, data: Uint8Array): FileInfo {
    if (!path)
      throw new Error('Path can\'t be null or empty');
    return api.grok_FileInfo_FromBytes(path, data);
  }

  static fromString(path: string, data: string): FileInfo {
    if (!path)
      throw new Error('Path can\'t be null or empty');
    return api.grok_FileInfo_FromString(path, data);
  }
}
