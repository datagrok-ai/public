/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {default as init} from 'parquet-wasm';
import {
  fromFeather as _fromFeather,
  toFeather as _toFeather,
  fromParquet as _fromParquet,
  toParquet as _toParquet} from './api/api';

export * from './package.g';
export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

export class PackageFunctions {
  @grok.decorators.autostart({})
  static parquetInit() {
  }

  @grok.decorators.init({})
  static async initPackage() {
    await init(_package.webRoot + 'dist/parquet_wasm_bg.wasm');
  }

  @grok.decorators.func({description: 'Converts DG.DataFrame to arrow'})
  static toFeather(table: DG.DataFrame,
    @grok.decorators.param({options: {initialValue: 'true'}}) asStream: boolean = true): Uint8Array | null {
    return _toFeather(table, asStream);
  }

  @grok.decorators.func({description: 'Converts arrow ipc stream to DG.DataFrame'})
  static fromFeather(bytes: Uint8Array): DG.DataFrame | null {
    return _fromFeather(bytes);
  }

  @grok.decorators.func({description: 'Converts DG.DataFrame to parquet'})
  static toParquet(table: DG.DataFrame,
    @grok.decorators.param({type: 'int', options: {nullable: true}}) compression?: number): Uint8Array | null {
    return _toParquet(table, compression);
  }

  @grok.decorators.func({description: 'Converts binary data in parquet format to DG.DataFrame'})
  static fromParquet(bytes: Uint8Array): DG.DataFrame | null {
    return _fromParquet(bytes);
  }

  @grok.decorators.fileHandler({ext: 'parquet'})
  static parquetFileHandler(@grok.decorators.param({type: 'list'}) bytes: WithImplicitCoercion<ArrayBuffer | SharedArrayBuffer>): DG.DataFrame[] | null {
    const res = this.fromParquet(new Uint8Array(bytes as ArrayBuffer));
    return res ? [res] : null;
  }

  @grok.decorators.fileHandler({ext: 'feather'})
  static featherFileHandler(@grok.decorators.param({type: 'list'}) bytes: WithImplicitCoercion<ArrayBuffer | SharedArrayBuffer>): DG.DataFrame[] | null {
    const df = this.fromFeather(new Uint8Array(bytes as ArrayBuffer));
    return df ? [df] : null;
  }

  @grok.decorators.fileExporter({description: 'Save as Parquet'})
  static saveAsParquet() {
    const table = grok.shell.t;
    const parquetUint8Array = this.toParquet(table);
    DG.Utils.download(table.name + '.parquet', (parquetUint8Array ?? new Uint8Array(0)) as Uint8Array<ArrayBuffer>);
  }

  @grok.decorators.fileExporter({description: 'Save as Feather'})
  static saveAsFeather() {
    const table = grok.shell.t;
    const arrowUint8Array = this.toFeather(table, false);
    DG.Utils.download(table.name + '.feather', (arrowUint8Array ?? new Uint8Array(0)) as Uint8Array<ArrayBuffer>);
  }
}
