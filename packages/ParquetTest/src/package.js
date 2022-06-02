/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import { tableFromIPC, tableFromArrays, tableToIPC} from 'apache-arrow';
import { readParquet, writeParquet, WriterPropertiesBuilder, Compression } from './arrow1';

export const _package = new DG.Package();



//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//name: saveAsParquet
//description: Save as Parquet
//tags: fileExporter
export function saveAsParquet(){
  let table = grok.shell.t;
  let column_names = table.columns.names();
  //let arrays = table.columns.toList();
  const res = tableFromArrays({
    names: table.columns.byName(column_names[2]).toList(),
  });
  grok.shell.info(res);
  const arrowUint8Array = tableToIPC(res, "stream");
  const writerProperties = new WriterPropertiesBuilder().setCompression(Compression.SNAPPY).build();
  const parquetUint8Array = writeParquet(arrowUint8Array, writerProperties);
} 