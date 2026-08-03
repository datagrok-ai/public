// Feather (Arrow IPC) ↔ DG.DataFrame conversion has moved to
// @datagrok-libraries/arrow. This module re-exports it so the package's
// public surface (`Arrow:toFeather`, `Arrow:fromFeather`) is unchanged, and
// keeps the parquet pair in place since parquet-wasm's WASM init lives with
// the package's webRoot.

import * as DG from 'datagrok-api/dg';
import {Compression, readParquet, Table, writeParquet, WriterPropertiesBuilder} from 'parquet-wasm';
import {toFeather, fromFeather} from '@datagrok-libraries/arrow';

export {toFeather, fromFeather};

export function toParquet(table: DG.DataFrame, compression?: Compression): Uint8Array | null {
  const arrowUint8Array = toFeather(table);
  if (arrowUint8Array == null) return null;
  const writerProperties = new WriterPropertiesBuilder().setCompression(compression ?? Compression.SNAPPY).build();
  return writeParquet(Table.fromIPCStream(arrowUint8Array), writerProperties);
}

export function fromParquet(bytes: Uint8Array): DG.DataFrame | null {
  if (!bytes) return null;
  const arrowUint8Array = readParquet(bytes).intoIPCStream();
  return fromFeather(arrowUint8Array);
}
