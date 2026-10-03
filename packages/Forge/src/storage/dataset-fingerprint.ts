import * as DG from 'datagrok-api/dg';

export interface ColumnFingerprint {
  name: string;
  type: string;
  missingCount: number;
  min?: number;
  max?: number;
  mean?: number;
  categories?: string[];
}

export interface DatasetFingerprint {
  rowCount: number;
  columnCount: number;
  hash: string;
  columns: ColumnFingerprint[];
}

const MAX_CATEGORIES = 20;
const FNV_OFFSET = 0x811c9dc5;
const FNV_PRIME = 0x01000193;

export function datasetFingerprint(features: DG.DataFrame, target: DG.Column): DatasetFingerprint {
  const columns = [...features.columns.toList(), target];
  const encoder = new TextEncoder();
  let hash = FNV_OFFSET;
  for (const col of columns) {
    hash = fnv1a(hash, encoder.encode(`${col.name}\u0000${col.type}\u0000`));
    // Bigint columns have no raw buffer.
    if (col.matches('numerical') && col.type !== DG.COLUMN_TYPE.BIG_INT)
      hash = fnv1a(hash, rawBytes(col));
    else {
      for (let i = 0; i < col.length; i++)
        hash = fnv1a(hash, encoder.encode(`${col.isNone(i) ? '' : String(col.get(i))}\u0000`));
    }
  }
  return {
    rowCount: target.length,
    columnCount: columns.length,
    hash: hash.toString(16).padStart(8, '0'),
    columns: columns.map(columnFingerprint),
  };
}

// The raw buffer may be longer than the column.
function rawBytes(col: DG.Column): Uint8Array {
  const data = col.getRawData();
  return new Uint8Array(data.buffer, data.byteOffset, col.length * data.BYTES_PER_ELEMENT);
}

function fnv1a(hash: number, bytes: Uint8Array): number {
  let h = hash;
  for (let i = 0; i < bytes.length; i++)
    h = Math.imul(h ^ bytes[i], FNV_PRIME) >>> 0;
  return h;
}

function columnFingerprint(col: DG.Column): ColumnFingerprint {
  const stats = col.stats;
  const fingerprint: ColumnFingerprint = {name: col.name, type: col.type, missingCount: stats.missingValueCount};
  if (col.matches('numerical'))
    return stats.valueCount > 0 ? {...fingerprint, min: stats.min, max: stats.max, mean: stats.avg} : fingerprint;
  const isCategorical = col.type === DG.COLUMN_TYPE.STRING || col.type === DG.COLUMN_TYPE.BOOL;
  return isCategorical && col.categories.length <= MAX_CATEGORIES ? {...fingerprint, categories: col.categories} :
    fingerprint;
}
