/**
 * DataFrame indexing and validation utilities.
 */

import {TimeSeriesDataFrame} from './types';

/* ================================================================== */
/*  Index                                                              */
/* ================================================================== */

/** Contiguous row range for a single sample within the flat arrays. */
export interface SampleRange {
  /** First row index (inclusive). */
  start: number;
  /** Number of rows in this sample. */
  length: number;
}

/**
 * Single-pass O(n) index builder.
 *
 * Scans `df.ids` and returns a map from sample id to its row range.
 * Assumes ids are contiguous (all rows for one id appear together).
 */
export function buildIndex(df: TimeSeriesDataFrame): Map<number, SampleRange> {
  const map = new Map<number, SampleRange>();
  const {ids, rowCount} = df;
  if (rowCount === 0) return map;

  let prevId = ids[0];
  let start = 0;

  for (let i = 1; i < rowCount; i++) {
    if (ids[i] !== prevId) {
      map.set(prevId, {start, length: i - start});
      prevId = ids[i];
      start = i;
    }
  }
  map.set(prevId, {start, length: rowCount - start});
  return map;
}

/* ================================================================== */
/*  Validation                                                         */
/* ================================================================== */

/**
 * Validates a TimeSeriesDataFrame for structural correctness.
 *
 * Checks:
 * 1. `ids.length === time.length === rowCount`
 * 2. Every `column.data.length === rowCount`
 * 3. No NaN or ±Inf in Float32/Float64 columns
 * 4. Ids are contiguous (sorted by id)
 *
 * @throws Error with a descriptive message on failure.
 */
export function validate(df: TimeSeriesDataFrame): void {
  const {ids, time, rowCount, columns} = df;

  if (ids.length !== rowCount)
    throw new Error(`ids.length (${ids.length}) !== rowCount (${rowCount})`);
  if (time.length !== rowCount)
    throw new Error(`time.length (${time.length}) !== rowCount (${rowCount})`);

  for (const col of columns) {
    if (col.data.length !== rowCount)
      throw new Error(`Column "${col.name}" data.length (${col.data.length}) !== rowCount (${rowCount})`);

    if (col.data instanceof Float32Array || col.data instanceof Float64Array) {
      for (let i = 0; i < rowCount; i++) {
        const v = col.data[i];
        if (Number.isNaN(v))
          throw new Error(`Column "${col.name}" contains NaN at index ${i}`);
        if (!Number.isFinite(v))
          throw new Error(`Column "${col.name}" contains non-finite value at index ${i}`);
      }
    }
  }

  // Check ids are contiguous
  if (rowCount > 1) {
    const seen = new Set<number>();
    let currentId = ids[0];
    seen.add(currentId);
    for (let i = 1; i < rowCount; i++) {
      if (ids[i] !== currentId) {
        currentId = ids[i];
        if (seen.has(currentId))
          throw new Error(`ids are not contiguous: id ${currentId} reappears after a different id`);
        seen.add(currentId);
      }
    }
  }
}
