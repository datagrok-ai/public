/**
 * Rows of a packed macromolecule batch (chem-meta's `AsmBatch`, Bio-on-Crux PLAN §3.2) as a batch of their own: what
 * CruxService.assembleMacromolecules sends a worker per chunk. A part reads, row by row, as the rows read in the whole
 * batch (`unpackRow`): its symbol table holds only the symbols its rows use, its offsets start at 0, and its typed
 * arrays are its own (transferable, and the whole batch's are not detached by a transfer).
 */
import {ASM_CONN_STRIDE, AsmBatch} from '@datagrok-libraries/chem-meta/src/macromolecule-assembler';

/** The number of rows of `batch`. */
export function batchRowCount(batch: AsmBatch): number {
  return batch.rowStart.length - 1;
}

/** The number of monomers (gaps counted) of row `i`. */
export function rowMonomerCount(batch: AsmBatch, i: number): number {
  return batch.rowStart[i + 1] - batch.rowStart[i];
}

/** Rows `start` to `end` (excluded) of `batch`, as a batch. */
export function sliceBatch(batch: AsmBatch, start: number, end: number): AsmBatch {
  const n = batchRowCount(batch);
  if (!Number.isInteger(start) || !Number.isInteger(end) || start < 0 || end > n || start > end)
    throw new Error(`sliceBatch: rows ${start} to ${end} are not rows of a batch of ${n}`);
  return pickRuns(batch, [[start, end]]);
}

/**
 * Rows `rows` of `batch` (indices in increasing order, as CruxService's chunks of a batch hold them), as a batch: the
 * rows of a run of consecutive indices are copied at once.
 */
export function pickRows(batch: AsmBatch, rows: ArrayLike<number>): AsmBatch {
  const n = batchRowCount(batch);
  const runs: [number, number][] = [];
  for (let k = 0; k < rows.length; k++) {
    const row = rows[k];
    if (!Number.isInteger(row) || row < 0 || row >= n || (k > 0 && row <= rows[k - 1]))
      throw new Error(`pickRows: row ${row} is not a next row of a batch of ${n}`);
    if (runs.length > 0 && runs[runs.length - 1][1] === row)
      runs[runs.length - 1][1] = row + 1;
    else
      runs.push([row, row + 1]);
  }
  return pickRuns(batch, runs);
}

/** The rows of `runs` (`[start, end)` each, in order) as a batch. */
function pickRuns(batch: AsmBatch, runs: [number, number][]): AsmBatch {
  let rowCount = 0;
  let monCount = 0;
  let chainCount = 0;
  let connCount = 0;
  for (const [s, e] of runs) {
    rowCount += e - s;
    monCount += batch.rowStart[e] - batch.rowStart[s];
    chainCount += batch.rowChainStart[e] - batch.rowChainStart[s];
    connCount += batch.rowConnStart[e] - batch.rowConnStart[s];
  }
  const rowStart = new Int32Array(rowCount + 1);
  const monSymbol = new Int32Array(monCount);
  const monRole = new Uint8Array(monCount);
  const monPosition = new Int32Array(monCount);
  const rowChainStart = new Int32Array(rowCount + 1);
  const chainStart = new Int32Array(chainCount);
  const rowConnStart = new Int32Array(rowCount + 1);
  const conn = new Int32Array(connCount * ASM_CONN_STRIDE);
  const rowError: (string | null)[] = new Array(rowCount);
  let r = 0;
  let m = 0;
  let c = 0;
  let k = 0;
  for (const [s, e] of runs) {
    const m0 = batch.rowStart[s];
    const m1 = batch.rowStart[e];
    monSymbol.set(batch.monSymbol.subarray(m0, m1), m);
    monRole.set(batch.monRole.subarray(m0, m1), m);
    monPosition.set(batch.monPosition.subarray(m0, m1), m);
    // chain starts are row-local: copied as they are
    chainStart.set(batch.chainStart.subarray(batch.rowChainStart[s], batch.rowChainStart[e]), c);
    conn.set(batch.conn.subarray(batch.rowConnStart[s] * ASM_CONN_STRIDE, batch.rowConnStart[e] * ASM_CONN_STRIDE),
      k * ASM_CONN_STRIDE);
    for (let i = s; i < e; i++, r++) {
      rowStart[r + 1] = m + batch.rowStart[i + 1] - m0;
      rowChainStart[r + 1] = c + batch.rowChainStart[i + 1] - batch.rowChainStart[s];
      rowConnStart[r + 1] = k + batch.rowConnStart[i + 1] - batch.rowConnStart[s];
      rowError[r] = batch.rowError[i] ?? null; // a hole (a sparse array) as no error
    }
    m += m1 - m0;
    c += batch.rowChainStart[e] - batch.rowChainStart[s];
    k += batch.rowConnStart[e] - batch.rowConnStart[s];
  }
  // the symbols the rows use, in the order of their first use; -1 (a gap) stays
  const local = new Int32Array(batch.symbols.length).fill(-1);
  const symbols: AsmBatch['symbols'] = [];
  for (let i = 0; i < monSymbol.length; i++) {
    const si = monSymbol[i];
    if (si < 0)
      continue;
    if (si >= local.length)
      throw new Error(`a monomer of the batch names symbol ${si} of ${local.length}`);
    if (local[si] < 0) {
      local[si] = symbols.length;
      symbols.push(batch.symbols[si]);
    }
    monSymbol[i] = local[si];
  }
  return {symbols, rowStart, monSymbol, monRole, monPosition, chainStart, rowChainStart, conn, rowConnStart, rowError};
}
