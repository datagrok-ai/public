import {getRdKitService} from '../../utils/chem-common-rdkit';
import {decomposeCluster, PositionRecord} from './sar-matrix-decompose';
import {CoreCluster, SarMatrix, SarMatrixCell, SarMatrixCellKind, SarMatrixColumn, SarMatrixRow}
  from './sar-matrix-types';

/** At most this many R-group positions become column groups (the two richest by default). */
const MAX_POSITIONS = 2;
/** At most this many substituent sub-columns are kept per position group. */
const MAX_COLUMNS_PER_POSITION = 8;
/** At most this many distinct cores become rows. */
const MAX_MATRIX_ROWS = 12;

/**
 * Fit a two-way additive (Free-Wilson) model over the observed cells of a matrix
 * and return a predictor for empty cells: `rowMean + columnMean - grandMean`.
 * A cell is only predictable when its row and its column each have at least one
 * observation.
 */
export function fitAdditiveModel(cells: SarMatrixCell[][], rowCount: number,
  columnCount: number): (rowIdx: number, colIdx: number) => {value: number, support: number} | null {
  const rowSum = new Float64Array(rowCount);
  const rowN = new Int32Array(rowCount);
  const colSum = new Float64Array(columnCount);
  const colN = new Int32Array(columnCount);
  let grandSum = 0;
  let grandN = 0;

  for (let ri = 0; ri < rowCount; ri++) {
    for (let ci = 0; ci < columnCount; ci++) {
      const cell = cells[ri][ci];
      if (cell.kind === 'real' && cell.value !== null) {
        rowSum[ri] += cell.value;
        rowN[ri]++;
        colSum[ci] += cell.value;
        colN[ci]++;
        grandSum += cell.value;
        grandN++;
      }
    }
  }

  const grandMean = grandN ? grandSum / grandN : 0;
  return (ri, ci) => (rowN[ri] && colN[ci]) ?
    {value: rowSum[ri] / rowN[ri] + colSum[ci] / colN[ci] - grandMean, support: Math.min(rowN[ri], colN[ci])} :
    null;
}

/** Most frequent value at a position across all records — the "held fixed" reference. */
function referenceValue(records: PositionRecord[], position: string): string {
  const counts = new Map<string, number>();
  for (const r of records) {
    const v = r.values[position];
    if (v)
      counts.set(v, (counts.get(v) ?? 0) + 1);
  }
  let best = '';
  let bestN = -1;
  for (const [v, n] of counts) {
    if (n > bestN) {
      best = v;
      bestN = n;
    }
  }
  return best;
}

/** Positions ranked by distinct-value count (richest first), capped at MAX_POSITIONS. */
function selectActivePositions(records: PositionRecord[], positions: string[]): string[] {
  return positions
    .map((p) => ({p, n: new Set(records.map((r) => r.values[p]).filter((v) => v)).size}))
    .filter((d) => d.n >= 2)
    .sort((a, b) => b.n - a.n)
    .slice(0, MAX_POSITIONS)
    .map((d) => d.p);
}

/** Rows: distinct core SMILES, most-populated first, capped at MAX_MATRIX_ROWS. */
function selectRows(records: PositionRecord[]): {coreSmiles: string, members: PositionRecord[]}[] {
  const byCore = new Map<string, PositionRecord[]>();
  for (const r of records) {
    let list = byCore.get(r.coreSmiles);
    if (!list)
      byCore.set(r.coreSmiles, list = []);
    list.push(r);
  }
  return [...byCore.entries()]
    .map(([coreSmiles, members]) => ({coreSmiles, members}))
    .sort((a, b) => b.members.length - a.members.length)
    .slice(0, MAX_MATRIX_ROWS);
}

/** Most frequent distinct values at a position, capped at MAX_COLUMNS_PER_POSITION. */
function topValues(records: PositionRecord[], position: string): string[] {
  const counts = new Map<string, number>();
  for (const r of records) {
    const v = r.values[position];
    if (v)
      counts.set(v, (counts.get(v) ?? 0) + 1);
  }
  return [...counts.entries()].sort((a, b) => b[1] - a[1])
    .slice(0, MAX_COLUMNS_PER_POSITION).map(([v]) => v);
}

/**
 * Step 6 + 7 (multi-position) — decompose every molecule in the cluster against one shared
 * anchor scaffold so R1, R2, ... are aligned across every row, then fill unmade
 * core × substituent combinations with Free-Wilson predictions. Falls back to
 * `assembleSinglePositionMatrix` when the cluster doesn't share a usable multi-position anchor.
 */
export async function assembleMultiPositionMatrix(cluster: CoreCluster, molecules: string[],
  activities: Float32Array, predict: boolean, scaffolds: string[] = []): Promise<SarMatrix> {
  const molIdx = [...new Set(cluster.series.flatMap((s) => s.members.map((m) => m.molIdx)))];
  const decomp = await decomposeCluster(molIdx, molecules, scaffolds);
  if (!decomp)
    return assembleSinglePositionMatrix(cluster, molecules, activities, predict);

  const activePositions = selectActivePositions(decomp.records, decomp.positions);
  if (activePositions.length === 0)
    return assembleSinglePositionMatrix(cluster, molecules, activities, predict);

  const refValues: {[position: string]: string} = {};
  for (const p of decomp.positions)
    refValues[p] = referenceValue(decomp.records, p);

  const rowGroups = selectRows(decomp.records);
  const rows: SarMatrixRow[] = rowGroups.map((g, i) => ({
    coreFragId: i, coreSmiles: g.coreSmiles, label: `Core ${i + 1}`,
  }));

  // Only records whose OTHER active positions sit at reference can ever populate a cell, so rank
  // each position's columns among those records — ranking over all records would spend the column
  // budget on substituents that leave the column empty.
  const eligibleFor = (position: string): PositionRecord[] => decomp.records.filter((r) =>
    activePositions.every((other) => other === position || r.values[other] === refValues[other]));

  const columns: SarMatrixColumn[] = [];
  const columnIndex = new Map<string, number>(); // `${position}|${value}` -> column index
  for (const position of activePositions) {
    for (const value of topValues(eligibleFor(position), position)) {
      columnIndex.set(`${position}|${value}`, columns.length);
      columns.push({position, substSmiles: value, count: 0});
    }
  }

  const cells: SarMatrixCell[][] = rows.map(() =>
    columns.map((): SarMatrixCell => ({kind: 'empty', value: null, molIdx: null, smiles: null})));

  // A record populates cell (row, position=value) only when every OTHER active position sits at
  // its reference value — the "vary one position, hold the rest at reference" slice.
  let minActivity = Infinity;
  let maxActivity = -Infinity;
  let realCount = 0;
  rowGroups.forEach((group, ri) => {
    for (const record of group.members) {
      const activity = activities[record.molIdx];
      if (!Number.isFinite(activity))
        continue; // missing activity — leave the combination unmade rather than fake a value
      for (const position of activePositions) {
        const heldAtRef = activePositions.every((other) =>
          other === position || record.values[other] === refValues[other]);
        if (!heldAtRef)
          continue;
        const ci = columnIndex.get(`${position}|${record.values[position]}`);
        if (ci === undefined)
          continue; // value fell outside the top-N columns kept for this position
        cells[ri][ci] = {kind: 'real', value: activity, molIdx: record.molIdx, smiles: molecules[record.molIdx]};
        columns[ci].count++;
        realCount++;
        if (activity < minActivity)
          minActivity = activity;
        if (activity > maxActivity)
          maxActivity = activity;
      }
    }
  });

  // Free-Wilson prediction runs independently per position group — each group is its own
  // rows × values 2-D slice.
  let virtualCount = 0;
  if (predict) {
    for (const position of activePositions) {
      const colIdxForPosition = columns
        .map((c, ci) => (c.position === position ? ci : -1)).filter((ci) => ci >= 0);
      const predictCell = fitAdditiveModel(
        cells.map((row) => colIdxForPosition.map((ci) => row[ci])),
        rows.length, colIdxForPosition.length);
      for (let ri = 0; ri < rows.length; ri++) {
        for (let k = 0; k < colIdxForPosition.length; k++) {
          const ci = colIdxForPosition[k];
          const predicted = predictCell(ri, k);
          if (predicted === null)
            continue;
          const existing = cells[ri][ci];
          if (existing.kind === 'empty') {
            cells[ri][ci] = {kind: 'virtual' as SarMatrixCellKind, value: predicted.value,
              molIdx: null, smiles: null, support: predicted.support};
            virtualCount++;
          } else if (existing.kind === 'real')
            existing.fit = predicted.value; // model's fitted value for an observed cell
        }
      }
    }
    await linkVirtualCellStructures(rows, columns, cells, refValues);
  }

  return {
    id: cluster.id, label: '', rows, columns, cells,
    minActivity: realCount ? minActivity : 0, maxActivity: realCount ? maxActivity : 0,
    realCount, virtualCount, scores: {}, positions: activePositions, refValues,
  };
}

/**
 * Assemble every virtual cell's SMILES by linking its row's multi-attachment-point core to one
 * fragment per active position — the varying position's column substituent, every other active
 * position's reference substituent — via the `linkRGroupFragments` worker call.
 */
async function linkVirtualCellStructures(rows: SarMatrixRow[], columns: SarMatrixColumn[],
  cells: SarMatrixCell[][], refValues: {[position: string]: string}): Promise<void> {
  // Every decomposed position must be filled, not just the active ones: the core carries an
  // attachment point per position, and any left open would leave a stray dummy atom.
  const allPositions = Object.keys(refValues);
  const attachIdx = allPositions.map((p) => Number.parseInt(p.replace(/^\D+/, ''), 10));
  if (attachIdx.some((n) => !Number.isFinite(n)))
    return;

  const cores: string[] = [];
  const fragmentColumns: string[][] = allPositions.map(() => []);
  const targets: SarMatrixCell[] = [];

  for (let ri = 0; ri < rows.length; ri++) {
    for (let ci = 0; ci < columns.length; ci++) {
      const cell = cells[ri][ci];
      if (cell.kind !== 'virtual' || cell.smiles !== null)
        continue;
      const varyingPosition = columns[ci].position;
      cores.push(rows[ri].coreSmiles);
      allPositions.forEach((position, k) => fragmentColumns[k].push(
        position === varyingPosition ? columns[ci].substSmiles : refValues[position]));
      targets.push(cell);
    }
  }
  if (targets.length === 0)
    return;
  const linked = await (await getRdKitService()).linkRGroupFragments(cores, fragmentColumns, attachIdx);
  for (let i = 0; i < targets.length; i++)
    targets[i].smiles = linked[i] || null;
}

/**
 * Step 6 + 7 (single-position fallback) — assemble one SAR Matrix from a single-cut cluster:
 * rows are the cluster's related cores, columns are the distinct varying substituents, and empty
 * cells are filled with Free-Wilson predictions. Used when multi-position decomposition can't
 * find a usable shared anchor.
 */
export function assembleSinglePositionMatrix(cluster: CoreCluster, molecules: string[], activities: Float32Array,
  predict: boolean): SarMatrix {
  const colIndex = new Map<string, number>();
  const columns: SarMatrixColumn[] = [];
  for (const series of cluster.series) {
    for (const member of series.members) {
      if (!colIndex.has(member.substSmiles)) {
        colIndex.set(member.substSmiles, columns.length);
        columns.push({position: 'R1', substSmiles: member.substSmiles, count: 0});
      }
    }
  }

  const rows: SarMatrixRow[] = cluster.series.map((series, i) => ({
    coreFragId: series.coreFragId,
    coreSmiles: series.coreSmiles,
    label: `Core ${i + 1}`,
  }));

  const cells: SarMatrixCell[][] = rows.map(() =>
    columns.map((): SarMatrixCell => ({kind: 'empty', value: null, molIdx: null, smiles: null})));

  let minActivity = Infinity;
  let maxActivity = -Infinity;
  let realCount = 0;
  cluster.series.forEach((series, ri) => {
    for (const member of series.members) {
      const activity = activities[member.molIdx];
      if (!Number.isFinite(activity))
        continue; // missing activity — leave the cell empty rather than fake a value
      const ci = colIndex.get(member.substSmiles)!;
      cells[ri][ci] = {kind: 'real', value: activity, molIdx: member.molIdx, smiles: molecules[member.molIdx]};
      columns[ci].count++;
      realCount++;
      if (activity < minActivity)
        minActivity = activity;
      if (activity > maxActivity)
        maxActivity = activity;
    }
  });

  let virtualCount = 0;
  if (predict) {
    const predictCell = fitAdditiveModel(cells, rows.length, columns.length);
    for (let ri = 0; ri < rows.length; ri++) {
      for (let ci = 0; ci < columns.length; ci++) {
        const predicted = predictCell(ri, ci);
        if (predicted === null)
          continue;
        const existing = cells[ri][ci];
        if (existing.kind === 'empty') {
          cells[ri][ci] = {kind: 'virtual' as SarMatrixCellKind, value: predicted.value,
            molIdx: null, smiles: null, support: predicted.support};
          virtualCount++;
        } else if (existing.kind === 'real')
          existing.fit = predicted.value; // model's fitted value for an observed cell
      }
    }
  }

  return {
    id: cluster.id,
    label: '',
    rows,
    columns,
    cells,
    minActivity: realCount ? minActivity : 0,
    maxActivity: realCount ? maxActivity : 0,
    realCount,
    virtualCount,
    scores: {},
    positions: ['R1'],
    refValues: {},
  };
}
