import {getRdKitService} from '../../utils/chem-common-rdkit';
import {ClusterDecomposition, PositionRecord} from './sar-matrix-decompose';
import {CoreCluster, SarMatrix, SarMatrixCell, SarMatrixCellKind, SarMatrixColumn, SarMatrixRow}
  from './sar-matrix-types';

/**
 * Fit a two-way additive (Free-Wilson) model and predict empty cells as
 * `rowMean + columnMean - grandMean`. A cell is predictable only when its row and
 * column each have at least one observation.
 */
export function fitAdditiveModel(cells: SarMatrixCell[][], rowCount: number,
  columnCount: number): (rowIdx: number, colIdx: number) =>
    {value: number, support: number, references: number} | null {
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
    // `support` is the weaker of the two arms (drives how faintly the cell draws); `references` is the
    // total measured compounds behind the estimate (what the panel lists and the filter counts).
    {value: rowSum[ri] / rowN[ri] + colSum[ci] / colN[ci] - grandMean,
      support: Math.min(rowN[ri], colN[ci]), references: rowN[ri] + colN[ci]} :
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

/** Positions that actually vary (>=2 distinct values), richest first. A single-valued position is not
 *  an axis. */
function selectActivePositions(records: PositionRecord[], positions: string[]): string[] {
  return positions
    .map((p) => ({p, n: new Set(records.map((r) => r.values[p]).filter((v) => v)).size}))
    .filter((d) => d.n >= 2)
    .sort((a, b) => b.n - a.n)
    .map((d) => d.p);
}

interface RowGroup {
  coreSmiles: string;
  folded: {[position: string]: string};
  members: PositionRecord[];
}

/**
 * One row per distinct (core, folded substituents). Folding every non-column-axis position into the
 * row identity is what gives each cell exactly one compound; keyed on the core alone, a cell would name
 * no single compound.
 */
function selectRows(records: PositionRecord[], foldedPositions: string[]): RowGroup[] {
  const byKey = new Map<string, RowGroup>();
  for (const r of records) {
    const folded: {[position: string]: string} = {};
    for (const position of foldedPositions)
      folded[position] = r.values[position] ?? '';
    const key = [r.coreSmiles, ...foldedPositions.map((p) => folded[p])].join('\0');
    let group = byKey.get(key);
    if (!group)
      byKey.set(key, group = {coreSmiles: r.coreSmiles, folded, members: []});
    group.members.push(r);
  }
  return [...byKey.values()];
}

/** Distinct values at a position, most frequent first. Every value becomes a column; thin columns are
 *  surfaced by the pane's reference-point filter rather than dropped. */
function topValues(records: PositionRecord[], position: string): string[] {
  const counts = new Map<string, number>();
  for (const r of records) {
    const v = r.values[position];
    if (v)
      counts.set(v, (counts.get(v) ?? 0) + 1);
  }
  return [...counts.entries()].sort((a, b) => b[1] - a[1]).map(([v]) => v);
}

/**
 * Build a multi-position matrix from the cluster's precomputed decomposition (shared anchor so
 * R1, R2, ... align across rows), then fill unmade combinations with Free-Wilson predictions.
 * A null decomposition falls back to `assembleSinglePositionMatrix`.
 */
export async function assembleMultiPositionMatrix(cluster: CoreCluster, molecules: string[],
  activities: Float32Array, predict: boolean, decomp: ClusterDecomposition | null): Promise<SarMatrix> {
  if (!decomp)
    return assembleSinglePositionMatrix(cluster, molecules, activities, predict);

  const activePositions = selectActivePositions(decomp.records, decomp.positions);
  if (activePositions.length === 0)
    return assembleSinglePositionMatrix(cluster, molecules, activities, predict);

  // Richest position is the column axis; EVERY other position folds into the row identity (not just
  // the richest few). A position left out of both axes is unconstrained, so two compounds differing
  // only there would collide in one cell.
  const columnPosition = activePositions[0];
  const foldedPositions = decomp.positions.filter((p) => p !== columnPosition);

  const refValues: {[position: string]: string} = {};
  for (const p of decomp.positions)
    refValues[p] = referenceValue(decomp.records, p);

  const candidates = new Set(topValues(decomp.records, columnPosition));
  const observed = (record: PositionRecord): string | null => {
    const value = record.values[columnPosition] ?? '';
    return candidates.has(value) && Number.isFinite(activities[record.molIdx]) ? value : null;
  };

  // Keep only rows with at least one observed cell (a row with none would draw blank and the additive
  // model can predict nothing in it either), densest first.
  const rowGroups = selectRows(decomp.records, foldedPositions)
    .map((group) => ({group, score: group.members.filter((r) => observed(r) !== null).length}))
    .filter((ranked) => ranked.score > 0)
    .sort((a, b) => b.score - a.score)
    .map((ranked) => ranked.group);

  // Mirror of the row rule on the column axis: drop columns whose compounds all sit in dropped rows,
  // as they would render blank top to bottom.
  const supported = new Set<string>();
  for (const group of rowGroups) {
    for (const record of group.members) {
      const value = observed(record);
      if (value !== null)
        supported.add(value);
    }
  }
  const columns: SarMatrixColumn[] = [];
  const columnIndex = new Map<string, number>();
  for (const value of topValues(decomp.records, columnPosition)) {
    if (!supported.has(value))
      continue;
    columnIndex.set(value, columns.length);
    columns.push({position: columnPosition, substSmiles: value});
  }

  const rows: SarMatrixRow[] = rowGroups.map((g, i) => ({
    coreSmiles: g.coreSmiles, keySmiles: g.coreSmiles, label: `Core ${i + 1}`,
    foldedValues: g.folded,
  }));
  await buildRowKeys(rows, foldedPositions);

  const cells: SarMatrixCell[][] = rows.map(() =>
    columns.map((): SarMatrixCell => ({kind: 'empty', value: null, molIdx: null, smiles: null})));

  let minActivity = Infinity;
  let maxActivity = -Infinity;
  const realMols = new Set<number>();
  rowGroups.forEach((group, ri) => {
    for (const record of group.members) {
      const activity = activities[record.molIdx];
      if (!Number.isFinite(activity))
        continue;
      const ci = columnIndex.get(record.values[columnPosition] ?? '');
      if (ci === undefined)
        continue;
      cells[ri][ci] = {kind: 'real', value: activity, molIdx: record.molIdx, smiles: molecules[record.molIdx]};
      realMols.add(record.molIdx);
      if (activity < minActivity)
        minActivity = activity;
      if (activity > maxActivity)
        maxActivity = activity;
    }
  });
  const realCount = realMols.size;

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
            molIdx: null, smiles: null, support: predicted.support, references: predicted.references};
          virtualCount++;
        }
        // Real cells get their fitted value from the leave-one-out pass in computeMatrixConfidence,
        // not from this in-sample fit — an in-sample residual would let a cliff cell hide.
      }
    }
    await linkVirtualCellStructures(rows, columns, cells, refValues);
  }

  return {
    id: cluster.id, label: '', rows, columns, cells,
    minActivity: realCount ? minActivity : 0, maxActivity: realCount ? maxActivity : 0,
    realCount, virtualCount, scores: {}, positions: [columnPosition], refValues,
    siteKey: cluster.siteKey,
    level: cluster.level,
    parentId: cluster.parentId,
  };
}

/**
 * Fill each row's `keySmiles` by attaching its folded substituents to the core, leaving the column
 * position open. `coreSmiles` is left untouched — the virtual-cell linker needs the fully open core.
 */
async function buildRowKeys(rows: SarMatrixRow[], foldedPositions: string[]): Promise<void> {
  if (foldedPositions.length === 0 || rows.length === 0)
    return;
  const attachIdx = foldedPositions.map((p) => Number.parseInt(p.replace(/^\D+/, ''), 10));
  if (attachIdx.some((n) => !Number.isFinite(n)))
    return;
  const cores = rows.map((row) => row.coreSmiles);
  const fragmentColumns = foldedPositions.map((p) => rows.map((row) => row.foldedValues[p] ?? ''));
  const linked = await (await getRdKitService()).linkRGroupFragments(cores, fragmentColumns, attachIdx);
  rows.forEach((row, i) => {
    if (linked[i])
      row.keySmiles = linked[i];
  });
}

/**
 * Assemble each virtual cell's SMILES by linking the row's core to one fragment per position: the
 * column substituent at the varying position, the row's own substituent at folded positions, the
 * cluster reference elsewhere. The folded substituent comes from the row (not the reference) so the
 * proposed structure stays a member of its row.
 */
async function linkVirtualCellStructures(rows: SarMatrixRow[], columns: SarMatrixColumn[],
  cells: SarMatrixCell[][], refValues: {[position: string]: string}): Promise<void> {
  // Fill every decomposed position, not just active ones: any attachment point left open would leave
  // a stray dummy atom.
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
      const folded = rows[ri].foldedValues;
      cores.push(rows[ri].coreSmiles);
      // Test membership (`in`), not truthiness: a folded hydrogen is stored as an empty string, and
      // treating that as "no value" would swap in the cluster reference and build a different row's
      // structure.
      allPositions.forEach((position, k) => fragmentColumns[k].push(
        position === varyingPosition ? columns[ci].substSmiles :
          (position in folded ? folded[position] : refValues[position])));
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
 * Single-position fallback — rows are the cluster's cores, columns the distinct varying substituents,
 * empty cells filled with Free-Wilson predictions. Used when no usable shared anchor is found.
 */
export function assembleSinglePositionMatrix(cluster: CoreCluster, molecules: string[], activities: Float32Array,
  predict: boolean): SarMatrix {
  const colIndex = new Map<string, number>();
  const columns: SarMatrixColumn[] = [];
  for (const series of cluster.series) {
    for (const member of series.members) {
      if (!colIndex.has(member.substSmiles)) {
        colIndex.set(member.substSmiles, columns.length);
        columns.push({position: 'R1', substSmiles: member.substSmiles});
      }
    }
  }

  const rows: SarMatrixRow[] = cluster.series.map((series, i) => ({
    coreSmiles: series.coreSmiles,
    keySmiles: series.coreSmiles,
    label: `Core ${i + 1}`,
    foldedValues: {},
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
        continue;
      const ci = colIndex.get(member.substSmiles)!;
      cells[ri][ci] = {kind: 'real', value: activity, molIdx: member.molIdx, smiles: molecules[member.molIdx]};
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
            molIdx: null, smiles: null, support: predicted.support, references: predicted.references};
          virtualCount++;
        }
        // Real cells' fitted value comes from the leave-one-out pass, not this in-sample fit.
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
    siteKey: cluster.siteKey,
    level: cluster.level,
    parentId: cluster.parentId,
  };
}
