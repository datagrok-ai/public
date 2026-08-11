import {getRdKitService} from '../../utils/chem-common-rdkit';
import {decomposeCluster, PositionRecord} from './sar-matrix-decompose';
import {CoreCluster, SarMatrix, SarMatrixCell, SarMatrixCellKind, SarMatrixColumn, SarMatrixRow}
  from './sar-matrix-types';

/**
 * Fit a two-way additive (Free-Wilson) model over the observed cells of a matrix
 * and return a predictor for empty cells: `rowMean + columnMean - grandMean`.
 * A cell is only predictable when its row and its column each have at least one
 * observation.
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
    // `support` is the weaker of the two arms — a column mean drawn from one compound is shaky however
    // many the row has — and drives how faintly the cell is drawn. `references` is how many measured
    // compounds went into the estimate at all, which is what the panel lists and the filter counts.
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

/** Positions that actually vary, richest first. A position with a single value across the cluster is
 *  not an axis — every row would carry the same substituent there. */
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
 * Every candidate row: one per distinct (core, substituents at the folded positions). Folding every
 * position that is not the column axis into the row identity is what gives a compound varying at
 * several positions a cell of its own; keyed on the core alone, a cell would name no single compound
 * and the extra variants have to be discarded. Ranking and capping happen in the caller, which knows
 * the column set — a row can only be judged once it is known which of its members can land in a cell.
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

/** Distinct values at a position, most frequent first. Every value becomes a column: the grouping
 *  already decides which compounds belong together, and thin columns are surfaced by the pane's
 *  reference-point filter rather than dropped before the user can see them. */
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

  // One varying site per matrix: the richest position becomes the column axis and EVERY other
  // position folds into the row identity, not just the richest few. That is what lets a cell name
  // exactly one compound without holding anything at a reference value and discarding the rest — a
  // position left out of both axes is unconstrained, so two compounds differing only there would
  // claim the same cell and one would overwrite the other. Positions that never vary contribute a
  // single value and so cost no extra rows.
  const columnPosition = activePositions[0];
  const foldedPositions = decomp.positions.filter((p) => p !== columnPosition);

  // Still needed to complete a virtual structure at positions that vary in the cluster but are
  // neither the column axis nor part of the row identity.
  const refValues: {[position: string]: string} = {};
  for (const p of decomp.positions)
    refValues[p] = referenceValue(decomp.records, p);

  // The most frequent substituents are the candidate columns; which of them survive depends on the
  // rows that end up kept, so the column list is finalised below.
  const candidates = new Set(topValues(decomp.records, columnPosition));
  const observed = (record: PositionRecord): string | null => {
    const value = record.values[columnPosition] ?? '';
    return candidates.has(value) && Number.isFinite(activities[record.molIdx]) ? value : null;
  };

  // Rows still have to place something: one with no observed cell draws as a blank line, and with no
  // observation the additive model cannot put a prediction in it either. Ordered by how many members
  // can occupy a cell so the densest rows lead, but none are dropped for being far down that order —
  // which of them belong together is the grouping's decision, not a row budget's.
  const rowGroups = selectRows(decomp.records, foldedPositions)
    .map((group) => ({group, score: group.members.filter((r) => observed(r) !== null).length}))
    .filter((ranked) => ranked.score > 0)
    .sort((a, b) => b.score - a.score)
    .map((ranked) => ranked.group);

  // The mirror of that rule, on the other axis. Candidates are ranked over every record, so a
  // substituent whose compounds all sit in rows that placed nothing has nothing left to show — no
  // observation to draw and, with no observation in the column, nothing the additive model can
  // predict either. Such a column would render blank from top to bottom.
  const supported = new Set<string>();
  for (const group of rowGroups) {
    for (const record of group.members) {
      const value = observed(record);
      if (value !== null)
        supported.add(value);
    }
  }
  const columns: SarMatrixColumn[] = [];
  const columnIndex = new Map<string, number>(); // substituent value -> column index
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

  // Each record maps to exactly one cell: its row carries the core and the folded substituents, its
  // column the substituent at the varying position.
  let minActivity = Infinity;
  let maxActivity = -Infinity;
  const realMols = new Set<number>();
  rowGroups.forEach((group, ri) => {
    for (const record of group.members) {
      const activity = activities[record.molIdx];
      if (!Number.isFinite(activity))
        continue; // missing activity — leave the combination unmade rather than fake a value
      const ci = columnIndex.get(record.values[columnPosition] ?? '');
      if (ci === undefined)
        continue; // value fell outside the top-N columns kept for the varying position
      cells[ri][ci] = {kind: 'real', value: activity, molIdx: record.molIdx, smiles: molecules[record.molIdx]};
      realMols.add(record.molIdx);
      if (activity < minActivity)
        minActivity = activity;
      if (activity > maxActivity)
        maxActivity = activity;
    }
  });
  const realCount = realMols.size;

  // One fit over the whole grid: with a single varying position the matrix is already the two-way
  // table the additive model assumes, so every observation informs every row and column effect.
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
 * Fill in each row's `keySmiles` by attaching its folded substituents to the core and leaving the
 * column position's attachment point open. Rows keep their core untouched, because the virtual-cell
 * linker needs the fully open core to build a structure at every position.
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
    // The linker writes attachment points as isotopes; restoring the map-number form keeps the open
    // position drawn as an R label rather than as a numbered star.
    if (linked[i])
      row.keySmiles = linked[i].replace(/\[(\d+)\*\]/g, '[*:$1]');
  });
}

/**
 * Assemble every virtual cell's SMILES by linking its row's multi-attachment-point core to one
 * fragment per position: the column substituent at the varying position, the row's own substituent at
 * every position folded into the row, and the cluster reference elsewhere. Taking the folded
 * substituent from the row rather than from the reference is what keeps the proposed structure a
 * member of the row it is drawn in.
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
      const folded = rows[ri].foldedValues;
      cores.push(rows[ri].coreSmiles);
      // Membership of `folded` decides, not truthiness: a row that carries hydrogen at a folded
      // position stores an empty substituent, and treating that as "no value" would silently swap in
      // the cluster reference and build a structure belonging to a different row.
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
        columns.push({position: 'R1', substSmiles: member.substSmiles});
      }
    }
  }

  // Nothing is folded here — one position means the core alone identifies the row.
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
        continue; // missing activity — leave the cell empty rather than fake a value
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
