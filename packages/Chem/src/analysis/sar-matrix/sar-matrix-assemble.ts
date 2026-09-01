import {_package} from '../../package';
import {getRdKitService} from '../../utils/chem-common-rdkit';
import {ClusterDecomposition, PositionRecord} from './sar-matrix-decompose';
import {CoreCluster, MatchedSeries, SarMatrix, SarMatrixCell, SarMatrixCellKind, SarMatrixColumn,
  SarMatrixRow}
  from './sar-matrix-types';

/**
 * Rows and columns joined into components by the cells actually observed, as a union-find over
 * `rowCount + columnCount` nodes: row `ri` is node `ri`, column `ci` is node `rowCount + ci`.
 */
function observedComponents(cells: SarMatrixCell[][], rowCount: number,
  columnCount: number): Int32Array {
  const parent = new Int32Array(rowCount + columnCount).map((_v, i) => i);
  const find = (x: number): number => {
    while (parent[x] !== x)
      x = parent[x] = parent[parent[x]];
    return x;
  };
  for (let ri = 0; ri < rowCount; ri++) {
    for (let ci = 0; ci < columnCount; ci++) {
      const cell = cells[ri][ci];
      if (cell.kind === 'real' && cell.value !== null)
        parent[find(ri)] = find(rowCount + ci);
    }
  }
  return parent.map((_v, i) => find(i));
}

/**
 * Fit a two-way additive (Free-Wilson) model and predict empty cells as
 * `rowMean + columnMean - grandMean`.
 *
 * A row and a column each having an observation is necessary but NOT sufficient: the model is only
 * identifiable within a connected component of the observed-cell graph. Two designs that share no
 * compound fix their row and column effects against unrelated baselines, so adding one's row effect to
 * the other's column effect is arithmetic over incomparable offsets — and it errs toward the grand
 * mean, which flatters exactly the core whose own measurements are worst.
 */
export function fitAdditiveModel(cells: SarMatrixCell[][], rowCount: number,
  columnCount: number): (rowIdx: number, colIdx: number) =>
    {value: number, support: number, references: number} | null {
  const rowAcc = new Float64Array(rowCount);
  const rowN = new Int32Array(rowCount);
  const colAcc = new Float64Array(columnCount);
  const colN = new Int32Array(columnCount);
  let grandSum = 0;
  // The fit sweeps the measured cells repeatedly and the grid is mostly holes, so they are collected
  // once here rather than rescanned.
  const obsRow: number[] = [];
  const obsCol: number[] = [];
  const obsVal: number[] = [];

  for (let ri = 0; ri < rowCount; ri++) {
    for (let ci = 0; ci < columnCount; ci++) {
      const cell = cells[ri][ci];
      if (cell.kind === 'real' && cell.value !== null) {
        rowN[ri]++;
        colN[ci]++;
        grandSum += cell.value;
        obsRow.push(ri);
        obsCol.push(ci);
        obsVal.push(cell.value);
      }
    }
  }

  const grandMean = obsVal.length ? grandSum / obsVal.length : 0;
  // Alternating least squares for value ≈ grand + row + column. Averaging the margins in one pass
  // solves that model only when every row is measured at the same columns; with holes it is a
  // different and worse estimator, and holes are the normal case here. Sweeping the two arms against
  // each other converges on the least-squares fit, and reproduces the one-pass answer exactly when
  // the design happens to be balanced.
  const rowEffect = new Float64Array(rowCount);
  const colEffect = new Float64Array(columnCount);
  for (let iteration = 0; iteration < MAX_FIT_SWEEPS; iteration++) {
    rowAcc.fill(0);
    for (let k = 0; k < obsVal.length; k++)
      rowAcc[obsRow[k]] += obsVal[k] - grandMean - colEffect[obsCol[k]];
    for (let ri = 0; ri < rowCount; ri++) {
      if (rowN[ri])
        rowEffect[ri] = rowAcc[ri] / rowN[ri];
    }
    colAcc.fill(0);
    for (let k = 0; k < obsVal.length; k++)
      colAcc[obsCol[k]] += obsVal[k] - grandMean - rowEffect[obsRow[k]];
    let shift = 0;
    for (let ci = 0; ci < columnCount; ci++) {
      if (!colN[ci])
        continue;
      const next = colAcc[ci] / colN[ci];
      shift = Math.max(shift, Math.abs(next - colEffect[ci]));
      colEffect[ci] = next;
    }
    if (shift < FIT_TOLERANCE)
      break;
  }

  const root = observedComponents(cells, rowCount, columnCount);
  return (ri, ci) => (rowN[ri] && colN[ci] && root[ri] === root[rowCount + ci]) ?
    // `support` is the weaker of the two arms (drives how faintly the cell draws); `references` is the
    // total measured compounds behind the estimate (what the panel lists and the filter counts).
    {value: grandMean + rowEffect[ri] + colEffect[ci],
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
    // Which of two equally rich positions becomes the column axis decides the whole matrix.
    .sort((a, b) => b.n - a.n || (a.p < b.p ? -1 : a.p > b.p ? 1 : 0))
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

/** Distinct values, most frequent first; the value itself breaks a tie so columns do not follow the
 *  order the records happened to arrive in. */
function rankByFrequency(values: Iterable<string>): string[] {
  const counts = new Map<string, number>();
  for (const value of values)
    counts.set(value, (counts.get(value) ?? 0) + 1);
  return [...counts.entries()]
    .sort((a, b) => b[1] - a[1] || (a[0] < b[0] ? -1 : a[0] > b[0] ? 1 : 0))
    .map(([value]) => value);
}

/** Distinct values at a position, most frequent first. Every value becomes a column; thin columns are
 *  surfaced by the pane's reference-point filter rather than dropped. */
function topValues(records: PositionRecord[], position: string): string[] {
  return rankByFrequency(records.map((r) => r.values[position]).filter((v) => v));
}

/** Ceiling on the cells one matrix may hold. Nothing upstream bounds a single-position matrix —
 *  `MAX_SAR_CLUSTER_SIZE` gates the decomposition, and a cluster past it arrives here with no anchor
 *  — so a runaway grouping would allocate millions of cells before anything could reject them. */
const MAX_MATRIX_CELLS = 250000;
/** Ceiling on rows; the cell ceiling then trims columns to fit around it. */
const MAX_MATRIX_ROWS = 2000;
/** Ceiling on columns, independent of the cell budget: the pane builds a grid column per substituent,
 *  so few rows must not buy unlimited width. */
const MAX_MATRIX_COLS = 500;
/** Sweep ceiling for the additive fit. It converges geometrically and reaches the tolerance in a
 *  handful of sweeps; the cap only bounds a pathologically ill-conditioned design. */
const MAX_FIT_SWEEPS = 50;
/** Largest column-effect move that still counts as converged, well under any activity unit in use. */
const FIT_TOLERANCE = 1e-9;

/** Claim a cell for a compound the set holds but has no activity for, so `fillVirtualCells` cannot
 *  predict it and offer it as one to make. Never over a measured cell: two members can share a slot,
 *  and a measurement outranks a blank. A compound that WAS assayed but whose value the chosen scale
 *  cannot represent is left alone — calling it untested would be a false claim about the data. */
function claimUnmeasured(cells: SarMatrixCell[][], ri: number, ci: number, molIdx: number,
  molecules: string[], assayed?: Uint8Array): void {
  if (assayed?.[molIdx])
    return;
  if (cells[ri][ci].kind === 'empty')
    cells[ri][ci] = {kind: 'unmeasured', value: null, molIdx, smiles: molecules[molIdx]};
}

/**
 * Cut the matrix down to the largest connected block of its observed cells, in place.
 *
 * The additive model reaches a cell only when its row and column sit in the same component, so
 * anything outside the largest one could never be filled and would draw as a permanent blank. A line
 * with no measurement at all is a singleton component and goes the same way.
 *
 * This DOES drop measured compounds, unlike a plain empty-line prune: a design sharing no compound
 * with the main block is not comparable to it, and showing both in one grid implies a comparison the
 * data cannot support. What survives is one grid where every cell is measured or predictable.
 */
function pruneUnobservedLines(rows: SarMatrixRow[], columns: SarMatrixColumn[],
  cells: SarMatrixCell[][]): void {
  const root = observedComponents(cells, rows.length, columns.length);
  const observedPerRoot = new Map<number, number>();
  for (let ri = 0; ri < rows.length; ri++) {
    for (let ci = 0; ci < columns.length; ci++) {
      if (cells[ri][ci].kind === 'real')
        observedPerRoot.set(root[ri], (observedPerRoot.get(root[ri]) ?? 0) + 1);
    }
  }
  let best = -1;
  let bestN = 0;
  // Ties break on the lower root so the surviving block does not depend on iteration order.
  for (const [r, n] of observedPerRoot) {
    if (n > bestN || (n === bestN && r < best)) {
      best = r;
      bestN = n;
    }
  }
  const keepRows = rows.map((_row, ri) => root[ri] === best);
  const keepCols = columns.map((_col, ci) => root[rows.length + ci] === best);
  if (keepRows.every(Boolean) && keepCols.every(Boolean))
    return;
  for (let ri = cells.length - 1; ri >= 0; ri--) {
    if (!keepRows[ri]) {
      cells.splice(ri, 1);
      rows.splice(ri, 1);
      continue;
    }
    for (let ci = keepCols.length - 1; ci >= 0; ci--) {
      if (!keepCols[ci])
        cells[ri].splice(ci, 1);
    }
  }
  for (let ci = keepCols.length - 1; ci >= 0; ci--) {
    if (!keepCols[ci])
      columns.splice(ci, 1);
  }
  // Labels number the rows as displayed, so they are reassigned rather than left with the gaps a
  // removed core would leave (Core 1, Core 3, Core 6).
  rows.forEach((row, ri) => row.label = `Core ${ri + 1}`);
}

/**
 * Fill every cell the additive model can reach, returning how many analogs were proposed.
 *
 * A real cell keeps its measured value: its fitted one comes from the leave-one-out pass in
 * `computeMatrixConfidence`, and an in-sample residual would let a cliff cell hide.
 */
function fillVirtualCells(cells: SarMatrixCell[][], rowCount: number, columnCount: number,
  predictUnmeasured: boolean): number {
  const predictCell = fitAdditiveModel(cells, rowCount, columnCount);
  let filled = 0;
  for (let ri = 0; ri < rowCount; ri++) {
    for (let ci = 0; ci < columnCount; ci++) {
      const predicted = predictCell(ri, ci);
      if (predicted === null)
        continue;
      const cell = cells[ri][ci];
      if (cell.kind === 'empty') {
        cells[ri][ci] = {kind: 'virtual' as SarMatrixCellKind, value: predicted.value,
          molIdx: null, smiles: null, support: predicted.support, references: predicted.references};
        filled++;
      } else if (predictUnmeasured && cell.kind === 'unmeasured') {
        // Keeps its kind, structure and row: the compound exists, so the prediction is an argument to
        // TEST it, and nothing downstream may offer it as one to synthesize.
        cell.value = predicted.value;
        cell.support = predicted.support;
        cell.references = predicted.references;
      }
    }
  }
  return filled;
}

/** The cluster's series cut to {@link MAX_MATRIX_ROWS}, densest first so a trim keeps the rows
 *  carrying the most data, broken by core so equal-sized series always cut the same way. */
function boundedSeries(cluster: CoreCluster): MatchedSeries[] {
  if (cluster.series.length <= MAX_MATRIX_ROWS)
    return cluster.series;
  const ranked = [...cluster.series].sort((a, b) => b.members.length - a.members.length ||
    (a.coreSmiles < b.coreSmiles ? -1 : a.coreSmiles > b.coreSmiles ? 1 : 0));
  _package.logger.warning(`SAR Matrix | a cluster of ${cluster.series.length} series was cut to ` +
    `${MAX_MATRIX_ROWS} rows; the rest carry the fewest compounds`);
  return ranked.slice(0, MAX_MATRIX_ROWS);
}

/**
 * Build a multi-position matrix from the cluster's precomputed decomposition (shared anchor so
 * R1, R2, ... align across rows), then fill unmade combinations with Free-Wilson predictions.
 * A null decomposition falls back to `assembleSinglePositionMatrix`.
 */
export async function assembleMultiPositionMatrix(cluster: CoreCluster, molecules: string[],
  activities: Float32Array, predict: boolean, decomp: ClusterDecomposition | null,
  predictUnmeasured = false, assayed?: Uint8Array): Promise<SarMatrix | null> {
  // Both ways into the fallback ask the same question, so the rule sits here rather than at the two
  // call sites: placeholder series have no real cores, and would read as one bare core per compound.
  const fallback = (): SarMatrix | null => cluster.requiresDecomposition ? null :
    assembleSinglePositionMatrix(cluster, molecules, activities, predict, predictUnmeasured, assayed);
  if (!decomp)
    return fallback();

  const activePositions = selectActivePositions(decomp.records, decomp.positions);
  if (activePositions.length === 0)
    return fallback();

  // Richest position is the column axis; EVERY other position folds into the row identity (not just
  // the richest few). A position left out of both axes is unconstrained, so two compounds differing
  // only there would collide in one cell.
  const columnPosition = activePositions[0];
  const foldedPositions = decomp.positions.filter((p) => p !== columnPosition);

  const refValues: {[position: string]: string} = {};
  for (const p of decomp.positions)
    refValues[p] = referenceValue(decomp.records, p);

  const columnValues = topValues(decomp.records, columnPosition);
  const candidates = new Set(columnValues);
  const observed = (record: PositionRecord): string | null => {
    const value = record.values[columnPosition] ?? '';
    return candidates.has(value) && Number.isFinite(activities[record.molIdx]) ? value : null;
  };

  // Keep only rows with at least one observed cell (a row with none would draw blank and the additive
  // model can predict nothing in it either), densest first.
  const rowGroups = selectRows(decomp.records, foldedPositions)
    .map((group) => ({group, score: group.members.filter((r) => observed(r) !== null).length}))
    .filter((ranked) => ranked.score > 0)
    .sort((a, b) => b.score - a.score || (a.group.coreSmiles < b.group.coreSmiles ? -1 :
      a.group.coreSmiles > b.group.coreSmiles ? 1 :
        JSON.stringify(a.group.folded) < JSON.stringify(b.group.folded) ? -1 : 1))
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
  for (const value of columnValues) {
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
      const ci = columnIndex.get(record.values[columnPosition] ?? '');
      if (ci === undefined)
        continue;
      if (!Number.isFinite(activity)) {
        claimUnmeasured(cells, ri, ci, record.molIdx, molecules, assayed);
        continue;
      }
      cells[ri][ci] = {kind: 'real', value: activity, molIdx: record.molIdx, smiles: molecules[record.molIdx]};
      realMols.add(record.molIdx);
      if (activity < minActivity)
        minActivity = activity;
      if (activity > maxActivity)
        maxActivity = activity;
    }
  });
  const realCount = realMols.size;
  pruneUnobservedLines(rows, columns, cells);

  let virtualCount = 0;
  if (predict) {
    virtualCount = fillVirtualCells(cells, rows.length, columns.length, predictUnmeasured);
    await linkVirtualCellStructures(rows, columns, cells, refValues);
  }

  return {
    id: cluster.id, label: cluster.label ?? '', rows, columns, cells,
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
  predict: boolean, predictUnmeasured = false, assayed?: Uint8Array): SarMatrix {
  const series = boundedSeries(cluster);

  const colIndex = new Map<string, number>();
  const columns: SarMatrixColumn[] = [];
  const substituentsOf = (): string[] => {
    // Unfiltered, unlike topValues: an empty substituent is a real column here, since a pooled series
    // carries placeholder members and MMP can mint an empty fragment name.
    const ranked = rankByFrequency(series.flatMap((m) => m.members.map((member) => member.substSmiles)));
    // Bounded on its own as well as against the row count: the cell budget alone would let a two-row
    // matrix have six figures of columns.
    const limit = Math.max(1, Math.min(MAX_MATRIX_COLS, Math.floor(MAX_MATRIX_CELLS / series.length)));
    if (ranked.length <= limit)
      return ranked;
    _package.logger.warning(`SAR Matrix | a matrix of ${ranked.length} substituents was cut to ${limit} columns; ` +
      'the rest occur least often');
    return ranked.slice(0, limit);
  };
  for (const substSmiles of substituentsOf()) {
    colIndex.set(substSmiles, columns.length);
    columns.push({position: 'R1', substSmiles});
  }

  const rows: SarMatrixRow[] = series.map((matched, i) => ({
    coreSmiles: matched.coreSmiles,
    keySmiles: matched.coreSmiles,
    label: `Core ${i + 1}`,
    foldedValues: {},
  }));

  const cells: SarMatrixCell[][] = rows.map(() =>
    columns.map((): SarMatrixCell => ({kind: 'empty', value: null, molIdx: null, smiles: null})));

  let minActivity = Infinity;
  let maxActivity = -Infinity;
  // Distinct compounds, not filled cells: one reached through two series sits in two cells, and
  // counting it twice overstates every "cpd" the UI shows.
  const realMols = new Set<number>();
  series.forEach((matched, ri) => {
    for (const member of matched.members) {
      const activity = activities[member.molIdx];
      const ci = colIndex.get(member.substSmiles);
      if (ci === undefined)
        continue;
      if (!Number.isFinite(activity)) {
        claimUnmeasured(cells, ri, ci, member.molIdx, molecules, assayed);
        continue;
      }
      cells[ri][ci] = {kind: 'real', value: activity, molIdx: member.molIdx, smiles: molecules[member.molIdx]};
      realMols.add(member.molIdx);
      if (activity < minActivity)
        minActivity = activity;
      if (activity > maxActivity)
        maxActivity = activity;
    }
  });
  const realCount = realMols.size;
  pruneUnobservedLines(rows, columns, cells);

  const virtualCount = predict ? fillVirtualCells(cells, rows.length, columns.length, predictUnmeasured) : 0;

  return {
    id: cluster.id,
    label: cluster.label ?? '',
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
