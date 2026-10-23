// Pure selection/view-model helpers for the RunComparison UI — index/split column pickers,
// target filtering, multi-value compatibility, per-entry statuses. No DG/platform dependencies.
// User documentation lives in help/compute/run-comparison.md

import {
  ColumnInfo, ComparisonEntryNodes, ComparisonTarget, ColumnTarget, isNumericType,
  TimeUnit, TIME_UNIT_MS, TIME_UNITS, AxisMode, AxisModeConfig, AxisModeSelection, AxisConfigMap, toTargetUnits,
  isTimeIndexType, isIndexCandidateType, TableNodeInfo,
} from './types';
import {
  normalizeName, nameSimilarity, nameMatchConfidence, normalizeUnits, unitsCompatibility, UnitConverter,
  FUZZY_NAME_THRESHOLD,
} from './matching';

/** The name, or with the first free ` (n)` suffix when another entry has it: charts split by run name. */
export function uniqueEntryName(name: string, taken: string[]): string {
  const names = new Set(taken);
  if (!names.has(name))
    return name;
  let n = 2;
  while (names.has(`${name} (${n})`))
    n++;
  return `${name} (${n})`;
}

/** Entry nodes without the outputs of steps whose results are outdated. */
export const withoutOutdatedOutputs = (nodes: ComparisonEntryNodes): ComparisonEntryNodes => ({
  ...nodes,
  scalars: nodes.scalars.filter((scalar) => !scalar.outdated),
  tables: nodes.tables.filter((table) => !table.outdated),
});

/** Substring match on the displayed text, with fuzzy per-token fallback for typos. */
export function matchesFilter(query: string, text: string): boolean {
  const q = normalizeName(query);
  if (!q)
    return true;
  const t = normalizeName(text);
  return t.includes(q) || t.split(' ').some((token) => nameSimilarity(token, q) >= FUZZY_NAME_THRESHOLD);
}

/**
 * Suggestion predicate for multi-value mode: targets of the anchor's kind, never mixed.
 * Every scalar is compatible (all scalars are numeric by construction); column targets
 * when both the anchor's and their own index are line-chartable (numeric or datetime)
 * everywhere — picks from different tables of a run chart as separate steps. Validity of
 * an already selected combination is enforced by builder padding, not here — selected
 * targets are never ejected.
 */
export function compatibleTargetsFor(
  anchor: ComparisonTarget | null,
  targets: ComparisonTarget[],
  getColumnType: (entryId: string, tablePath: string, columnName: string) => string | undefined,
): ComparisonTarget[] {
  if (!anchor)
    return [];
  if (anchor.kind === 'scalar')
    return targets.filter((target) => target.kind === 'scalar');
  const lineIndexed = (target: ColumnTarget) => target.bindings.every((b) =>
    isTimeIndexType(getColumnType(b.entryId, b.tablePath, b.indexColumnName)));
  if (!lineIndexed(anchor))
    return [];
  return targets.filter((target) => target.kind === 'column' && lineIndexed(target));
}

/**
 * Per selected target, the charted runs it has no data for: runs any selected target charts
 * (union), minus its own picks.
 */
export function selectionGaps(selected: ComparisonTarget[]): Map<string, string[]> {
  const charted = new Set(selected.flatMap((target) =>
    (target.bindings as {entryId: string}[]).map((b) => b.entryId)));
  return new Map(selected.map((target) => {
    const own = new Set((target.bindings as {entryId: string}[]).map((b) => b.entryId));
    return [target.key, [...charted].filter((entryId) => !own.has(entryId))];
  }));
}

export type ColumnValuesGetter =
  (entryId: string, tablePath: string, columnName: string) => ArrayLike<unknown> | undefined;

export const EQUALITY_REL_TOL = 1e-3;

// PEP 485 (math.isclose) with zero absolute tolerance: relative closeness only,
// so near-zero values are never conflated with zero at an arbitrary scale
const isClose = (a: number, b: number) =>
  Number.isFinite(a) && Number.isFinite(b) &&
  Math.abs(a - b) <= EQUALITY_REL_TOL * Math.max(Math.abs(a), Math.abs(b));

const sameValue = (x: unknown, y: unknown) => {
  if (x === y || (x !== x && y !== y))
    return true;
  return typeof x === 'number' && typeof y === 'number' && isClose(x, y);
};

const exactValue = (x: unknown, y: unknown) => x === y || (x !== x && y !== y);

const toTargetValues = (values: ArrayLike<unknown> | undefined, factor?: number) =>
  values == null || factor == null ? values : Array.from(values, (value) => toTargetUnits(value, factor));

/**
 * True when every bound run carries the same data for the target: scalar values or
 * element-wise value/index/split columns, with numbers compared to within
 * EQUALITY_REL_TOL relative tolerance. Datetime index columns compare exactly:
 * their values are epoch ms, where relative tolerance would swallow real shifts
 * (0.1% of a present-day timestamp is weeks). Unfetchable column data counts as differing.
 */
export function isTargetEqualAcrossRuns(
  target: ComparisonTarget,
  getColumnValues: ColumnValuesGetter,
  getColumnType?: (entryId: string, tablePath: string, columnName: string) => string | undefined,
): boolean {
  if (target.bindings.length < 2)
    return false;
  if (target.kind === 'scalar') {
    const first = toTargetUnits(target.bindings[0].value, target.bindings[0].factor);
    return target.bindings.every((b) => sameValue(toTargetUnits(b.value, b.factor), first));
  }
  const indexCompare = getColumnType != null && target.bindings.every((b) =>
    getColumnType(b.entryId, b.tablePath, b.indexColumnName) === 'datetime') ? exactValue : sameValue;
  const sameValues = (
    a: ArrayLike<unknown> | null | undefined, b: ArrayLike<unknown> | null | undefined,
    compare = sameValue,
  ) => {
    if (a == null || b == null || a.length !== b.length)
      return false;
    for (let i = 0; i < a.length; i++) {
      if (!compare(a[i], b[i]))
        return false;
    }
    return true;
  };
  const series = target.bindings.map((b) => ({
    value: toTargetValues(getColumnValues(b.entryId, b.tablePath, b.columnName), b.factor),
    index: getColumnValues(b.entryId, b.tablePath, b.indexColumnName),
    split: b.splitColumnName ? getColumnValues(b.entryId, b.tablePath, b.splitColumnName) : null,
  }));
  const [first, ...rest] = series;
  if (!first.value || !first.index || first.split === undefined)
    return false;
  return rest.every((s) => sameValues(s.value, first.value) &&
    sameValues(s.index, first.index, indexCompare) &&
    (s.split === first.split || sameValues(s.split, first.split)));
}

// bare 'm' is deliberately not recognized — ambiguous with meters
export function parseTimeUnit(units?: string): TimeUnit | undefined {
  const unit = normalizeUnits(units);
  return TIME_UNITS.find((item) => item === unit);
}

/**
 * Non-default view of the stored axis-mode picks: tables whose current index column
 * is numeric or datetime. Stale picks (index switched away) don't resolve but are kept.
 */
export function resolveAxisModes(
  entries: ComparisonEntryNodes[],
  indexColumns: Map<string, Map<string, string>>,
  axisModeSelection: AxisModeSelection,
): AxisConfigMap {
  const map: AxisConfigMap = new Map();
  for (const entry of entries) {
    for (const table of entry.tables) {
      const stored = axisModeSelection[entry.entryId]?.[table.path];
      if (stored == null || stored.mode === 'series')
        continue;
      const indexName = indexColumns.get(entry.entryId)?.get(table.path);
      const column = indexName ? table.columns.find((col) => col.name === indexName) : undefined;
      if (!column || !isTimeIndexType(column.type))
        continue;
      const tables = map.get(entry.entryId) ?? new Map<string, {mode: 'timeseries' | 'points', units?: TimeUnit}>();
      tables.set(table.path, stored.mode === 'timeseries' && isNumericType(column.type) ?
        {mode: stored.mode, units: stored.units ?? parseTimeUnit(column.units) ?? 's'} : {mode: stored.mode});
      map.set(entry.entryId, tables);
    }
  }
  return map;
}

/** Require-all gate: bindings chart in the given mode only when every bound table is configured. */
export function targetAxisMode(
  target: {bindings: {entryId: string, tablePath: string}[]},
  axisModes: AxisConfigMap,
  mode: 'timeseries' | 'points',
): 'full' | 'partial' | 'none' {
  if (target.bindings.length === 0)
    return 'none';
  const configured = target.bindings
    .filter((b) => axisModes.get(b.entryId)?.get(b.tablePath)?.mode === mode).length;
  if (configured === 0)
    return 'none';
  return configured === target.bindings.length ? 'full' : 'partial';
}

// display unit of the elapsed axis: the first numeric binding's units in participation
// order; datetime-only targets have no units selector anywhere, so the unit is auto-picked
export function resolveDisplayUnit(
  bindings: {entryId: string, tablePath: string}[],
  axisModes: AxisConfigMap,
): TimeUnit | 'auto' {
  for (const binding of bindings) {
    const config = axisModes.get(binding.entryId)?.get(binding.tablePath);
    if (config?.mode === 'timeseries' && config.units)
      return config.units;
  }
  return 'auto';
}

/** Largest of days/h/min/s spanning at least two whole units; sub-2s ranges fall to ms. */
export function pickAutoUnit(maxElapsedMs: number): TimeUnit {
  for (const unit of ['days', 'h', 'min', 's'] as TimeUnit[]) {
    if (maxElapsedMs >= 2 * TIME_UNIT_MS[unit])
      return unit;
  }
  return 'ms';
}

export const isSplitCandidate = (column: ColumnInfo, indexColumnName: string) =>
  column.type === 'string' && column.name !== indexColumnName;

export interface PickerSelections {
  index: Record<string, Record<string, string>>;
  split: Record<string, Record<string, string>>;
  axis: AxisModeSelection;
}

// tables merged into one picker row: the same output of the same function
const sameFuncKey = (table: TableNodeInfo) => table.nqName ? `${table.nqName}|${table.path.split('/').pop()}` : null;

export interface TableRef {
  entryId: string;
  tablePath: string;
}

/** Tables of the other runs that share the table's merged picker row; none for a table without nqName. */
export function sameFuncPeers(table: TableNodeInfo, others: ComparisonEntryNodes[]): TableRef[] {
  const key = sameFuncKey(table);
  return key ? others.flatMap((other) => other.tables.filter((peer) => sameFuncKey(peer) === key)
    .map((peer) => ({entryId: other.entryId, tablePath: peer.path}))) : [];
}

/** The pick when all picks agree; undefined for no picks or a disagreement. */
export function agreedPick<T>(picks: T[], same: (a: T, b: T) => boolean = (a, b) => a === b): T | undefined {
  return picks.length > 0 && picks.every((pick) => same(pick, picks[0])) ? picks[0] : undefined;
}

const SERIES: AxisModeConfig = {mode: 'series'};
const sameAxis = (a: AxisModeConfig, b: AxisModeConfig) => a.mode === b.mode && a.units === b.units;

function annotatedPicks(table: TableNodeInfo): {index?: string, split?: string, axis?: AxisModeConfig} {
  const mode = table.defaultAxisMode;
  return {
    index: table.defaultIndexColumn,
    split: table.defaultSplitColumn,
    // units matter only on the elapsed axis; storing them elsewhere would let a
    // unit the user never saw outrank the column-meta prefill later
    axis: mode && {mode, ...mode === 'timeseries' && table.defaultTimeUnits ? {units: table.defaultTimeUnits} : {}},
  };
}

/**
 * Pickers pre-filled wherever the user has no pick yet: with the pick the same function's tables
 * in the other runs agree on, else with the tables' {comparison: {...}} annotations. The axis mode
 * comes from the runs only along with their index. Unchanged selections are returned as the same objects.
 */
export function withDefaultPicks(
  entry: ComparisonEntryNodes, selections: PickerSelections, others: ComparisonEntryNodes[] = [],
): PickerSelections {
  let {index, split, axis} = selections;
  const set = <T, >(selection: Record<string, Record<string, T>>, tablePath: string, value: T) =>
    ({...selection, [entry.entryId]: {...selection[entry.entryId], [tablePath]: value}});
  for (const table of entry.tables) {
    const columnOf = (name?: string) => table.columns.find((col) => col.name === name);
    const annotated = annotatedPicks(table);
    const peers = sameFuncPeers(table, others);
    const picksOf = <T, >(selection: Record<string, Record<string, T>>, refs = peers) =>
      refs.map(({entryId, tablePath}) => selection[entryId]?.[tablePath]).filter((pick): pick is T => !!pick);

    const peerIndex = columnOf(agreedPick(picksOf(index)))?.name;
    if (!index[entry.entryId]?.[table.path]) {
      const name = peerIndex ?? columnOf(annotated.index)?.name;
      if (name)
        index = set(index, table.path, name);
    }
    const currentIndex = index[entry.entryId]?.[table.path] ?? '';

    if (!split[entry.entryId]?.[table.path]) {
      const column = [agreedPick(picksOf(split)), annotated.split].map(columnOf)
        .find((col) => col && isSplitCandidate(col, currentIndex));
      if (column)
        split = set(split, table.path, column.name);
    }

    if (!axis[entry.entryId]?.[table.path]) {
      const indexPeers = peers.filter(({entryId, tablePath}) => index[entryId]?.[tablePath] === currentIndex);
      const peerAxis = peerIndex && currentIndex === peerIndex ?
        agreedPick(indexPeers.map(({entryId, tablePath}) => axis[entryId]?.[tablePath] ?? SERIES), sameAxis) :
        undefined;
      const config = peerAxis ?? annotated.axis;
      if (config)
        axis = set(axis, table.path, config);
    }
  }
  return {index, split, axis};
}

export function selectionToMap(
  selection: Record<string, Record<string, string>>,
  isValid: (entryId: string, tablePath: string, columnName: string) => boolean,
): Map<string, Map<string, string>> {
  const map = new Map<string, Map<string, string>>();
  for (const [entryId, tables] of Object.entries(selection)) {
    const tableMap = new Map<string, string>();
    for (const [tablePath, columnName] of Object.entries(tables)) {
      if (columnName && isValid(entryId, tablePath, columnName))
        tableMap.set(tablePath, columnName);
    }
    map.set(entryId, tableMap);
  }
  return map;
}

export interface IndexRow {
  key: string;
  entryName?: string;
  // merged rows show compare-style coverage instead of a run name
  coverage?: {count: number, total: number};
  label: string;
  title: string;
  members: {entryId: string, tablePath: string}[];
  candidates: ColumnInfo[];
  current: string;
  splitCandidates: ColumnInfo[];
  currentSplit: string;
  // set when the current index column is numeric or datetime; units for numeric only
  axis?: {mode: AxisMode, units?: TimeUnit};
}

/**
 * Rows for the index/split column pickers: one per table, or one per same-function
 * output group (by nqName + output name) when merging is on. Stored selections that
 * no longer match a candidate column are treated as unset.
 */
export function computeIndexRows(
  entries: ComparisonEntryNodes[],
  indexSelection: Record<string, Record<string, string>>,
  splitSelection: Record<string, Record<string, string>>,
  mergeSameFuncs: boolean,
  axisModeSelection: AxisModeSelection,
): IndexRow[] {
  const perTable = entries.flatMap((entry) =>
    entry.tables.map((table) => ({entry, table})));

  const validCurrent = (
    selection: Record<string, Record<string, string>>,
    entryId: string, tablePath: string, candidates: {name: string}[],
  ) => {
    const stored = selection[entryId]?.[tablePath] ?? '';
    return candidates.some((col) => col.name === stored) ? stored : '';
  };

  const splitCandidatesOf = (columns: ColumnInfo[], currentIndex: string) =>
    columns.filter((col) => isSplitCandidate(col, currentIndex));

  // an explicitly annotated default index is always offered, even for exotic types
  const indexCandidatesOf = (columns: ColumnInfo[], defaultIndexColumn?: string) =>
    columns.filter((col) => isIndexCandidateType(col.type) || col.name === defaultIndexColumn);

  // merged rows agree like current/currentSplit: a non-default mode only when every
  // member agrees on it, shared units when members agree, else metadata prefill
  // (index column of the seed table)
  const axisOf = (
    members: {entryId: string, tablePath: string}[],
    indexColumn?: ColumnInfo,
  ): IndexRow['axis'] => {
    if (!indexColumn || !isTimeIndexType(indexColumn.type))
      return undefined;
    const stored = members.map(({entryId, tablePath}) => axisModeSelection[entryId]?.[tablePath]);
    const modes = new Set(stored.map((config) => config?.mode ?? 'series'));
    const mode: AxisMode = modes.size === 1 ? [...modes][0] : 'series';
    if (!isNumericType(indexColumn.type))
      return {mode};
    const units = new Set(stored.map((config) => config?.units)
      .filter((unit): unit is TimeUnit => unit != null));
    return {
      mode,
      units: units.size === 1 ? [...units][0] : parseTimeUnit(indexColumn.units) ?? 's',
    };
  };

  const singleRow = ({entry, table}: typeof perTable[number]): IndexRow => {
    const candidates = indexCandidatesOf(table.columns, table.defaultIndexColumn);
    const current = validCurrent(indexSelection, entry.entryId, table.path, candidates);
    const splitCandidates = splitCandidatesOf(table.columns, current);
    const members = [{entryId: entry.entryId, tablePath: table.path}];
    return {
      key: `${entry.entryId}:${table.path}`,
      entryName: entry.entryName,
      label: table.friendlyPath ?? table.path,
      title: `${entry.entryName} · ${table.path}`,
      members,
      candidates,
      current,
      splitCandidates,
      currentSplit: validCurrent(splitSelection, entry.entryId, table.path, splitCandidates),
      axis: axisOf(members, candidates.find((col) => col.name === current)),
    };
  };

  if (!mergeSameFuncs)
    return perTable.map(singleRow);

  const groupKey = ({table}: typeof perTable[number]) => sameFuncKey(table);
  const groups = new Map<string, typeof perTable>();
  for (const item of perTable) {
    const key = groupKey(item);
    if (key)
      groups.set(key, [...groups.get(key) ?? [], item]);
  }

  const rows: IndexRow[] = [];
  const emitted = new Set<string>();
  for (const item of perTable) {
    const key = groupKey(item);
    const group = key ? groups.get(key)! : [item];
    if (group.length < 2) {
      rows.push(singleRow(item));
      continue;
    }
    if (emitted.has(key!))
      continue;
    emitted.add(key!);
    const sharedColumns = group[0].table.columns
      .filter((col) => group.every(({table}) =>
        table.columns.some((c) => c.name === col.name && c.type === col.type)));
    const candidates = indexCandidatesOf(sharedColumns, group[0].table.defaultIndexColumn);
    const paths = new Set(group.map(({table}) => table.path));
    const label = paths.size === 1 ?
      (group[0].table.friendlyPath ?? group[0].table.path) :
      group[0].table.name;
    const entryIds = new Set(group.map(({entry}) => entry.entryId));
    const currents = new Set(group.map(({entry, table}) =>
      validCurrent(indexSelection, entry.entryId, table.path, candidates)));
    const current = currents.size === 1 ? [...currents][0] : '';
    const splitCandidates = splitCandidatesOf(sharedColumns, current);
    const currentSplits = new Set(group.map(({entry, table}) =>
      validCurrent(splitSelection, entry.entryId, table.path, splitCandidates)));
    const members = group.map(({entry, table}) => ({entryId: entry.entryId, tablePath: table.path}));
    rows.push({
      key: `merged:${key}`,
      coverage: {count: entryIds.size, total: entries.length},
      label,
      title: group.map(({entry, table}) => `${entry.entryName} · ${table.path}`).join('\n'),
      members,
      candidates,
      current,
      splitCandidates,
      currentSplit: currentSplits.size === 1 ? [...currentSplits][0] : '',
      axis: axisOf(members, candidates.find((col) => col.name === current)),
    });
  }
  return rows;
}

export type ExclusionReason =
  | 'no similar data'
  | 'index not set'
  | 'disabled'
  | 'units differ';

export interface EntryStatus {
  entryId: string;
  matched: boolean;
  reason?: ExclusionReason;
  // the run's units when they differ from the target's
  units?: string;
  // the run still charts — its table just lacks the axis-mode settings on a partial target
  warning?: 'relative timeseries not set' | 'independent points not set';
}

export interface StatusOptions {
  axisModes?: AxisConfigMap;
  aliases?: Map<string, number>;
  convert?: UnitConverter;
}

// the run's value of the target's name (exact, normalized or mapped; never fuzzy) in units
// that neither match nor convert — the hint that a unit mapping would let it join
function differingUnits(
  entry: ComparisonEntryNodes, target: ComparisonTarget, indexColumns: Map<string, Map<string, string>>,
  {aliases, convert}: StatusOptions,
): string | undefined {
  const candidates = target.kind === 'scalar' ? entry.scalars :
    entry.tables.filter((table) => indexColumns.get(entry.entryId)?.get(table.path)).flatMap((table) => table.columns);
  return candidates.find((item) => {
    const confidence = nameMatchConfidence(target.displayName, item.name, aliases);
    return confidence && confidence !== 'fuzzy' && unitsCompatibility(target.units, item.units, convert) === 'mismatch';
  })?.units;
}

/** Per-entry participation status for a selected target. */
export function getEntryStatuses(
  entries: ComparisonEntryNodes[],
  target: ComparisonTarget | null,
  indexColumns: Map<string, Map<string, string>>,
  options: StatusOptions = {},
): EntryStatus[] {
  const {axisModes} = options;
  const modeOf = (mode: 'timeseries' | 'points') =>
    target?.kind === 'column' && axisModes ? targetAxisMode(target, axisModes, mode) : 'none';
  const partialMode = modeOf('timeseries') === 'partial' ? 'timeseries' :
    modeOf('points') === 'partial' ? 'points' : null;
  return entries.map((entry) => {
    if (!target)
      return {entryId: entry.entryId, matched: false, reason: 'no similar data'};
    const binding = (target.bindings as {entryId: string, tablePath?: string}[])
      .find((b) => b.entryId === entry.entryId);
    if (binding) {
      if (partialMode != null && binding.tablePath != null &&
        axisModes!.get(entry.entryId)?.get(binding.tablePath)?.mode !== partialMode) {
        return {entryId: entry.entryId, matched: true,
          warning: partialMode === 'timeseries' ? 'relative timeseries not set' : 'independent points not set'};
      }
      return {entryId: entry.entryId, matched: true};
    }
    if (target.kind === 'column' &&
      target.candidates.some((candidate) => candidate.binding.entryId === entry.entryId))
      return {entryId: entry.entryId, matched: false, reason: 'disabled'};
    if (target.kind === 'column' &&
      entry.tables.length > 0 &&
      entry.tables.every((table) => !indexColumns.get(entry.entryId)?.get(table.path)))
      return {entryId: entry.entryId, matched: false, reason: 'index not set'};
    const units = differingUnits(entry, target, indexColumns, options);
    if (units)
      return {entryId: entry.entryId, matched: false, reason: 'units differ', units};
    return {entryId: entry.entryId, matched: false, reason: 'no similar data'};
  });
}

/**
 * Step context for targets sharing a display name (repeated workflow steps, same-named
 * tables): the seed's step path for scalars, its table path for columns.
 */
export function targetContexts(targets: ComparisonTarget[]): Map<string, string> {
  const counts = new Map<string, number>();
  for (const target of targets)
    counts.set(target.displayName, (counts.get(target.displayName) ?? 0) + 1);
  const contexts = new Map<string, string>();
  for (const target of targets) {
    if (counts.get(target.displayName)! < 2)
      continue;
    let context: string | undefined;
    if (target.kind === 'column') {
      const seed = target.candidates.find((candidate) => candidate.auto)?.binding;
      context = seed?.tableFriendlyPath ?? seed?.tableName;
    } else {
      const seed = target.bindings[0];
      context = seed?.friendlyPath?.replace(new RegExp(` · ${escapeRegExp(seed.name)}$`), '');
    }
    if (context && context !== target.displayName)
      contexts.set(target.key, context);
  }
  return contexts;
}

const escapeRegExp = (text: string) => text.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');

export interface UnmatchedItem {
  name: string;
  units?: string;
  runCount: number;
}

/**
 * Values matching considers (numeric scalars, numeric columns of indexed tables) that no
 * compare row picks or lists as a candidate, one item per name and units. Values named like
 * a row in its units are left out: a mapping cannot place them (e.g. a step only one run has).
 */
export function unmatchedItems(
  entries: ComparisonEntryNodes[],
  targets: ComparisonTarget[],
  indexColumns: Map<string, Map<string, string>>,
): UnmatchedItem[] {
  const scalarId = (entryId: string, path: string) => `${entryId}|${path}`;
  const columnId = (entryId: string, tablePath: string, columnName: string) =>
    `${entryId}|${tablePath}|${columnName}`;
  const present = new Set(targets.flatMap((target) => target.kind === 'scalar' ?
    target.bindings.map((b) => scalarId(b.entryId, b.path)) :
    target.candidates.map(({binding: b}) => columnId(b.entryId, b.tablePath, b.columnName))));
  const nameKey = (name: string, units?: string) => `${normalizeName(name)}|${normalizeUnits(units)}`;
  const rowKeys = new Set(targets.map((target) => nameKey(target.displayName, target.units)));
  const items = new Map<string, {name: string, units?: string, runs: Set<string>}>();
  const add = (entryId: string, id: string, name: string, units?: string) => {
    const key = nameKey(name, units);
    if (present.has(id) || rowKeys.has(key))
      return;
    const item = items.get(key) ?? {name, units, runs: new Set<string>()};
    item.runs.add(entryId);
    items.set(key, item);
  };
  for (const entry of entries) {
    for (const scalar of entry.scalars) {
      if (isNumericType(scalar.valueType))
        add(entry.entryId, scalarId(entry.entryId, scalar.path), scalar.name, scalar.units);
    }
    for (const table of entry.tables) {
      const indexColumnName = indexColumns.get(entry.entryId)?.get(table.path);
      if (!indexColumnName)
        continue;
      for (const column of table.columns) {
        if (column.name !== indexColumnName && isNumericType(column.type))
          add(entry.entryId, columnId(entry.entryId, table.path, column.name), column.name, column.units);
      }
    }
  }
  return [...items.values()]
    .map(({name, units, runs}) => ({name, units, runCount: runs.size}))
    .sort((a, b) => a.name.localeCompare(b.name));
}
