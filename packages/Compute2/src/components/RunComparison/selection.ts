// Pure selection/view-model helpers for the RunComparison UI — index/split column pickers,
// target filtering, multi-value compatibility, per-entry statuses. No DG/platform dependencies.
// User documentation lives in help/compute/run-comparison.md

import {
  ColumnInfo, ComparisonEntryNodes, ComparisonTarget, ColumnTarget, isNumericType,
  TimeUnit, TIME_UNIT_MS, AxisMode, AxisModeSelection, AxisConfigMap,
  isTimeIndexType, isIndexCandidateType,
} from './types';
import {normalizeName, nameSimilarity, FUZZY_NAME_THRESHOLD} from './matching';

/** Substring match on the displayed text, with fuzzy per-token fallback for typos. */
export function matchesFilter(query: string, text: string): boolean {
  const q = normalizeName(query);
  if (!q)
    return true;
  const t = normalizeName(text);
  return t.includes(q) || t.split(' ').some((token) => nameSimilarity(token, q) >= FUZZY_NAME_THRESHOLD);
}

export interface MultiValueOverlap {
  // anchor runs the other target charts from the same table
  aligned: number;
  // anchor runs (entryIds) the other target has no pick for
  missing: string[];
  // anchor runs picked from a different table
  conflicting: string[];
}

interface OverlapTarget {
  bindings: {entryId: string, tablePath?: string}[];
}

// per anchor run: the other target's pick either shares the table (rows can be shared),
// is absent (padded), or points elsewhere (padded, flagged); index/split agreement is
// automatic for same-table picks since pickers are per (run, table); scalar bindings
// have no tablePath, so they are aligned or missing, never conflicting
export function multiValueOverlap(anchor: OverlapTarget, other: OverlapTarget): MultiValueOverlap {
  const otherByRun = new Map(other.bindings.map((b) => [b.entryId, b]));
  const result: MultiValueOverlap = {aligned: 0, missing: [], conflicting: []};
  for (const binding of anchor.bindings) {
    const otherBinding = otherByRun.get(binding.entryId);
    if (!otherBinding)
      result.missing.push(binding.entryId);
    else if (otherBinding.tablePath !== binding.tablePath)
      result.conflicting.push(binding.entryId);
    else
      result.aligned++;
  }
  return result;
}

/**
 * Suggestion predicate for multi-value mode. Scalar anchors: every scalar target is
 * compatible (all scalars are numeric by construction), never mixed with columns.
 * Column anchors: column targets that chart at least one of the anchor's runs from the
 * same table and conflict on none, provided the anchor's index is line-chartable
 * (numeric or datetime) everywhere. Validity of an already selected combination is
 * enforced by builder padding, not here — selected targets are never ejected.
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
  const lineIndexed = anchor.bindings.every((b) => {
    const type = getColumnType(b.entryId, b.tablePath, b.indexColumnName);
    return type != null && (isNumericType(type) || type === 'datetime');
  });
  if (!lineIndexed)
    return [];
  return targets.filter((target): target is ColumnTarget => {
    if (target.kind !== 'column')
      return false;
    const overlap = multiValueOverlap(anchor, target);
    return overlap.conflicting.length === 0 && overlap.aligned >= 1;
  });
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

/**
 * True when every bound run carries the same data for the target: scalar values or
 * element-wise value/index/split columns, with numbers compared to within
 * EQUALITY_REL_TOL relative tolerance. Unfetchable column data counts as differing.
 */
export function isTargetEqualAcrossRuns(target: ComparisonTarget, getColumnValues: ColumnValuesGetter): boolean {
  if (target.bindings.length < 2)
    return false;
  if (target.kind === 'scalar')
    return target.bindings.every((b) => sameValue(b.value, target.bindings[0].value));
  const sameValues = (a: ArrayLike<unknown> | null | undefined, b: ArrayLike<unknown> | null | undefined) => {
    if (a == null || b == null || a.length !== b.length)
      return false;
    for (let i = 0; i < a.length; i++) {
      if (!sameValue(a[i], b[i]))
        return false;
    }
    return true;
  };
  const series = target.bindings.map((b) => ({
    value: getColumnValues(b.entryId, b.tablePath, b.columnName),
    index: getColumnValues(b.entryId, b.tablePath, b.indexColumnName),
    split: b.splitColumnName ? getColumnValues(b.entryId, b.tablePath, b.splitColumnName) : null,
  }));
  const [first, ...rest] = series;
  if (!first.value || !first.index || first.split === undefined)
    return false;
  return rest.every((s) => sameValues(s.value, first.value) && sameValues(s.index, first.index) &&
    (s.split === first.split || sameValues(s.split, first.split)));
}

const UNIT_ALIASES: Record<string, TimeUnit> = {
  ms: 'ms', millis: 'ms', millisecond: 'ms', milliseconds: 'ms',
  s: 's', sec: 's', secs: 's', second: 's', seconds: 's',
  min: 'min', mins: 'min', minute: 'min', minutes: 'min',
  h: 'h', hr: 'h', hrs: 'h', hour: 'h', hours: 'h',
  d: 'days', day: 'days', days: 'days',
};

// bare 'm' is deliberately not recognized — ambiguous with meters
export function parseTimeUnit(units?: string): TimeUnit | undefined {
  return units ? UNIT_ALIASES[units.trim().toLowerCase()] : undefined;
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

/** Require-all gate: a target charts in the given mode only when every bound table is configured. */
export function targetAxisMode(
  target: ColumnTarget,
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

  const groupKey = ({table}: typeof perTable[number]) =>
    table.nqName ? `${table.nqName}|${table.path.split('/').pop()}` : null;
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
      (group[0].table.name ?? group[0].table.nqName!);
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
  | 'disabled';

export interface EntryStatus {
  entryId: string;
  matched: boolean;
  reason?: ExclusionReason;
  // the run still charts — its table just lacks the axis-mode settings on a partial target
  warning?: 'relative timeseries not set' | 'independent points not set';
}

/** Per-entry participation status for a selected target. */
export function getEntryStatuses(
  entries: ComparisonEntryNodes[],
  target: ComparisonTarget | null,
  indexColumns: Map<string, Map<string, string>>,
  axisModes?: AxisConfigMap,
): EntryStatus[] {
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
    return {entryId: entry.entryId, matched: false, reason: 'no similar data'};
  });
}
