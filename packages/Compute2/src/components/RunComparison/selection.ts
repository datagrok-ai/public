// Pure selection/view-model helpers for the RunComparison UI — index/split column pickers,
// target filtering, multi-value compatibility, per-entry statuses. No DG/platform dependencies.
// User documentation lives in help/compute/run-comparison.md

import {
  ColumnInfo, ComparisonEntryNodes, ComparisonTarget, ColumnTarget, isNumericType,
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

// over enabled bindings only (matching guarantees at most one per run)
export function bindingSignature(target: ColumnTarget): string {
  return target.bindings
    .map((b) => `${b.entryId}|${b.tablePath}|${b.indexColumnName}|${b.splitColumnName ?? ''}`)
    .sort().join(';');
}

/**
 * Column targets sharing the anchor's bindings signature (same runs, tables, index and
 * split columns), provided the anchor's index is line-chartable (numeric or datetime).
 * These are the targets a multi-value comparison can combine.
 */
export function compatibleTargetsFor(
  anchor: ComparisonTarget | null,
  targets: ComparisonTarget[],
  getColumnType: (entryId: string, tablePath: string, columnName: string) => string | undefined,
): ColumnTarget[] {
  if (!anchor || anchor.kind !== 'column')
    return [];
  const lineIndexed = anchor.bindings.every((b) => {
    const type = getColumnType(b.entryId, b.tablePath, b.indexColumnName);
    return type != null && (isNumericType(type) || type === 'datetime');
  });
  if (!lineIndexed)
    return [];
  const signature = bindingSignature(anchor);
  return targets.filter((target): target is ColumnTarget =>
    target.kind === 'column' && bindingSignature(target) === signature);
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
  isAllowedIndexType: (type?: string) => boolean,
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

  // an explicitly annotated default index is always offered, even when its type is toggled off
  const indexCandidatesOf = (columns: ColumnInfo[], defaultIndexColumn?: string) =>
    columns.filter((col) => isAllowedIndexType(col.type) || col.name === defaultIndexColumn);

  const singleRow = ({entry, table}: typeof perTable[number]): IndexRow => {
    const candidates = indexCandidatesOf(table.columns, table.defaultIndexColumn);
    const current = validCurrent(indexSelection, entry.entryId, table.path, candidates);
    const splitCandidates = splitCandidatesOf(table.columns, current);
    return {
      key: `${entry.entryId}:${table.path}`,
      entryName: entry.entryName,
      label: table.friendlyPath ?? table.path,
      title: `${entry.entryName} · ${table.path}`,
      members: [{entryId: entry.entryId, tablePath: table.path}],
      candidates,
      current,
      splitCandidates,
      currentSplit: validCurrent(splitSelection, entry.entryId, table.path, splitCandidates),
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
    rows.push({
      key: `merged:${key}`,
      coverage: {count: entryIds.size, total: entries.length},
      label,
      title: group.map(({entry, table}) => `${entry.entryName} · ${table.path}`).join('\n'),
      members: group.map(({entry, table}) => ({entryId: entry.entryId, tablePath: table.path})),
      candidates,
      current,
      splitCandidates,
      currentSplit: currentSplits.size === 1 ? [...currentSplits][0] : '',
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
}

/** Per-entry participation status for a selected target. */
export function getEntryStatuses(
  entries: ComparisonEntryNodes[],
  target: ComparisonTarget | null,
  indexColumns: Map<string, Map<string, string>>,
): EntryStatus[] {
  return entries.map((entry) => {
    if (!target)
      return {entryId: entry.entryId, matched: false, reason: 'no similar data'};
    const binding = (target.bindings as {entryId: string}[]).find((b) => b.entryId === entry.entryId);
    if (binding)
      return {entryId: entry.entryId, matched: true};
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
