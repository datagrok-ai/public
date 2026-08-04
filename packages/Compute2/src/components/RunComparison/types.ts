// Shared types for the run comparison pipeline (extraction -> matching -> selection -> builders).
// User documentation lives in help/compute/run-comparison.md

import type * as DG from 'datagrok-api/dg';

export const RUN_COLUMN = 'Run';

export type MatchConfidence = 'exact' | 'normalized' | 'fuzzy';

export interface ScalarNodeInfo {
  path: string;
  name: string;
  // human-readable path: step friendly names + IO caption; `path` stays the stable id
  friendlyPath?: string;
  valueType: string;
  units?: string;
  value: number | null;
}

export interface ColumnInfo {
  name: string;
  type: string;
  units?: string;
}

export interface TableNodeInfo {
  path: string;
  name: string;
  friendlyPath?: string;
  // nqName of the producing function; used to merge same-function tables in index selection
  nqName?: string;
  // from the {comparisonIndex: ...} / {comparisonSplit: ...} output annotations
  defaultIndexColumn?: string;
  defaultSplitColumn?: string;
  columns: ColumnInfo[];
  rowCount: number;
}

export type EntrySourceKind = 'workflow' | 'function' | 'raw';

export interface ComparisonEntryNodes {
  entryId: string;
  entryName: string;
  // 'raw' (standalone) tables get enabled in every compatible cluster, not just the best one
  sourceKind?: EntrySourceKind;
  scalars: ScalarNodeInfo[];
  tables: TableNodeInfo[];
}

export interface ComparisonEntry {
  id: string;
  name: string;
  sourceKind: EntrySourceKind;
  modelName: string;
  nodes: ComparisonEntryNodes;
  dataFrames: Map<string, DG.DataFrame>;
}

export interface ScalarBinding {
  entryId: string;
  path: string;
  name: string;
  friendlyPath?: string;
  units?: string;
  value: number | null;
}

export interface ColumnBinding {
  entryId: string;
  tablePath: string;
  tableName: string;
  tableFriendlyPath?: string;
  columnName: string;
  units?: string;
  indexColumnName: string;
  splitColumnName?: string;
}

export interface TargetBase {
  key: string;
  displayName: string;
  confidence: MatchConfidence;
  unitsWarning: boolean;
  coverage: number;
  // coverage under default enablement — a toggle-independent sort key, so rows
  // don't jump while the user (un)checks candidates
  defaultCoverage: number;
  total: number;
}

export interface ScalarTarget extends TargetBase {
  kind: 'scalar';
  bindings: ScalarBinding[];
}

export interface ColumnCandidate {
  binding: ColumnBinding;
  // confidence/unitsWarn are relative to the cluster canonical name and seed units
  confidence: MatchConfidence;
  unitsWarn: boolean;
  auto: boolean;
  enabled: boolean;
}

export interface ColumnTarget extends TargetBase {
  kind: 'column';
  candidates: ColumnCandidate[];
  // enabled candidates' bindings; what statuses, signatures and builders consume
  bindings: ColumnBinding[];
}

// index/split column names are deliberately excluded so toggles survive picker changes
export const candidateId = (b: ColumnBinding) => `${b.entryId}|${b.tablePath}|${b.columnName}`;

// targetKey -> candidateId -> enabled
export type CandidateOverrides = Record<string, Record<string, boolean>>;

export type ComparisonTarget = ScalarTarget | ColumnTarget;

const NUMERIC_TYPES = new Set(['int', 'double', 'float', 'number', 'bigint']);

export function isNumericType(type: string): boolean {
  return NUMERIC_TYPES.has(type);
}

export type TimeUnit = 'ms' | 's' | 'min' | 'h' | 'days';

export const TIME_UNITS: TimeUnit[] = ['ms', 's', 'min', 'h', 'days'];

export const TIME_UNIT_MS: Record<TimeUnit, number> =
  {ms: 1, s: 1000, min: 60_000, h: 3_600_000, days: 86_400_000};

// stored per-table pick (same Record shape as index/split selections);
// units are absent when enabled on a datetime index — a unit the user never saw
// must not outrank the metadata prefill if the index later moves to a numeric column
export interface TimeSeriesConfig {
  enabled: boolean;
  units?: TimeUnit;
}

export type TimeSeriesSelection = Record<string, Record<string, TimeSeriesConfig>>;

// resolved, enabled-only view: entryId -> tablePath; units absent <=> datetime index
export type TimeSeriesConfigMap = Map<string, Map<string, {units?: TimeUnit}>>;

export const isTimeIndexType = (type?: string) =>
  type != null && (isNumericType(type) || type === 'datetime');

export const isIndexCandidateType = (type?: string) =>
  type != null && (isNumericType(type) || type === 'string' || type === 'datetime');
