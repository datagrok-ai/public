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

export interface ComparisonEntryNodes {
  entryId: string;
  entryName: string;
  scalars: ScalarNodeInfo[];
  tables: TableNodeInfo[];
}

export type EntrySourceKind = 'workflow' | 'function' | 'raw';

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
  total: number;
}

export interface ScalarTarget extends TargetBase {
  kind: 'scalar';
  bindings: ScalarBinding[];
}

export interface ColumnTarget extends TargetBase {
  kind: 'column';
  bindings: ColumnBinding[];
}

export type ComparisonTarget = ScalarTarget | ColumnTarget;

const NUMERIC_TYPES = new Set(['int', 'double', 'float', 'number', 'bigint']);

export function isNumericType(type: string): boolean {
  return NUMERIC_TYPES.has(type);
}
