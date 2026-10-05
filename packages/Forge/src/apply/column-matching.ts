import * as DG from 'datagrok-api/dg';
import {bigIntProblem, isReadableNumber} from '../training/default-features';
import {ColumnSchema} from '../training/train-model';

export const MAX_NAME_DISTANCE = 0.3;
/** Feature name -> table column name; a pre-encoded one-hot column maps to itself. */
export type ColumnMapping = Map<string, string>;
export type Compatibility = {kind: 'ok'} | {kind: 'hint'; message: string} | {kind: 'error'; message: string};
export interface MappingProblem { feature: string; message: string; isUnmapped: boolean }

const NUMERICAL_TYPES = ['int', 'double', 'float', 'num', 'qnum', 'bigint'];
const KIND_TEXTS: {[type: string]: string} = {
  string: 'text',
  bool: 'yes/no values',
  datetime: 'dates',
};
const INF = 1e9;

export function compatibility(feature: ColumnSchema, col: DG.Column, method: string): Compatibility {
  const semType = feature.semType;
  const isSameSemType = semType !== undefined && semType.toLowerCase() === (col.semType ?? '').toLowerCase();
  if (NUMERICAL_TYPES.includes(feature.type)) {
    if (!isReadableNumber(col)) {
      return {kind: 'error', message: col.type === DG.COLUMN_TYPE.BIG_INT ?
        bigIntProblem(col.name, method, 'choose another column') :
        `'${col.name}' holds ${kindText(col.type)} but '${feature.name}' needs numbers. Choose a column with numbers.`};
    }
    if (semType !== undefined && !isSameSemType) {
      return {kind: 'hint',
        message: `'${col.name}' is not marked as ${semType}; check that it is the same kind of value.`};
    }
    return {kind: 'ok'};
  }
  if (col.type !== feature.type) {
    const needed = kindText(feature.type);
    return {kind: 'error', message: `'${col.name}' holds ${kindText(col.type)} but '${feature.name}' needs ` +
      `${needed}. Choose a column with ${needed}.`};
  }
  if (semType !== undefined && !isSameSemType) {
    return {kind: 'error', message: `'${col.name}' is not a ${semType} column, which '${feature.name}' needs. ` +
      `Choose a ${semType} column.`};
  }
  return {kind: 'ok'};
}

export function nameDistance(feature: string, column: string): number {
  return lowerCaseDistance(feature.toLowerCase(), column.toLowerCase());
}

/** Columns named like the features (case-insensitive); without one, the feature's pre-encoded one-hot columns. */
export function exactMapping(features: ColumnSchema[], table: DG.DataFrame): ColumnMapping {
  const mapping: ColumnMapping = new Map();
  const names = table.columns.names();
  for (const feature of features) {
    const column = table.col(feature.name);
    if (column !== null)
      mapping.set(feature.name, column.name);
    else {
      const prefix = `${feature.name.toLowerCase()}=`;
      for (const encoded of names.filter((n) => n.toLowerCase().includes(prefix)))
        mapping.set(encoded, encoded);
    }
  }
  return mapping;
}

/** Exact names first, then the closest compatible names within MAX_NAME_DISTANCE; features without one stay out. */
export function suggestMapping(features: ColumnSchema[], table: DG.DataFrame): ColumnMapping {
  const mapping: ColumnMapping = new Map();
  const fits = (feature: ColumnSchema, col: DG.Column) => compatibility(feature, col, '').kind !== 'error';
  const used = new Set<string>();
  for (const feature of features) {
    const col = table.col(feature.name);
    if (col !== null && !used.has(col.name) && fits(feature, col)) {
      mapping.set(feature.name, col.name);
      used.add(col.name);
    }
  }

  const rest = features.filter((f) => !mapping.has(f.name));
  const free = table.columns.toList().filter((c) => !used.has(c.name));
  const freeNames = free.map((c) => c.name.toLowerCase());
  const cost = rest.map((f) => {
    const name = f.name.toLowerCase();
    return free.map((c, j) => {
      const distance = lowerCaseDistance(name, freeNames[j]);
      return distance <= MAX_NAME_DISTANCE && fits(f, c) ? distance : INF;
    });
  });
  const assignment = assign(cost) ?? closestFirst(cost);
  assignment.forEach((j, i) => {
    if (j >= 0)
      mapping.set(rest[i].name, free[j].name);
  });
  return mapping;
}

export function isSuggested(features: ColumnSchema[], table: DG.DataFrame): boolean {
  return suggestMapping(features, table).size === features.length;
}

/** The table columns of [feature] under [mapping]: its mapped column, or its pre-encoded one-hot columns. */
export function mappedColumns(feature: ColumnSchema, mapping: ColumnMapping): [string, string][] {
  const column = mapping.get(feature.name);
  if (column !== undefined)
    return [[feature.name, column]];
  const prefix = `${feature.name.toLowerCase()}=`;
  return [...mapping.entries()].filter(([key]) => key.toLowerCase().includes(prefix));
}

/** Per feature, in order: unmapped, missing from the table, incompatible, or a column used by an earlier feature. */
export function mappingProblems(features: ColumnSchema[], mapping: ColumnMapping, table: DG.DataFrame,
  method: string): MappingProblem[] {
  const problems: MappingProblem[] = [];
  const usedBy = new Map<string, string>();
  const add = (feature: ColumnSchema, message: string, isUnmapped: boolean = false) =>
    problems.push({feature: feature.name, message, isUnmapped});
  for (const feature of features) {
    const columns = mappedColumns(feature, mapping);
    if (columns.length === 0) {
      add(feature, `Choose a column for '${feature.name}'.`, true);
      continue;
    }
    const missing = columns.find(([, name]) => table.col(name) === null);
    if (missing !== undefined) {
      add(feature, `The column '${missing[1]}' is no longer in the table. ` +
        `Choose another column for '${feature.name}'.`);
      continue;
    }
    const [[key, name]] = columns;
    if (key !== feature.name)
      continue;
    const fit = compatibility(feature, table.getCol(name), method);
    const other = usedBy.get(name);
    if (fit.kind === 'error')
      add(feature, fit.message);
    else if (other !== undefined)
      add(feature, `'${name}' is also used for '${other}'. Choose a different column.`);
    if (other === undefined)
      usedBy.set(name, feature.name);
  }
  return problems;
}

/** Hungarian assignment of rows to columns (rows <= columns) minimizing the cost; null without a finite solution. */
function assign(cost: number[][]): number[] | null {
  const n = cost.length;
  const m = n === 0 ? 0 : cost[0].length;
  if (n > m)
    return null;
  const at = (i: number, j: number) => i === 0 || j === 0 ? INF : cost[i - 1][j - 1];
  const u = new Float64Array(n + 1);
  const v = new Float64Array(m + 1);
  const rowOf = new Int32Array(m + 1);
  const way = new Int32Array(m + 1);
  for (let i = 1; i <= n; i++) {
    rowOf[0] = i;
    let j0 = 0;
    const minDelta = new Float64Array(m + 1).fill(INF);
    const isUsed = new Uint8Array(m + 1);
    do {
      isUsed[j0] = 1;
      const i0 = rowOf[j0];
      let delta = INF;
      let j1 = -1;
      for (let j = 1; j <= m; j++) {
        if (isUsed[j])
          continue;
        const current = at(i0, j) - u[i0] - v[j];
        if (current < minDelta[j]) {
          minDelta[j] = current;
          way[j] = j0;
        }
        if (minDelta[j] < delta) {
          delta = minDelta[j];
          j1 = j;
        }
      }
      if (j1 < 0)
        return null;
      for (let j = 0; j <= m; j++) {
        if (!isUsed[j])
          minDelta[j] -= delta;
        else {
          u[rowOf[j]] += delta;
          v[j] -= delta;
        }
      }
      j0 = j1;
    } while (rowOf[j0] !== 0);
    do {
      const j1 = way[j0];
      rowOf[j0] = rowOf[j1];
      j0 = j1;
    } while (j0 !== 0);
  }
  const assignment = new Array<number>(n).fill(-1);
  for (let j = 1; j <= m; j++) {
    if (rowOf[j] === 0)
      continue;
    if (at(rowOf[j], j) >= INF)
      return null;
    assignment[rowOf[j] - 1] = j - 1;
  }
  return assignment;
}

/** Pairs by increasing cost: each row takes the closest column still free; rows without a finite cost get -1. */
function closestFirst(cost: number[][]): number[] {
  const pairs: [number, number][] = [];
  cost.forEach((row, i) => row.forEach((c, j) => {
    if (c < INF)
      pairs.push([i, j]);
  }));
  pairs.sort(([i1, j1], [i2, j2]) => cost[i1][j1] - cost[i2][j2]);
  const assignment = new Array<number>(cost.length).fill(-1);
  const taken = new Set<number>();
  for (const [i, j] of pairs) {
    if (assignment[i] < 0 && !taken.has(j)) {
      assignment[i] = j;
      taken.add(j);
    }
  }
  return assignment;
}

function lowerCaseDistance(a: string, b: string): number {
  return a.length === 1 ? DG.StringUtils.levenshteinDistance(a, b) : DG.StringUtils.jaroWinklerDistance(a, b);
}

/** What a column of [type] holds, in words: numbers, text, yes/no values, dates. */
export function kindText(type: string): string {
  return NUMERICAL_TYPES.includes(type) ? 'numbers' : KIND_TEXTS[type] ?? `${type} values`;
}
