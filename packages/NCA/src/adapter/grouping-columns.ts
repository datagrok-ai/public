/**
 * Grouping-column propagation (NOTE-11 gap 1).
 *
 * The per-subject results DataFrame (`writeResultsToDataFrame`) carries only
 * SUBJECT / ROUTE / STATUS + NCA parameters — no CompoundID / Dose / Sex to
 * group by. The per-compound aggregator (`aggregateNcaResults`) needs the group
 * key present as a column. This module propagates constant-per-subject metadata
 * from the source observations DataFrame into the results DataFrame.
 *
 * **Why index-lookup, not a join.** The results frame rows are in `profiles`
 * order (row i ↔ `profiles[i]`, lex-sorted by subjectKey — see
 * `writeResultsToDataFrame`). Each grouping value is constant within a subject,
 * so reading `sourceCol` at `profiles[i].observationRows[0]` is exact and O(1)
 * per subject — no key join, no ambiguity.
 *
 * **What is carried.** Every source column that is (a) NOT an analysis input
 * (subjectId / time / concentration / bodyWeight / lloq / blqFlag — these are
 * either redundant with SUBJECT or vary within a subject), (b) constant within
 * every subject, and (c) not already present in the results frame. Dose and
 * route ARE carried (valid group keys: dose proportionality, route comparison).
 *
 * **Constancy gate is load-bearing.** Without it, a sparse design with one row
 * per subject would make time / concentration falsely "constant" and they'd be
 * carried as group keys. Inputs are excluded by name AND a per-row-varying
 * metadata column (rare, but possible) is rejected by the constancy check.
 *
 * Pure w.r.t. the source DataFrame; mutates `resultsDf` in place (appends
 * columns) — the results frame is the caller's to own. Returns the names added.
 */

import type * as DG from 'datagrok-api/dg';
import * as DGNS from 'datagrok-api/dg';

import type {ProfileIndex} from './profile-index';
import type {ColumnContract} from './column-contract';

/** Build the set of source column names that are analysis inputs (never group
 *  keys). subjectId is redundant with SUBJECT; time/concentration vary within a
 *  subject; bodyWeight/lloq/blqFlag are per-row analysis inputs. Dose and route
 *  are intentionally NOT excluded — they are valid grouping keys. */
function inputColumnNames(contract: ColumnContract): ReadonlySet<string> {
  const obs = contract.observations;
  const names = [obs.subjectId.name, obs.time.name, obs.concentration.name];
  if (obs.bodyWeight !== undefined) names.push(obs.bodyWeight.name);
  if (obs.lloq !== undefined) names.push(obs.lloq.name);
  if (obs.blqFlag !== undefined) names.push(obs.blqFlag.name);
  // The encoded RouteCode column (per-row-route Option A) is an analysis input,
  // redundant with the results frame's string `ROUTE` column — never a group key.
  if (obs.route !== undefined) names.push(obs.route.name);
  return new Set(names);
}

/** True when `col` holds the same value at every row in `rows` (NA-aware: two
 *  NA cells are "equal"). Uses `getString` for a type-uniform comparison. */
function isConstantOverRows(col: DG.Column, rows: Int32Array): boolean {
  if (rows.length <= 1) return true;
  const first = col.isNone(rows[0]) ? null : col.getString(rows[0]);
  for (let j = 1; j < rows.length; j++) {
    const v = col.isNone(rows[j]) ? null : col.getString(rows[j]);
    if (v !== first) return false;
  }
  return true;
}

/** Build a per-subject column (one value per profile, read at
 *  `observationRows[0]`) preserving the source column's type. */
function buildPerSubjectColumn(
  srcCol: DG.Column, profiles: ReadonlyArray<ProfileIndex>,
): DG.Column {
  const n = profiles.length;
  const firstRowOf = (i: number) => profiles[i].observationRows[0];

  switch (srcCol.type) {
  case DGNS.COLUMN_TYPE.INT: {
    const out = new Int32Array(n);
    for (let i = 0; i < n; i++) {
      const r = firstRowOf(i);
      out[i] = srcCol.isNone(r) ? DGNS.INT_NULL : (srcCol.get(r) as number);
    }
    return DGNS.Column.fromInt32Array(srcCol.name, out);
  }
  case DGNS.COLUMN_TYPE.FLOAT: {
    const out = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      const r = firstRowOf(i);
      out[i] = srcCol.isNone(r) ? DGNS.FLOAT_NULL : (srcCol.get(r) as number);
    }
    return DGNS.Column.fromFloat64Array(srcCol.name, out);
  }
  case DGNS.COLUMN_TYPE.BOOL: {
    const out = new Array<boolean | null>(n);
    for (let i = 0; i < n; i++) {
      const r = firstRowOf(i);
      out[i] = srcCol.isNone(r) ? null : (srcCol.get(r) as boolean);
    }
    return DGNS.Column.fromList(DGNS.COLUMN_TYPE.BOOL, srcCol.name, out);
  }
  default: {
    // STRING (and any other category) → stringify per subject.
    const out = new Array<string | null>(n);
    for (let i = 0; i < n; i++) {
      const r = firstRowOf(i);
      out[i] = srcCol.isNone(r) ? null : srcCol.getString(r);
    }
    return DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, srcCol.name, out);
  }
  }
}

/**
 * Append constant-per-subject metadata columns from `sourceDf` onto `resultsDf`
 * so downstream `aggregateNcaResults` can group by them. Row alignment is by
 * index: `resultsDf` row i ↔ `profiles[i]`.
 *
 * @returns names of the columns appended (for UI "Group by" population / tests).
 * @throws if `resultsDf.rowCount !== profiles.length` (alignment broken).
 */
export function appendGroupingColumns(
  resultsDf: DG.DataFrame,
  sourceDf: DG.DataFrame,
  profiles: ReadonlyArray<ProfileIndex>,
  contract: ColumnContract,
): string[] {
  if (resultsDf.rowCount !== profiles.length) {
    throw new Error(
      `appendGroupingColumns: row/profile mismatch — results=${resultsDf.rowCount}, ` +
      `profiles=${profiles.length}. Results frame must be in profiles order.`,
    );
  }
  const excluded = inputColumnNames(contract);
  const added: string[] = [];

  for (const srcCol of sourceDf.columns as Iterable<DG.Column>) {
    if (excluded.has(srcCol.name)) continue;
    if (resultsDf.col(srcCol.name) !== null) continue; // name collision (SUBJECT/ROUTE/...)
    // Reject columns that vary within any subject — not valid group keys.
    const constant = profiles.every((p) => isConstantOverRows(srcCol, p.observationRows));
    if (!constant) continue;
    resultsDf.columns.add(buildPerSubjectColumn(srcCol, profiles));
    added.push(srcCol.name);
  }
  return added;
}

/**
 * The Report-stage Group-by candidates of a computed analysis — the ONE rule the mapping pane's
 * Group-by ribbon offers and a headless save (the demo seed) validates its `groupByCols` against:
 * the carried per-subject columns present in the observations frame, then `ROUTE` when both frames
 * can split by it (the results frame always carries the canonical `ROUTE`; `obs.col('ROUTE')`
 * resolves the source route column, DG `col()` being case-insensitive), then `CL_BASIS` when the
 * results span more than one basis and the observations carry the derived route codes.
 */
export function groupByCandidates(
  carried: ReadonlyArray<string>, results: DG.DataFrame, obs: DG.DataFrame, routeCodeCol: string,
): string[] {
  const out = carried.filter((nm) => obs.col(nm) !== null);
  if (results.col('ROUTE') !== null && obs.col('ROUTE') !== null && !out.includes('ROUTE')) out.push('ROUTE');
  const clBasis = results.col('CL_BASIS');
  if (clBasis !== null && clBasis.categories.length > 1 && obs.col(routeCodeCol) !== null &&
    !out.includes('CL_BASIS'))
    out.push('CL_BASIS');
  return out;
}
