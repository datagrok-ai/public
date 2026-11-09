/**
 * The passthrough (annotation) candidate gate — which source columns the Add-to-study
 * `Annotations` picker offers (v10, study-grid-passthrough-columns F1). Moved from
 * `views/nca-add-to-study.ts` (nca-hub-demo-studies) so the headless demo seed applies the SAME
 * gate without the pipeline importing the view layer. Pure.
 */

import type * as DG from 'datagrok-api/dg';

import {FIXED_RESULT_COLUMNS} from './results-writer';
import {DEFAULT_AGG_PARAMS} from './results-aggregator';
import {ACCEPTANCE_COL, ACCEPT_REASONS_COL, FLAG_COLS} from './acceptance';
import type {AnalysisRecipe, RecipeColumns} from './study-schema';

/**
 * Reserved results-frame column names — v10 passthrough candidate gate, collision half (AC-F1f).
 * `appendGroupingColumns` SILENTLY SKIPS a source column whose name already exists on the results
 * frame (`grouping-columns.ts` name-collision guard), so a source annotation column named e.g.
 * `ROUTE`/`CMAX` is NOT carried — `buildPassthroughByTier` would then read the RESERVED column
 * instead of the annotation. Excluding such names at the gate (rather than rename-on-carry, which
 * would break the study grid's typographic cross-analysis merge, R2-4) is the single mechanism.
 * DERIVED from {@link FIXED_RESULT_COLUMNS} (the `writeResultsToDataFrame` single source of truth,
 * itself build-time drift-guarded) + {@link DEFAULT_AGG_PARAMS} — NOT hand-typed, so a new result
 * column can never silently become offerable (architect finding, 2026-07-14). NOT derived from the
 * live `input.resultsDf`: that frame is POST-`appendGroupingColumns`, so a validly-carried
 * passthrough column IS already present on it — excluding present-on-resultsDf names would drop
 * every real annotation column. Case-insensitive membership (DG column names are). */
export const RESERVED_RESULT_COLS: ReadonlySet<string> = new Set<string>(
  [...FIXED_RESULT_COLUMNS, ...DEFAULT_AGG_PARAMS, ACCEPTANCE_COL, ACCEPT_REASONS_COL, ...FLAG_COLS]
    .map((n) => n.toLowerCase()));

/** The mapped-role / grouping-key source column names for this recipe — the passthrough
 *  candidate gate's exclusion set, built DIRECTLY from {@link RecipeColumns} (R1/LBC-1). NOT
 *  from `grouping-columns.ts` `inputColumnNames`, which draws from `ColumnContract.observations`
 *  and omits `doseCol`/`sexCol`/`infusionDurationCol` — so a mapped-but-ungrouped `sexCol` would
 *  leak as a passthrough candidate and double-render (once as its E4 keyInventory counted-join,
 *  once as a bare-when-unique passthrough). Every populated role field + `groupByCols`. */
export function recipeRoleNames(recipe: AnalysisRecipe): ReadonlySet<string> {
  const c: RecipeColumns = recipe.columns;
  const names = [
    c.subject, c.time, c.concentration, c.doseCol, c.routeCol, c.infusionDurationCol,
    c.bodyWeightCol, c.lloqCol, c.blqFlagCol, c.sexCol,
  ].filter((x): x is string => x !== null);
  return new Set<string>([...names, ...recipe.groupByCols]);
}

/** True when `col` holds ≤1 distinct non-null value within every subject (NA-aware), grouping
 *  `df` rows by `subjectCol`. The per-subject constancy gate (reuses the shipped
 *  `isConstantOverRows` semantics, `grouping-columns.ts:55` — but over `df` subject groups, since
 *  the dialog holds only the source frame). A per-observation column (a per-timepoint flag) that
 *  varies within a subject is NOT a valid passthrough (it has no single per-subpopulation value).
 *  When `subjectCol` is null/absent the whole frame is one group → the column must be globally
 *  constant (conservative). **NO cardinality cap** — a SMILES column (~1 distinct/subject, ≫20
 *  across a cassette) MUST stay offerable (RG-D). */
export function isConstantWithinSubject(
  col: DG.Column, subjectCol: DG.Column | null,
): boolean {
  const n = col.length;
  if (subjectCol === null) {
    let first: string | null = null; let seen = false;
    for (let i = 0; i < n; i++) {
      if (col.isNone(i)) continue;
      const v = col.getString(i);
      if (!seen) {first = v; seen = true;} else if (v !== first) return false;
    }
    return true;
  }
  const perSubject = new Map<string, string | null>();
  for (let i = 0; i < n; i++) {
    if (col.isNone(i)) continue;
    const s = subjectCol.isNone(i) ? '' : subjectCol.getString(i);
    const v = col.getString(i);
    if (!perSubject.has(s)) perSubject.set(s, v);
    else if (perSubject.get(s) !== v) return false;
  }
  return true;
}

/**
 * The passthrough candidate predicate (v10, study-grid-passthrough-columns F1) — a source
 * column is offerable in the add-to-study "Additional columns" multiselect iff ALL hold, with
 * **NO cardinality cap** (RG-D, the one place passthrough diverges from E4):
 *  1. NOT a mapped role / grouping key (exclusion set from {@link recipeRoleNames} — RecipeColumns
 *     roles + groupByCols, R1/LBC-1) and NOT an NCA result parameter;
 *  2. constant-within-subject ({@link isConstantWithinSubject}) — a per-observation column that
 *     varies within a subject is not a valid annotation;
 *  3. NOT a name-collision with a reserved results-frame column ({@link RESERVED_RESULT_COLS},
 *     R1/LBC-3 + R2-4) — such a source column is silently dropped by `appendGroupingColumns`, so
 *     it is EXCLUDED (never renamed-on-carry, which would break the grid's typographic merge).
 * Pure + exported so the cap-free gate is unit-testable in isolation (the RG-D proof).
 */
export function passthroughCandidate(
  col: DG.Column, recipe: AnalysisRecipe, subjectCol: DG.Column | null,
): boolean {
  const name = col.name;
  const lower = name.toLowerCase();
  if (recipeRoleNames(recipe).has(name)) return false;
  if (RESERVED_RESULT_COLS.has(lower)) return false;
  if (DEFAULT_AGG_PARAMS.some((p) => p.toLowerCase() === lower)) return false;
  return isConstantWithinSubject(col, subjectCol);
}

