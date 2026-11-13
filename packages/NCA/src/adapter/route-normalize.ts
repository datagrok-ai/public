/**
 * Tolerant route-string → `RouteCode` normalizer.
 *
 * The data's per-row `Route` column carries free-ish strings (`IV` / `PO` /
 * `SC`, `i.v.`, NCI C66729 preferred terms, CRO abbreviations). This module
 * encodes them to the closed `nca.RouteCode` set the contract expects, with
 * three safety behaviours the closed-set {@link routeStringToCode} cannot give:
 *
 *  1. **Tolerant matching** — case / whitespace / punctuation are normalized
 *     before lookup, and a curated synonym set (CRO abbreviations + NCI
 *     Thesaurus C66729 *Route of Administration* preferred terms) is matched.
 *  2. **IV-adjacency guard (O-6)** — a string that *looks* intravascular
 *     (`iv` word, `intraven`/`infus`/`bolus`/`drip` substring) but does not map
 *     is a hard **block**, never a silent downgrade to `other`. The else branch
 *     in sci-comp `computeNca` prepends extravascular `(0,0)`; silently routing
 *     an IV-adjacent string to `other` would understate early AUC + corrupt CL
 *     with no warning. A clearly-non-IV unknown (`ip`, `topical`, …) is a
 *     *warning* + `other` (extravascular augmentation is correct there).
 *  3. **Bare-IV disambiguation** — bare `IV` /
 *     `intravenous` defaults to **IV-bolus** (small-molecule preclinical PK);
 *     when a positive infusion duration is present (per-row column or scalar)
 *     it resolves to **IV-infusion** instead. The default pushes a
 *     `BARE_IV_DEFAULTED_TO_BOLUS` warning into the per-profile warnings stream.
 *
 * `routeStringToCode` (the throwing, exact-match closed-set encoder in
 * `profile-inputs.ts`) is retained for internal closed-set use (the scalar
 * dropdown's items are always valid); only the **column** path — raw user data
 * — runs through this tolerant layer.
 */

import type * as DG from 'datagrok-api/dg';
import * as DGNS from 'datagrok-api/dg';
import {nca} from '@datagrok-libraries/sci-comp';

import type {ValidationError} from './column-contract';
import {stringifyKey} from './adapter-utils';

/** ParameterWarning code pushed when bare `IV` defaults to IV-bolus (AD-3). It
 *  rides the OPEN `nca.ParameterWarning.code` union — no sci-comp change. */
export const BARE_IV_DEFAULTED_TO_BOLUS = 'BARE_IV_DEFAULTED_TO_BOLUS';

/** Name of the transient int32 RouteCode column the mapping pane derives from
 *  the picked STRING route column and hands to `buildContract` as
 *  `observations.route` (Option A). Tilde-prefixed so it is visibly internal
 *  and unlikely to collide with a user column; excluded from group-by carry-
 *  through (see `grouping-columns.ts`). */
export const DERIVED_ROUTE_CODE_COL = '~RouteCode';

/** Normalize a raw route string for synonym lookup: lowercase, **drop dots**
 *  (so abbreviations collapse: `i.v.` → `iv`, `P.O.` → `po`, `s.c.` → `sc`),
 *  map `_`/`-`/`/` to spaces (so compound forms split: `IV-bolus` → `iv bolus`),
 *  collapse runs of whitespace, trim. Dots must NOT become spaces — otherwise
 *  `p.o.` → `p o` would miss the `po` synonym (caught by grok test 2026-06-05). */
export function normalizeRouteToken(raw: string): string {
  return raw.toLowerCase().replace(/\./g, '').replace(/[_\-/]+/g, ' ').replace(/\s+/g, ' ').trim();
}

/** Bare-IV synonyms — intravascular but bolus-vs-infusion ambiguous (AD-3). */
const BARE_IV: ReadonlySet<string> = new Set(['iv', 'intravenous', 'iv injection']);
/** Explicit IV-bolus synonyms (incl. NCI C66729 multi-word preferred terms). */
const IV_BOLUS: ReadonlySet<string> = new Set([
  'iv bolus', 'bolus', 'intravenous bolus', 'intravenous bolus injection',
  'intravenous bolus admin', 'intravenous bolus administration',
]);
/** Explicit IV-infusion synonyms. */
const IV_INFUSION: ReadonlySet<string> = new Set([
  'iv infusion', 'intravenous infusion', 'intravenous drip', 'infusion', 'ivi',
]);
const PO: ReadonlySet<string> = new Set([
  'po', 'oral', 'gavage', 'by mouth', 'oral administration', 'oral dosing', 'oral gavage',
]);
const SC: ReadonlySet<string> = new Set([
  'sc', 'subcutaneous', 'subq', 'sub q', 'subcutaneous injection', 'subcutaneous administration',
]);
const IM: ReadonlySet<string> = new Set([
  'im', 'intramuscular', 'intramuscular injection', 'intramuscular administration',
]);

/** IV-adjacency substrings (whole-word `iv` handled separately). `intraven`
 *  rather than `intra` so `intraperitoneal` / `intramuscular` do NOT false-
 *  positive into the BLOCK branch (AC-A5). */
const IV_ADJACENT_SUBSTR: ReadonlyArray<string> = ['intraven', 'infus', 'bolus', 'drip'];

function isIvAdjacent(norm: string): boolean {
  if (norm.split(' ').includes('iv')) return true;
  return IV_ADJACENT_SUBSTR.some((s) => norm.includes(s));
}

/**
 * Outcome of normalizing one route string.
 *  - `mapped`  — resolved to a `RouteCode`. `bareIvDefaulted` is true when a
 *    bare-IV string fell back to IV-bolus (warning); `infusionMissingDuration`
 *    is true when it resolved to IV-infusion but no positive duration was given.
 *  - `other`   — clearly-non-IV unknown → `ROUTE_OTHER` (warning).
 *  - `block`   — IV-adjacent unmappable → blocking error.
 *  - `empty`   — blank/whitespace cell → caller leaves `INT_NULL` (scalar fallback).
 */
export type RouteNormalization =
  | {readonly status: 'mapped'; readonly code: nca.RouteCode;
      readonly bareIvDefaulted: boolean; readonly infusionMissingDuration: boolean}
  | {readonly status: 'other'; readonly code: typeof nca.ROUTE_OTHER}
  | {readonly status: 'block'}
  | {readonly status: 'empty'};

/**
 * Normalize one raw route string. `infusionPositive` is true when a positive
 * infusion duration accompanies this row (per-row column value `> 0` or a
 * positive scalar) — it disambiguates bare IV → IV-infusion (AD-3) and clears
 * the "explicit infusion with no duration" flag.
 */
export function normalizeRouteString(
  raw: string | null | undefined,
  infusionPositive: boolean,
): RouteNormalization {
  if (raw === null || raw === undefined) return {status: 'empty'};
  const n = normalizeRouteToken(raw);
  if (n === '') return {status: 'empty'};

  if (BARE_IV.has(n)) {
    return infusionPositive ?
      {status: 'mapped', code: nca.ROUTE_IV_INFUSION, bareIvDefaulted: false, infusionMissingDuration: false} :
      {status: 'mapped', code: nca.ROUTE_IV_BOLUS, bareIvDefaulted: true, infusionMissingDuration: false};
  }
  if (IV_BOLUS.has(n))
    return {status: 'mapped', code: nca.ROUTE_IV_BOLUS, bareIvDefaulted: false, infusionMissingDuration: false};
  if (IV_INFUSION.has(n)) {
    return {status: 'mapped', code: nca.ROUTE_IV_INFUSION, bareIvDefaulted: false,
      infusionMissingDuration: !infusionPositive};
  }
  if (PO.has(n))
    return {status: 'mapped', code: nca.ROUTE_PO, bareIvDefaulted: false, infusionMissingDuration: false};
  if (SC.has(n))
    return {status: 'mapped', code: nca.ROUTE_SC, bareIvDefaulted: false, infusionMissingDuration: false};
  if (IM.has(n))
    return {status: 'mapped', code: nca.ROUTE_IM, bareIvDefaulted: false, infusionMissingDuration: false};

  // No synonym hit. IV-adjacent → BLOCK (O-6 guard); else → other + warning.
  if (isIvAdjacent(n)) return {status: 'block'};
  return {status: 'other', code: nca.ROUTE_OTHER};
}

/** Result of encoding a whole route column (AD-1 Option A). */
export interface RouteColumnResult {
  /** Per source-row `RouteCode`, `INT_NULL` for empty / blocked rows. Becomes
   *  the derived int32 column passed to `buildContract` as `observations.route`. */
  readonly codes: Int32Array;
  /** Blocking validation errors (severity = error): column-quality pre-scan
   *  failure + IV-adjacent unmappable rows. Non-empty ⇒ compute is aborted. */
  readonly errors: ValidationError[];
  /** Non-blocking per-subject warnings (keyed by stringified subject id), merged
   *  into that profile's `provenance.warnings` by the controller. */
  readonly warningsBySubject: Map<string, nca.ParameterWarning[]>;
}

/** Column-quality pre-scan threshold (AD-5, user-ratified): a route column with
 *  strictly fewer than this fraction of unique non-null values mapping to a
 *  known route synonym is blocked (`UNMAPPABLE_ROUTE_COLUMN`). Exactly at the
 *  threshold is NOT blocked (strict `<`). */
export const ROUTE_COLUMN_MAP_THRESHOLD = 0.5;

/**
 * Encode a picked STRING/INT route column to per-row `RouteCode`s, collecting
 * blocking errors and per-subject warnings (AD-1/2/3/5).
 *
 * Steps:
 *  1. **Column-quality pre-scan** — fraction of unique non-null values that
 *     resolve to a *known* route (`mapped`) must be ≥ {@link ROUTE_COLUMN_MAP_THRESHOLD};
 *     otherwise the whole column is blocked (`UNMAPPABLE_ROUTE_COLUMN`) and no
 *     per-row processing happens (the user almost certainly mapped the wrong
 *     column, e.g. SEX). `other` / `block` / `empty` do NOT count as known.
 *  2. **Per-row encode** — `normalizeRouteString` with the row's infusion-positive
 *     state. IV-adjacent unmappable → `UNMAPPABLE_ROUTE` error (row anchor);
 *     non-IV unknown → `other` + per-subject `UNMAPPABLE_ROUTE` warning;
 *     bare-IV default → `BARE_IV_DEFAULTED_TO_BOLUS` warning; explicit infusion
 *     with no positive duration → `ROUTE_INFUSION_DURATION_MISSING` warning.
 *  3. **Within-subject constancy** — a subject whose rows carry > 1 distinct
 *     resolved code gets a `ROUTE_INCONSISTENT_WITHIN_SUBJECT` warning (the
 *     extractor keeps its first-row read; this surfaces the inconsistency).
 *
 * Per-subject warnings are de-duplicated by code (one of each kind per subject).
 */
export function normalizeRouteColumn(
  routeCol: DG.Column,
  subjectCol: DG.Column,
  opts: {
    readonly infusionDurationCol?: DG.Column | null;
    readonly infusionScalarPositive?: boolean;
  } = {},
): RouteColumnResult {
  const n = routeCol.length;
  const codes = new Int32Array(n).fill(DGNS.INT_NULL);
  const errors: ValidationError[] = [];
  const warningsBySubject = new Map<string, nca.ParameterWarning[]>();

  const pushWarn = (subjectKey: string, w: nca.ParameterWarning): void => {
    let list = warningsBySubject.get(subjectKey);
    if (list === undefined) {
      list = [];
      warningsBySubject.set(subjectKey, list);
    }
    if (!list.some((e) => e.code === w.code)) list.push(w);
  };

  const durCol = opts.infusionDurationCol ?? null;
  const scalarPos = opts.infusionScalarPositive ?? false;
  const infusionPositiveAt = (row: number): boolean => {
    if (durCol !== null) {
      const v = durCol.get(row);
      if (typeof v === 'number' && Number.isFinite(v) && v > 0) return true;
    }
    return scalarPos;
  };

  // --- Step 1: column-quality pre-scan over unique non-null values. ---------
  const uniqueNorm = new Set<string>();
  for (let i = 0; i < n; i++) {
    if (routeCol.isNone(i)) continue;
    const t = normalizeRouteToken(String(routeCol.get(i)));
    if (t !== '') uniqueNorm.add(t);
  }
  if (uniqueNorm.size > 0) {
    let known = 0;
    for (const t of uniqueNorm)
      if (normalizeRouteString(t, false).status === 'mapped') known++;
    if (known / uniqueNorm.size < ROUTE_COLUMN_MAP_THRESHOLD) {
      errors.push({
        code: 'UNMAPPABLE_ROUTE_COLUMN',
        message:
          `Route column '${routeCol.name}': only ${known} of ${uniqueNorm.size} distinct ` +
          `values map to a known route. Below the ${Math.round(ROUTE_COLUMN_MAP_THRESHOLD * 100)}% ` +
          `threshold — this is probably not a route column. Map the correct column.`,
        column: routeCol.name,
      });
      return {codes, errors, warningsBySubject};
    }
  }

  // --- Step 2 + 3: per-row encode + within-subject constancy. ---------------
  const codesBySubject = new Map<string, Set<number>>();
  for (let i = 0; i < n; i++) {
    const subjectKey = stringifyKey(subjectCol.get(i));
    if (routeCol.isNone(i)) continue;
    const res = normalizeRouteString(String(routeCol.get(i)), infusionPositiveAt(i));
    switch (res.status) {
    case 'empty':
      break; // leave INT_NULL → scalar fallback
    case 'block':
      errors.push({
        code: 'UNMAPPABLE_ROUTE',
        message:
          `Route value '${routeCol.get(i)}' at row ${i} looks intravascular but does not ` +
          `map to a known route. Blocked to avoid wrong (extravascular) AUC augmentation — ` +
          `use a recognised route string (e.g. IV-bolus, IV-infusion).`,
        column: routeCol.name,
        row: i,
      });
      break;
    case 'other':
      codes[i] = res.code;
      addCode(codesBySubject, subjectKey, res.code);
      pushWarn(subjectKey, {
        code: 'UNMAPPABLE_ROUTE',
        severity: 'warning',
        message:
          `Route value '${routeCol.get(i)}' is not a recognised route; treated as ` +
          `'other' (extravascular). Verify the route string.`,
      });
      break;
    case 'mapped':
      codes[i] = res.code;
      addCode(codesBySubject, subjectKey, res.code);
      if (res.bareIvDefaulted) {
        pushWarn(subjectKey, {
          code: BARE_IV_DEFAULTED_TO_BOLUS,
          severity: 'warning',
          message:
            `Bare 'IV' route defaulted to IV-bolus (no positive infusion duration). ` +
            `Set an infusion duration to treat it as IV-infusion.`,
        });
      }
      if (res.infusionMissingDuration) {
        pushWarn(subjectKey, {
          code: 'ROUTE_INFUSION_DURATION_MISSING',
          severity: 'warning',
          message:
            `Route resolves to IV-infusion but no positive infusion duration was provided. ` +
            `Computed with extravascular augmentation (sci-comp v0.6 does not yet model ` +
            `infusion duration); supply a duration when it becomes available upstream.`,
        });
      }
      break;
    }
  }

  // Within-subject route constancy (G5): > 1 distinct resolved code per subject.
  for (const [subjectKey, set] of codesBySubject) {
    if (set.size > 1) {
      pushWarn(subjectKey, {
        code: 'ROUTE_INCONSISTENT_WITHIN_SUBJECT',
        severity: 'warning',
        message:
          `Subject '${subjectKey}' has more than one route across its rows; the profile ` +
          `uses its first row's route. Per-occasion analysis is out of scope (v1).`,
      });
    }
  }

  return {codes, errors, warningsBySubject};
}

function addCode(map: Map<string, Set<number>>, subjectKey: string, code: number): void {
  let set = map.get(subjectKey);
  if (set === undefined) {
    set = new Set<number>();
    map.set(subjectKey, set);
  }
  set.add(code);
}

/**
 * Defensive scalar-vs-column conflict check (AD-5, `ROUTE_SCALAR_COLUMN_CONFLICT`).
 *
 * The mapping pane disables the scalar route dropdown when a route column is
 * mapped, so the steady state of a correct configuration never triggers this.
 * It exists for the programmatic / inherited-scalar path: a contract that
 * arrives with *both* a route column and a disagreeing non-null scalar. The
 * column wins; this returns a one-shot warning naming the conflict, or `null`
 * when there is no conflict.
 */
export function routeScalarColumnConflict(
  columnRouteCode: nca.RouteCode | null,
  scalarRouteCode: nca.RouteCode | null,
): nca.ParameterWarning | null {
  if (columnRouteCode === null || scalarRouteCode === null) return null;
  if (columnRouteCode === scalarRouteCode) return null;
  return {
    code: 'ROUTE_SCALAR_COLUMN_CONFLICT',
    severity: 'warning',
    message:
      `Route column (code ${columnRouteCode}) disagrees with the scalar route ` +
      `(code ${scalarRouteCode}); the column wins.`,
  };
}
