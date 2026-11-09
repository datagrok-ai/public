/**
 * Row ↔ domain mapping of the EMS schema `nca` 0.2.0: pure functions, each in
 * both directions, so {@link EmsNcaRepository} holds no mapping of its own and the round trips are
 * unit-tested without a server.
 *
 * The domain types (`AnalysisEntry`, `SummaryCache`, `LzOverride`, `ReviewLedger`, `PipelineFlag`,
 * `AnalysisOrigin`) stay as they are; only their storage changes:
 *  - fixed-field objects → typed columns on `analysis` (signature, stage, origin, confirmed, tier);
 *  - the results cache → `analysis_result` rows (one per `CompoundTierRow`);
 *  - λz overrides, review decisions, pipeline flags → their own tables;
 *  - the run part of a pipeline origin → `pipeline_run` (one row per run and study).
 *
 * Read-back normalisation, measured on the stand (`NCA EMS: platform`): a datetime comes back as a
 * Dayjs, an unset `string_list` as null, NaN never (JSON carries it as null), a float as the exact
 * double. Every reader here maps those back, and every value a diff compares goes through
 * {@link canonical} with the same normalisation on both sides.
 */

import type {
  AnalysisEntry, AnalysisRecipe, AnalysisOrigin, CompoundCurve, CompoundTierRow, ConfirmedRecipe, LzOverride,
  MeanCv, PipelineFlag, PipelineOrigin, ReviewDecisionRecord, ReviewLedger, SourceSignature, StageCache,
  SummaryCache,
} from './study-schema';
import {recipeOf} from './study-schema';
import {computeRecipeFingerprint} from './fingerprint';

// eslint-disable-next-line @typescript-eslint/no-explicit-any
export type Bag = Record<string, any>;

// ─── Scalars ────────────────────────────────────────────────────────────────────

/** A stored float: null (how JSON carries NaN) reads back as NaN. */
const num = (v: unknown): number => (v === null || v === undefined ? NaN : Number(v));
/** A float to store: a non-finite value is null (JSON cannot carry it). */
const finite = (x: number): number | null => (Number.isFinite(x) ? x : null);
/** A stored `string_list`: an unset list reads back as null. */
const list = (v: unknown): string[] => (Array.isArray(v) ? v.map((x) => String(x)) : []);

/** A datetime as ISO-8601 with milliseconds (`toISOString`): a Dayjs, a Date or a parseable string;
 *  an unparseable string is kept verbatim (shown, never dropped); null / undefined / '' → null. */
export function isoOf(d: unknown): string | null {
  if (d === null || d === undefined || d === '') return null;
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const t = (d as any).toISOString;
  if (typeof t === 'function') return t.call(d) as string;
  const parsed = new Date(String(d));
  return Number.isNaN(parsed.getTime()) ? String(d) : parsed.toISOString();
}

/** JSON with sorted keys and `undefined` dropped — json columns come back with their keys reordered. */
export function canonical(x: unknown): string {
  const walk = (v: unknown): unknown => {
    if (Array.isArray(v)) return v.map(walk);
    if (v !== null && typeof v === 'object') {
      return Object.fromEntries(Object.keys(v as Bag).sort().filter((k) => (v as Bag)[k] !== undefined)
        .map((k) => [k, walk((v as Bag)[k])]));
    }
    return v;
  };
  return JSON.stringify(walk(x) ?? null);
}

// ─── Column sets and their read-back normalisation ──────────────────────────────

const DATETIME_COLUMNS = new Set(['results_computed_at', 'source_mtime', 'source_computed_at', 'confirmed_at',
  'updated_at', 'set_at', 'decided_at', 'superseded_at', 'run_at']);
const LIST_COLUMNS = new Set(['group_keys', 'group_values', 'warning_codes', 'warning_messages', 'reason_codes',
  'reason_texts', 'reasons', 'info', 'pipeline_only']);

/** One column value as a diff compares it: a datetime as ISO, a list with null ≡ [], absent ≡ null,
 *  and an empty string ≡ null — EMS stores a scalar `''` as null (measured: a written
 *  `tier_flag_summary: ''` reads back null), so without this every reopen would see a change. Items
 *  of a list keep `''` (the `string_list` probe). */
function norm(column: string, v: unknown): unknown {
  if (DATETIME_COLUMNS.has(column)) return isoOf(v);
  if (LIST_COLUMNS.has(column)) return list(v);
  return v === undefined || v === '' ? null : v;
}

/** `columns` of `row`, each normalised and canonical — the stored state a later write diffs against. */
export function canonicalColumns(row: Bag, columns: readonly string[]): Record<string, string> {
  return Object.fromEntries(columns.map((c) => [c, canonical(norm(c, row[c]))]));
}

/** The canonical form of a whole row over `columns` (one string). */
export function canonicalRow(row: Bag, columns: readonly string[]): string {
  return canonical(columns.map((c) => norm(c, row[c])));
}

/**
 * The `analysis` columns a domain entry fully determines — every write diffs these against the
 * stored row and sends only what changed. Not here: `study_id`, `moved_from_label`, `pipeline_run`
 * (set by the write that owns them), `updated_at` and `results_computed_at` (set by the writer).
 */
export const ENTRY_COLUMNS = [
  'label', 'engine_version', 'results_basis', 'has_cache',
  'tier_n_subjects', 'tier_n_compounds', 'tier_route', 'tier_dose', 'tier_status', 'tier_flag_summary',
  'results_snapshot_hash', 'group_keys',
  'row_count', 'source_kind', 'source_file', 'source_ref',
  'source_sig_kind', 'source_mtime', 'source_hash', 'source_computed_at',
  'recipe', 'recipe_version',
  'is_sparse', 'sparse_sampling_type', 'auc_method', 'blq_pre_first', 'blq_embedded', 'blq_after_last',
  'blq_consecutive',
  'confirmed_fingerprint', 'confirmed_at', 'confirmed_recipe',
  'stage', 'stage_pending', 'stage_not_applicable',
  'origin_kind', 'origin_label', 'origin_container_id',
  'pipeline_partition', 'pipeline_recipe_hash', 'pipeline_recipe_hash_version', 'pipeline_supersedes',
  'pipeline_flag_count',
] as const;

/** The columns of a results write: changing any of them changes the results (decision 2). */
export const RESULT_ENTRY_COLUMNS: ReadonlySet<string> = new Set(['engine_version', 'results_basis', 'has_cache',
  'tier_n_subjects', 'tier_n_compounds', 'tier_route', 'tier_dose', 'tier_status', 'tier_flag_summary',
  'results_snapshot_hash', 'group_keys']);

// ─── Recipe search columns (hybrid, decision 7) ─────────────────────────────────

/** The derived search columns of a recipe; never read back into it. Optional chaining: a test or a
 *  pre-v13 recipe may lack a block, which then stores null. */
export function recipeSearchCols(recipe: AnalysisRecipe): Bag {
  const sparse = recipe.modes?.sparse;
  const blq = recipe.ncaRules?.blq;
  return {
    is_sparse: sparse?.on === true,
    sparse_sampling_type: sparse?.on === true ? (sparse.samplingType ?? null) : null,
    auc_method: recipe.ncaRules?.aucMethod ?? null,
    blq_pre_first: blq?.preFirstMeasurable ?? null,
    blq_embedded: blq?.embedded ?? null,
    blq_after_last: blq?.afterLast ?? null,
    blq_consecutive: blq?.consecutiveAfterLast ?? null,
  };
}

// ─── Source signature ───────────────────────────────────────────────────────────

export function signatureCols(s: SourceSignature): Bag {
  return {source_sig_kind: s.kind, source_mtime: isoOf(s.mtime), source_hash: s.hash,
    source_computed_at: isoOf(s.computedAt), row_count: s.rowCount};
}

/** A row without a signature reads as the weakest one, stamped at the row's creation (the reading
 *  of an absent signature before 0.2.0). */
export function signatureOf(r: Bag): SourceSignature {
  return {
    kind: r.source_sig_kind ?? 'computed-at-only', rowCount: r.row_count ?? 0, mtime: isoOf(r.source_mtime),
    hash: r.source_hash ?? null, computedAt: isoOf(r.source_computed_at) ?? isoOf(r.created_on) ?? '',
  };
}

// ─── Stage ──────────────────────────────────────────────────────────────────────

export function stageCols(s: StageCache | null): Bag {
  return {stage: s?.stage ?? null, stage_pending: s?.pending ?? null, stage_not_applicable: s?.notApplicable ?? null};
}

export function stageOf(r: Bag): StageCache | null {
  if (r.stage === null || r.stage === undefined) return null;
  return {stage: r.stage, pending: r.stage_pending ?? 0, notApplicable: r.stage_not_applicable ?? 0};
}

// ─── Confirmed recipe (decision 4) ──────────────────────────────────────────────

/** `confirmed_recipe` is written only when the snapshot's fingerprint differs from the stored
 *  recipe's (an unconfirmed Setup edit that was saved, P6); otherwise null and rebuilt on read. */
export function confirmedCols(c: ConfirmedRecipe | null, recipeFingerprint: string): Bag {
  if (c === null) return {confirmed_fingerprint: null, confirmed_at: null, confirmed_recipe: null};
  return {
    confirmed_fingerprint: c.fingerprint, confirmed_at: isoOf(c.at),
    confirmed_recipe: c.fingerprint === recipeFingerprint ? null : c.recipe,
  };
}

export function confirmedOf(r: Bag, recipe: AnalysisRecipe): ConfirmedRecipe | null {
  if (r.confirmed_fingerprint === null || r.confirmed_fingerprint === undefined) return null;
  return {recipe: (r.confirmed_recipe ?? recipe) as AnalysisRecipe, fingerprint: r.confirmed_fingerprint,
    at: isoOf(r.confirmed_at) ?? ''};
}

// ─── Origin (standalone / pipeline) + pipeline_run ──────────────────────────────

/** The per-analysis origin columns. A pipeline origin's run part goes to `pipeline_run`
 *  ({@link pipelineRunCols}); the `pipeline_run` ref is set by the write that owns the run row. */
export function originCols(o: AnalysisOrigin | null): Bag {
  const none = {origin_kind: null, origin_label: null, origin_container_id: null, pipeline_partition: null,
    pipeline_recipe_hash: null, pipeline_recipe_hash_version: null, pipeline_supersedes: null};
  if (o === null) return none;
  if (o.kind === 'standalone') {
    return {...none, origin_kind: 'standalone', origin_label: o.label,
      origin_container_id: o.containerId};
  }
  return {...none, origin_kind: 'pipeline', pipeline_partition: o.partition, pipeline_recipe_hash: o.recipeHash,
    pipeline_recipe_hash_version: o.recipeHashVersion, pipeline_supersedes: o.supersedes ?? null};
}

/** The `pipeline_run` row of a pipeline origin in `studyId`. `mapping` is stored as JSON TEXT (a
 *  string column keeps it verbatim; jsonb would reorder its keys — decision 13). */
export function pipelineRunCols(o: PipelineOrigin, studyId: string): Bag {
  return {study_id: studyId, run_id: o.runId, run_at: isoOf(o.runAt), run_user: o.user,
    engine_version: o.engineVersion, rules: o.rules, mapping: JSON.stringify(o.mapping)};
}

/** The origin of an analysis row read with `expand: ['pipeline_run']` (the run columns arrive as
 *  `pipeline_run.<column>`). A pipeline row whose run is gone (setnull) keeps its own fields and
 *  reads empty run fields. */
export function originOf(r: Bag): AnalysisOrigin | null {
  if (r.origin_kind === 'standalone')
    return {kind: 'standalone', label: r.origin_label ?? '', containerId: r.origin_container_id ?? ''};
  if (r.origin_kind !== 'pipeline') return null;
  const run = (c: string): unknown => r[`pipeline_run.${c}`];
  let mapping: Record<string, unknown> = {};
  const text = run('mapping');
  if (typeof text === 'string' && text !== '')
    try {mapping = JSON.parse(text) as Record<string, unknown>;} catch {/* not ours to repair */}

  return {
    kind: 'pipeline', runId: (run('run_id') as string | null) ?? '', runAt: isoOf(run('run_at')) ?? '',
    user: (run('run_user') as string | null) ?? '', engineVersion: (run('engine_version') as string | null) ?? '',
    mapping, rules: (run('rules') as string | null) ?? '', partition: r.pipeline_partition ?? '',
    recipeHash: r.pipeline_recipe_hash ?? '', recipeHashVersion: r.pipeline_recipe_hash_version ?? 0,
    ...(r.pipeline_supersedes != null ? {supersedes: r.pipeline_supersedes as string} : {}),
  };
}

// ─── Results summary: analysis tier + analysis_result rows (decisions 2, 5, 8, 15) ────

/** The MeanCv parameters of a tier row and their column prefixes (`%F` has its own null rule). */
const MEAN_CV_PARAMS: ReadonlyArray<readonly [keyof CompoundTierRow, string]> = [
  ['cmax', 'cmax'], ['tmax', 'tmax'], ['aucInf', 'auc_inf'], ['thalf', 'thalf'], ['cl', 'cl'], ['vz', 'vz'],
  ['mrt', 'mrt'], ['vss', 'vss'], ['pctExtrap', 'pct_extrap'],
];

/** Every `analysis_result` column the app writes, except `analysis_id` (the diff of a results write). */
export const RESULT_COLUMNS: readonly string[] = [
  'study_id', 'ord', 'compound_id', 'is_whole_analysis', 'route', 'dose_norm', 'n', 'engine_version',
  'results_basis',
  ...[...MEAN_CV_PARAMS.map(([, p]) => p), 'f_pct'].flatMap((p) => [`${p}_mean`, `${p}_cv`, `${p}_n`,
    `${p}_n_excluded`]),
  'f_pct_ne', 'has_curve', 'group_values', 'warning_codes', 'warning_messages', 'structure', 'composition', 'curve',
];

/** The columns the study overview reads (all but the curve, loaded lazily). */
export const RESULT_COLUMNS_WITHOUT_CURVE: readonly string[] = ['analysis_id',
  ...RESULT_COLUMNS.filter((c) => c !== 'curve')];

function meanCvCols(p: string, m: MeanCv | null): Bag {
  return {[`${p}_mean`]: m === null ? null : finite(m.mean), [`${p}_cv`]: m === null ? null : finite(m.cvPct),
    [`${p}_n`]: m === null ? null : m.n, [`${p}_n_excluded`]: m?.nExcluded ?? null};
}

/** A stored MeanCv: a null mean / %CV reads back as NaN — the not-estimated value, and with `n = 0`
 *  the all-excluded case FIELD-02 renders `not reportable`. */
function meanCvOf(r: Bag, p: string): MeanCv {
  const ex = r[`${p}_n_excluded`];
  return {mean: num(r[`${p}_mean`]), cvPct: num(r[`${p}_cv`]), n: r[`${p}_n`] ?? 0,
    ...(ex !== null && ex !== undefined ? {nExcluded: ex as number} : {})};
}

/** The group key labels of a summary: one key list per analysis (every grouped row carries the same
 *  keys in pick order; a whole-analysis row carries none, `[]`). Rows that disagree are a broken invariant
 *  of `buildSummaryCache` and refuse to store rather than mislabel a column. */
export function groupKeysOf(cache: SummaryCache): string[] {
  let keys: string[] | null = null;
  for (const r of cache.compoundTier) {
    if (r.groupValues.length === 0) continue;
    const k = r.groupValues.map((g) => g.key);
    if (keys === null) keys = k;
    else if (canonical(keys) !== canonical(k))
      throw new Error(`The result rows carry different group keys (${keys.join(', ')} vs ${k.join(', ')}).`);
  }
  return keys ?? [];
}

/** The analysis-tier columns of a summary (null cache → all null, `has_cache` false). */
export function tierCols(cache: SummaryCache | null): Bag {
  const t = cache?.analysisTier ?? null;
  return {
    has_cache: cache !== null,
    tier_n_subjects: t?.nSubjects ?? null, tier_n_compounds: t?.nCompounds ?? null, tier_route: t?.route ?? null,
    tier_dose: t?.dose ?? null, tier_status: t?.status ?? null, tier_flag_summary: t?.flagSummary ?? null,
    group_keys: cache === null ? [] : groupKeysOf(cache),
  };
}

/** The SMILES of a row with exactly one value in a structure annotation (decision 8); else null. */
function structureOf(r: CompoundTierRow): string | null {
  const e = r.passthrough.find((p) => p.isStructure);
  return e !== undefined && e.values.length === 1 ? e.values[0].value : null;
}

/** The `analysis_result` rows of a summary, in tier order. The staleness stamps are copies of the
 *  analysis's, written in the same transaction (decision 15). */
export function resultRowsOf(analysisId: string, studyId: string, cache: SummaryCache,
  stamps: {readonly engineVersion: string | null; readonly resultsBasis: string | null}): Bag[] {
  return cache.compoundTier.map((r, ord) => ({
    analysis_id: analysisId, study_id: studyId, ord,
    compound_id: r.compoundId, is_whole_analysis: r.isWholeAnalysis, route: r.route, dose_norm: r.doseNorm, n: r.n,
    engine_version: stamps.engineVersion, results_basis: stamps.resultsBasis,
    ...Object.assign({}, ...MEAN_CV_PARAMS.map(([k, p]) => meanCvCols(p, r[k] as MeanCv))),
    ...meanCvCols('f_pct', r.fPct),
    f_pct_ne: r.fPct !== null && Number.isNaN(r.fPct.mean),
    has_curve: r.curve !== null,
    group_values: r.groupValues.map((g) => g.value),
    warning_codes: r.warnings.map((w) => w.code), warning_messages: r.warnings.map((w) => w.message),
    structure: structureOf(r),
    composition: {keyInventory: r.keyInventory, passthrough: r.passthrough},
    curve: r.curve,
  }));
}

/** One tier row of a stored result row. `curve`: the row's own column when it was read (a light
 *  read leaves it `undefined`; a stored no-curve is `null`), else the
 *  lazily loaded one by row id, else null. */
function tierRowOf(r: Bag, groupKeys: readonly string[],
  curves: ReadonlyMap<string, CompoundCurve | null> | undefined): CompoundTierRow {
  const values = list(r.group_values);
  const codes = list(r.warning_codes);
  const messages = list(r.warning_messages);
  const composition = (r.composition ?? {}) as Bag;
  const fPctApplies = r.f_pct_n !== null && r.f_pct_n !== undefined;
  return {
    compoundId: r.compound_id ?? '', isWholeAnalysis: r.is_whole_analysis === true,
    groupValues: values.map((value, i) => ({key: groupKeys[i] ?? '', value})),
    keyInventory: composition.keyInventory ?? [], passthrough: composition.passthrough ?? [],
    n: r.n ?? 0, route: r.route ?? '', doseNorm: r.dose_norm ?? 'absolute',
    cmax: meanCvOf(r, 'cmax'), tmax: meanCvOf(r, 'tmax'), aucInf: meanCvOf(r, 'auc_inf'),
    thalf: meanCvOf(r, 'thalf'), cl: meanCvOf(r, 'cl'), vz: meanCvOf(r, 'vz'), mrt: meanCvOf(r, 'mrt'),
    vss: meanCvOf(r, 'vss'),
    fPct: !fPctApplies ? null : r.f_pct_ne === true ? {...meanCvOf(r, 'f_pct'), mean: NaN} : meanCvOf(r, 'f_pct'),
    pctExtrap: meanCvOf(r, 'pct_extrap'),
    curve: r.curve !== undefined ? ((r.curve ?? null) as CompoundCurve | null) : (curves?.get(r.id) ?? null),
    warnings: codes.map((code, i) => ({code, message: messages[i] ?? code})),
  };
}

/**
 * The summary of an analysis row and its result rows (any order; sorted by `ord`). `has_cache`
 * false → null (never confirmed, FIELD-04); true with zero rows → an empty tier, not null.
 */
export function summaryCacheOf(a: Bag, rows: readonly Bag[],
  curves?: ReadonlyMap<string, CompoundCurve | null>): SummaryCache | null {
  if (a.has_cache !== true) return null;
  const keys = list(a.group_keys);
  return {
    analysisTier: {
      nCompounds: a.tier_n_compounds ?? 0, nSubjects: a.tier_n_subjects ?? 0, route: a.tier_route ?? '',
      dose: a.tier_dose ?? '', status: a.tier_status ?? 'ok', flagSummary: a.tier_flag_summary ?? '',
    },
    compoundTier: [...rows].sort((x, y) => (x.ord ?? 0) - (y.ord ?? 0)).map((r) => tierRowOf(r, keys, curves)),
  };
}

// ─── λz overrides → lz_override (decision 3) ────────────────────────────────────

/** The business key of one override point. */
export const lzKey = (profileKey: string, timeKey: string, ordinal: number): string =>
  `${profileKey}\u0001${timeKey}\u0001${ordinal}`;

/** The columns of an `lz_override` row the app owns (the diff of an override write). */
export const LZ_COLUMNS: readonly string[] = ['profile_key', 'time_key', 'ordinal', 'conc_key', 'action',
  'set_by', 'set_by_label', 'set_at'];

/** `set_by` / `set_at` are explicit (not the system author / created_on): a save beside copies the
 *  rows and must keep the analyst and the time of each point (decision 10). `set_by` is the user
 *  reference (null without a platform user), `set_by_label` the name as recorded. */
export function lzRowsOf(analysisId: string, overrides: ReadonlyArray<LzOverride>): Bag[] {
  return overrides.map((o) => ({analysis_id: analysisId, profile_key: o.profileKey, time_key: o.timeKey,
    ordinal: o.ordinal, conc_key: o.concKey, action: o.action, set_by: o.userId ?? null, set_by_label: o.user,
    set_at: isoOf(o.at)}));
}

/** The domain form of a set of overrides, normalised (ordered by key, times as ISO) — what a stored
 *  set reads back as, so a content comparison does not depend on order or time format. */
export function normLz(overrides: ReadonlyArray<LzOverride>): LzOverride[] {
  return overrides.map((o) => ({...o, at: isoOf(o.at) ?? ''}))
    .sort((x, y) => cmp(lzKey(x.profileKey, x.timeKey, x.ordinal), lzKey(y.profileKey, y.timeKey, y.ordinal)));
}

export function lzOverridesOf(rows: readonly Bag[]): LzOverride[] {
  return normLz(rows.map((r) => ({profileKey: r.profile_key, timeKey: r.time_key, ordinal: r.ordinal ?? 0,
    concKey: r.conc_key ?? '', action: r.action, user: r.set_by_label ?? '',
    ...(r.set_by != null ? {userId: r.set_by as string} : {}), at: isoOf(r.set_at) ?? ''})));
}

// ─── Review decisions → review_decision, append-only (decisions 3, 10) ──────────

const cmp = (a: string, b: string): number => (a < b ? -1 : a > b ? 1 : 0);

/** The immutable core of one decision record — its identity in the append-only table. */
function decisionCore(unitKey: string, d: ReviewDecisionRecord): Bag {
  return {unit_key: unitKey, action: d.action, reason_codes: d.reasons.map((x) => x.code),
    reason_texts: d.reasons.map((x) => x.text), decided_by: d.userId ?? null, decided_by_label: d.user,
    decided_at: isoOf(d.at),
    fingerprint: d.fingerprint};
}

const DECISION_CORE: readonly string[] = ['unit_key', 'action', 'reason_codes', 'reason_texts', 'decided_by',
  'decided_by_label', 'decided_at', 'fingerprint'];

/** The identity of a stored decision row (or of the row a record would be). */
export const decisionKey = (row: Bag): string => canonicalRow(row, DECISION_CORE);

/** Every record of every ledger as a row: the current one without `superseded_*`, each history
 *  record with them. */
export function decisionRowsOf(analysisId: string, studyId: string, ledgers: ReadonlyArray<ReviewLedger>): Bag[] {
  const out: Bag[] = [];
  for (const l of ledgers) {
    for (const d of [...l.history, ...(l.current !== null ? [l.current] : [])]) {
      out.push({analysis_id: analysisId, study_id: studyId, ...decisionCore(l.unitKey, d),
        superseded_by: d.supersededBy ?? null, superseded_at: isoOf(d.supersededAt)});
    }
  }
  return out;
}

const byTime = (x: ReviewDecisionRecord, y: ReviewDecisionRecord): number =>
  cmp(x.at, y.at) || cmp(x.supersededAt ?? '', y.supersededAt ?? '') || cmp(canonical(x), canonical(y));

/** Ledgers normalised as a stored set reads back: ordered by unit, history oldest first, times as
 *  ISO, units with no record dropped (a ledger with nothing in it stores nothing). */
export function normLedgers(ledgers: ReadonlyArray<ReviewLedger>): ReviewLedger[] {
  const rec = (d: ReviewDecisionRecord): ReviewDecisionRecord => ({
    action: d.action, reasons: d.reasons.map((x) => ({code: x.code, text: x.text})), user: d.user,
    ...(d.userId !== undefined ? {userId: d.userId} : {}), at: isoOf(d.at) ?? '', fingerprint: d.fingerprint,
    ...(d.supersededBy !== undefined ? {supersededBy: d.supersededBy} : {}),
    ...(d.supersededAt !== undefined ? {supersededAt: isoOf(d.supersededAt) ?? ''} : {}),
  });
  return ledgers.filter((l) => l.current !== null || l.history.length > 0)
    .map((l) => ({unitKey: l.unitKey, current: l.current === null ? null : rec(l.current),
      history: l.history.map(rec).sort(byTime)}))
    .sort((x, y) => cmp(x.unitKey, y.unitKey));
}

/** The ledgers of stored rows: grouped by unit; `current` = the row without `superseded_by` (the
 *  latest, should there ever be two), `history` = the others oldest first. */
export function ledgersOf(rows: readonly Bag[]): ReviewLedger[] {
  const byUnit = new Map<string, ReviewDecisionRecord[]>();
  const current = new Map<string, ReviewDecisionRecord[]>();
  for (const r of rows) {
    const codes = list(r.reason_codes);
    const texts = list(r.reason_texts);
    const d: ReviewDecisionRecord = {
      action: r.action, reasons: codes.map((code, i) => ({code, text: texts[i] ?? ''})),
      user: r.decided_by_label ?? '', ...(r.decided_by != null ? {userId: r.decided_by as string} : {}),
      at: isoOf(r.decided_at) ?? '', fingerprint: r.fingerprint ?? '',
      ...(r.superseded_by != null ? {supersededBy: r.superseded_by} : {}),
      ...(r.superseded_by != null && r.superseded_at != null ? {supersededAt: isoOf(r.superseded_at) ?? ''} : {}),
    };
    const into = r.superseded_by == null ? current : byUnit;
    into.set(r.unit_key, [...(into.get(r.unit_key) ?? []), d]);
  }
  const units = [...new Set([...byUnit.keys(), ...current.keys()])];
  return normLedgers(units.map((unitKey) => {
    const cur = [...(current.get(unitKey) ?? [])].sort(byTime);
    return {unitKey, current: cur.length > 0 ? cur[cur.length - 1] : null,
      history: [...(byUnit.get(unitKey) ?? []), ...cur.slice(0, -1)]};
  }));
}

/** A stored decision row as the tracker keeps it. */
export interface StoredDecision {
  readonly id: string;
  readonly key: string;
  readonly supersededBy: string | null;
  readonly supersededAt: string | null;
}

export function storedDecisionOf(r: Bag): StoredDecision {
  return {id: r.id, key: decisionKey(r), supersededBy: r.superseded_by ?? null, supersededAt: isoOf(r.superseded_at)};
}

/** One record of the target ledgers against the stored rows: its row, the stored row it matched
 *  (by the immutable core) and the op it needs (null = stored as is). */
export interface DecisionStep {
  readonly row: Bag;
  readonly match: StoredDecision | null;
  readonly op: Bag | null;
}

/**
 * The steps that bring the stored decisions to `ledgers`, one per target record: an insert for a
 * new record, an update of `superseded_*` for a stored record that was superseded since, else no
 * op. Never a delete, never a change of the immutable core (append-only, AC-8.1): the stored rows
 * no target record matched stay in the table and come back in `unmatched` (the tracker keeps them).
 */
export function decisionSteps(stored: ReadonlyArray<StoredDecision>, analysisId: string, studyId: string,
  ledgers: ReadonlyArray<ReviewLedger>): {steps: DecisionStep[]; unmatched: StoredDecision[]} {
  const pool = new Map<string, StoredDecision[]>();
  for (const s of stored) pool.set(s.key, [...(pool.get(s.key) ?? []), s]);
  const steps = decisionRowsOf(analysisId, studyId, ledgers).map((row): DecisionStep => {
    const match = pool.get(decisionKey(row))?.shift() ?? null;
    if (match === null) return {row, match, op: {op: 'insert', table: 'review_decision', values: row}};
    const same = match.supersededBy === row.superseded_by && match.supersededAt === row.superseded_at;
    return {row, match, op: same ? null : {op: 'update', table: 'review_decision', id: match.id,
      values: {superseded_by: row.superseded_by, superseded_at: row.superseded_at}}};
  });
  return {steps, unmatched: [...pool.values()].flat()};
}

// ─── Pipeline flags → pipeline_flag (decision 11) ───────────────────────────────

/** The columns of a `pipeline_flag` row the app owns. */
export const FLAG_COLUMNS: readonly string[] = ['profile_key', 'verdict', 'reasons', 'info', 'pipeline_only'];

export function flagRowsOf(analysisId: string, flags: ReadonlyArray<PipelineFlag>): Bag[] {
  return flags.map((f) => ({analysis_id: analysisId, profile_key: f.profileKey, verdict: f.verdict,
    reasons: [...f.reasons], info: [...f.info], pipeline_only: [...f.pipelineOnly]}));
}

/** FIELD-11 null semantics: rows alone cannot tell "app entry" from "pipeline entry with nothing
 *  flagged", so the origin decides — `[]` for a pipeline entry, `undefined` for any other. */
export function flagsOf(rows: readonly Bag[], originKind: string | null | undefined): PipelineFlag[] | undefined {
  if (rows.length === 0) return originKind === 'pipeline' ? [] : undefined;
  return [...rows].sort((x, y) => cmp(x.profile_key, y.profile_key)).map((r) => ({
    profileKey: r.profile_key, verdict: r.verdict, reasons: list(r.reasons), info: list(r.info),
    pipelineOnly: list(r.pipeline_only)}));
}

/** FIELD-11: the profiles with a pipeline-only reason (the count the lists show). */
const flagCountOf = (flags: ReadonlyArray<PipelineFlag> | undefined): number =>
  (flags ?? []).filter((f) => f.pipelineOnly.length > 0).length;

// ─── The whole entry ────────────────────────────────────────────────────────────

/** The {@link ENTRY_COLUMNS} of a domain entry. `resultsHash` = the snapshot hash its results were
 *  computed from (null without results). Both stamps round-trip `?? null` (architect B4): only a
 *  freshly computed summary sets them, in the entry it builds. */
export function entryCols(e: AnalysisEntry, resultsHash: string | null, recipeVersion: number): Bag {
  const recipe = recipeOf(e);
  const stored = e.source.kind === 'stored';
  return {
    label: e.label, engine_version: e.engineVersion ?? null, results_basis: e.resultsBasis ?? null,
    ...tierCols(e.summaryCache),
    results_snapshot_hash: e.summaryCache === null ? null : resultsHash,
    ...signatureCols(e.sourceSignature),
    source_kind: e.source.kind,
    source_file: stored ? `file://${e.source.ref}` : null, source_ref: stored ? null : e.source.ref,
    recipe, recipe_version: recipeVersion, ...recipeSearchCols(recipe),
    ...confirmedCols(e.confirmed ?? null, computeRecipeFingerprint(recipe)),
    ...stageCols(e.stageCache ?? null),
    ...originCols(e.origin ?? null),
    pipeline_flag_count: e.pipelineFlags !== undefined ? flagCountOf(e.pipelineFlags) : (e.pipelineFlagCount ?? 0),
  };
}

/** A domain entry from its stored rows. `children` that were not read stay as the entry's "absent"
 *  readings: no overrides / decisions read → `[]`; no flags read → the stored count only. */
export function entryOf(a: Bag, results: readonly Bag[], children: {
  readonly lz?: readonly Bag[]; readonly decisions?: readonly Bag[]; readonly flags?: readonly Bag[];
  readonly curves?: ReadonlyMap<string, CompoundCurve | null>;
} = {}): AnalysisEntry {
  const recipe = a.recipe as AnalysisRecipe;
  const stored = a.source_kind === 'stored' || a.source_kind === undefined || a.source_kind === null;
  const source = {kind: stored ? 'stored' as const : a.source_kind,
    ref: stored ? fromFileUrl(a.source_file) : (a.source_ref ?? '')};
  const flags = children.flags !== undefined ? flagsOf(children.flags, a.origin_kind) : undefined;
  return {
    ...recipe,
    id: a.id, label: a.label, source, sourceSignature: signatureOf(a),
    summaryCache: summaryCacheOf(a, results, children.curves),
    engineVersion: a.engine_version ?? null, resultsBasis: a.results_basis ?? null,
    resultsComputedAt: isoOf(a.results_computed_at),
    confirmed: confirmedOf(a, recipe),
    lzOverrides: lzOverridesOf(children.lz ?? []),
    reviewDecisions: ledgersOf(children.decisions ?? []),
    stageCache: stageOf(a),
    origin: originOf(a),
    ...(flags !== undefined ? {pipelineFlags: flags} : {}),
    pipelineFlagCount: a.pipeline_flag_count ?? 0,
  };
}

/** `file://<path>` → `<path>`. */
export const fromFileUrl = (url: string | undefined | null): string =>
  url === undefined || url === null ? '' : url.startsWith('file://') ? url.slice('file://'.length) : url;
