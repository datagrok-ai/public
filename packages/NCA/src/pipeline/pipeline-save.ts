/**
 * `ncaSaveToStudy`: save a pipeline run into NCA Studio so a scientist can open a
 * flagged analysis in the app (λz review, exclusions, Confirm).
 *
 * Layout:
 * - One NCA study per delivered study (the first `groupBy` key, e.g. `Study`) and owner. (owner,
 *   study key) is the study's EMS business key, so re-runs and concurrent runs find the same study,
 *   and the pipeline only ever writes its OWN study (never another user's study with the same label).
 * - One analysis per pipeline partition (compound, period, auto-split column); a sparse partition
 *   is one analysis PER ROUTE, labelled `<partition> · <route>`, because the app's sparse panel
 *   computes with one route per analysis. The app keys profiles on the subject alone (BUG-07), so a
 *   partition never shares an analysis with another one.
 * - The source CSV is the delivered data only. The pipeline's verdict and reasons live on the entry
 *   (`pipelineFlags`), the run's provenance in `origin` (schema v16).
 * - The entry is saved confirmed, with the summary cache and the review stage the APP computes on
 *   open (same `NcaController` compute, the app's own `reviewSnapshot`, the app's acceptance rule),
 *   so the hub and the study grid show what the scientist will see.
 *
 * Re-running never overwrites a scientist's work:
 * - every run writes its data to NEW snapshots (a snapshot an analysis already points at is never
 *   rewritten);
 * - an analysis is updated in place only when it is the pipeline's own entry for that partition and
 *   is untouched since the run that wrote it (`pipelineUntouched`); otherwise the run is saved beside;
 * - the plan is written in ONE transaction under the versions of the base it was planned on: a base
 *   that moved meanwhile (a scientist's save, another run) re-plans (≤ 3), never overwrites.
 */

import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {
  INSTALLED_ENGINE_VERSION, STANDALONE_ROUTE_SEGMENT,
  computeRecipeFingerprint, computeRecipeHash, RECIPE_HASH_VERSION, buildProfileIndices, appendGroupingColumns,
  buildSummaryCache, newId, uniqueLabel, NCA_APP_ROOT, makeAnalysisEntry, recipeOf, AuditLog, deriveStage,
  denseProfilesViable, computeSparseForRecipe, sparseCompositeUnits, sparseUnitKey, numAt,
  DERIVED_ROUTE_CODE_COL, evaluateAcceptance, ensureBoolBlq, sparseCompositeKey, sparseDoseGroups,
  getNcaRepository, KeepsChangingError, reportableResults, writeReviewColumns, groupByCandidates,
} from '../adapter';
import type {
  AcceptanceThresholds,
  AnalysisEntry, AnalysisRecipe, NcaStudy, PipelineFlag, PipelineOrigin, StageCache,
  SummaryCache, ConfirmedRecipe, ObservationsCtx, SamplingType, SparseUnits, NcaRepository, SourcePointer,
  CompositeUnitInput,
} from '../adapter';
import {NcaController} from '../controller/nca-controller';
import type {NcaRunOutcome} from '../controller/nca-controller';
import {
  COMPOSITE_SUBJECT, VERDICT_BLOCKED, VERDICT_REVIEW, groupRows, parseRules, prepareAppRun,
  rowRoutes, strAt, subFrame, REASON_NCA_FAILED, REASON_LZ_NOT_ESTIMABLE, type PipelineRules,
} from './pipeline-compute';
import {mappedColumns, parseMapping, sexOf, unitsOf, type PipelineMapping} from './pipeline-mapping';

// ─── The analysis units of a run ───────────────────────────────────────────────

interface Partition {
  study: string;
  label: string;
  keys: Map<string, string>;
  /** Observation rows of the partition. */
  rows: number[];
  /** Result rows of the partition (all, blocked included). */
  resultRows: number[];
  sparse: boolean;
}

/** One analysis to save: a dense partition, or one route of a sparse partition. */
interface SaveUnit {
  study: string;
  /** The partition label (no route). */
  base: string;
  /** The analysis identity — `base`, or `base · <route>` for a sparse analysis (AC-4.0). */
  partition: string;
  sparse: boolean;
  /** Observation rows of `table`. */
  rows: number[];
  /** Non-blocked result rows. */
  resultRows: number[];
}

/** One saved / skipped analysis in the function's output. */
interface SavedRow {
  study: string; analysis: string; verdict: string; profiles: number; review: number;
  studyId: string; analysisId: string; url: string; status: string;
}

/** Pipeline partitions (study × the other key columns). */
function partitionsOf(table: DG.DataFrame, results: DG.DataFrame, keyNames: string[]): Partition[] {
  const subj = results.col('SUBJECT')!;
  const verdict = results.col('VERDICT')!;
  const keyOf = (df: DG.DataFrame, i: number): string => keyNames.map((k) => strAt(df.col(k), i)).join('\u0001');
  const resByKey = groupRows(Array.from({length: results.rowCount}, (_, i) => i), (i) => keyOf(results, i));
  const obsByKey = groupRows(Array.from({length: table.rowCount}, (_, i) => i),
    (i) => keyNames.map((k) => table.col(k) === null ? '' : strAt(table.col(k), i)).join('\u0001'));
  const parts: Partition[] = [];
  for (const [key, resultRows] of resByKey) {
    const values = key.split('\u0001');
    const live = resultRows.filter((r) => verdict.get(r) !== VERDICT_BLOCKED);
    parts.push({
      study: values[0] || '(no study)',
      label: values.slice(1).filter((v) => v !== '').join(' · ') || 'All',
      keys: new Map(keyNames.map((k, i) => [k, values[i]])),
      // An auto-split partition's key includes the split column, which the observations carry too.
      rows: obsByKey.get(key) ?? [],
      resultRows,
      sparse: live.length > 0 && live.every((r) => strAt(subj, r) === COMPOSITE_SUBJECT),
    });
  }
  return parts;
}

/**
 * The analyses of the run: a dense partition is one unit; a sparse partition is one unit per route
 * (F6). A unit with no non-blocked result is reported, not saved.
 */
function saveUnitsOf(table: DG.DataFrame, results: DG.DataFrame, parts: Partition[], m: PipelineMapping): {
  units: SaveUnit[]; blocked: Array<{study: string; label: string; reason: string}>;
} {
  const verdict = results.col('VERDICT')!;
  const reasons = results.col('REASONS');
  const route = results.col('ROUTE')!;
  const units: SaveUnit[] = [];
  const blocked: Array<{study: string; label: string; reason: string}> = [];
  for (const p of parts) {
    const firstReason = (rows: number[]): string => strAt(reasons, rows[0]) || 'every profile is blocked';
    if (!p.sparse) {
      const live = p.resultRows.filter((r) => verdict.get(r) !== VERDICT_BLOCKED);
      if (live.length === 0 || p.rows.length === 0) {
        blocked.push({study: p.study, label: p.label, reason: firstReason(p.resultRows)});
        continue;
      }
      units.push({study: p.study, base: p.label, partition: p.label, sparse: false, rows: p.rows, resultRows: live});
      continue;
    }
    const routes = rowRoutes(subFrame(table, p.rows), m);
    if (!routes.ok) {
      blocked.push({study: p.study, label: p.label, reason: routes.reason});
      continue;
    }
    const rowsByRoute = groupRows(p.rows.map((_, i) => i), (i) => routes.routes[i]);
    const resultsByRoute = groupRows(p.resultRows, (r) => strAt(route, r));
    for (const [rt, idx] of rowsByRoute) {
      const label = `${p.label} · ${rt}`;
      const res = resultsByRoute.get(rt) ?? [];
      const live = res.filter((r) => verdict.get(r) !== VERDICT_BLOCKED);
      if (live.length === 0) {
        blocked.push({study: p.study, label, reason: res.length > 0 ? firstReason(res) : 'no composite computed'});
        continue;
      }
      units.push({
        study: p.study, base: p.label, partition: label, sparse: true,
        rows: idx.map((i) => p.rows[i]), resultRows: live,
      });
    }
  }
  return {units, blocked};
}

// ─── The recipe and the source CSV ─────────────────────────────────────────────

/** The sampling-type declaration of a sparse recipe from the mapping's `design` (AC-6.2): `auto`
 *  stays `(auto)`, so the app derives the topology from the data exactly as `ncaCompute` did. */
function samplingTypeOf(m: PipelineMapping): SamplingType {
  switch (m.design) {
  case 'sparse-destructive': return 'destructive';
  case 'sparse-batch': return 'batch';
  default: return '(auto)';
  }
}

/** The recipe the app would hold for this mapping (every role a column, nothing scalar). */
export function recipeFor(m: PipelineMapping, rulesJson: string, sparse: boolean,
  sex: {col: string; map: Record<string, string>} | null = null): AnalysisRecipe {
  const pr = parseRules(rulesJson);
  const u = unitsOf(m);
  const col = (name: string | undefined): string | null => name ?? null;
  // An unmapped optional role is `'value'` with a 0 scalar, never `'none'`: the pane's replay maps
  // `'none'` to `'value'` (`applyRecipe`), so a `'none'` here would fingerprint differently from the
  // reopened Setup, and the reopen would land on Setup with no compute (AC-5.7, measured in the live
  // test). The app computes an unmapped role exactly this way, so nothing numerical changes.
  return {
    columns: {
      subject: m.subject, time: m.time, concentration: m.concentration,
      doseSource: m.dose !== undefined ? 'column' : 'value', doseCol: col(m.dose), doseScalar: 0,
      // A mapping without a route column has a defaultRoute (`mappingGaps`); 'PO' is never reached.
      // A sparse analysis holds ONE route (F6), so the app's "first route cell" read is exact.
      routeSource: m.route !== undefined ? 'column' : 'value', routeScalar: m.defaultRoute ?? 'PO',
      routeCol: col(m.route),
      // The app's "no infusion" is a 0 scalar, not 'none' (it is part of the recipe fingerprint).
      infusionSource: 'value', infusionDurationCol: null, infusionDurationScalar: 0,
      bodyWeightSource: m.bodyWeight !== undefined ? 'column' : 'value', bodyWeightCol: col(m.bodyWeight),
      bodyWeightScalar: 0,
      blqFlagCol: col(m.blq),
      lloqSource: m.lloq !== undefined ? 'column' : 'value', lloqCol: col(m.lloq), lloqScalar: 0,
      sexCol: sex?.col ?? null, sexValueMap: sex?.map ?? {},
    },
    modes: {
      nonmem: false, confirmedRoute: null,
      sparse: {
        on: sparse,
        nominalTimeCol: sparse ? (m.nominalTime ?? m.time) : null,
        animalIdCol: sparse ? m.subject : null,
        samplingType: sparse ? samplingTypeOf(m) : '(auto)',
      },
    },
    ncaRules: pr.rules,
    units: {time: u.time, concentration: u.concentration, dose: u.dose, bodyWeight: u.bodyWeight},
    groupByCols: [],
    passthroughCols: [],
    acceptanceRSquared: pr.acceptance.lambdaAcceptanceRSquared,
  };
}

/** The analysis's source frame: its rows and the delivered columns. An all-empty column is dropped
 *  only when the mapping does not name it — a mapped LLOQ / BLQ / groupBy column stays even when
 *  empty, or the recipe would name a column the CSV lacks (review 6b). */
export function sourceFrame(table: DG.DataFrame, rows: ReadonlyArray<number>, m: PipelineMapping): DG.DataFrame {
  const sub = subFrame(table, rows);
  const keep = mappedColumns(m);
  for (const c of sub.columns.toList()) {
    if (c.name.startsWith('~') || (!keep.has(c.name) && c.stats.missingValueCount === sub.rowCount))
      sub.columns.remove(c.name);
  }
  return sub;
}

// ─── What the app computes on open ─────────────────────────────────────────────

/** The cache, stage and confirmed recipe of one analysis — computed the way the app computes them. */
export interface AppState {
  readonly summaryCache: SummaryCache;
  readonly stageCache: StageCache;
  readonly confirmed: ConfirmedRecipe;
}

/** One analysis computed the way the app computes it on open (see {@link appRunOf}). */
export interface AppRun {
  /** The SAME `NcaController` the app drives — overrides and decisions go through it. */
  readonly controller: NcaController;
  /** Its completed outcome: `results` and `perProfile` are patched in place by an override. */
  readonly outcome: Extract<NcaRunOutcome, {status: 'completed'}>;
  /** The frame the controller computes on (a clone of the source, with the derived route codes). */
  readonly work: DG.DataFrame;
  readonly observations: ObservationsCtx;
  /** The app's UC-04 composites for this recipe (empty for a dense recipe). */
  readonly composites: ReadonlyArray<CompositeUnitInput & {fingerprint: string}>;
  /** `denseProfilesViable` of the work frame — the per-profile review units' gate. */
  readonly viable: boolean;
  /** The Report-stage Group-by candidates the app's ribbon offers for this analysis. */
  readonly groupCandidates: ReadonlyArray<string>;
}

/**
 * The app's compute of one analysis (review 5A): the SAME `NcaController` compute the app runs on
 * open — no `acceptanceIgnore`, outlier, terminal-rise, comment or engine routing (those are
 * pipeline policy, in `pipelineFlags`). A sparse analysis's controller runs the dense per-animal
 * compute, exactly as the app does; its review units are the app's composites (computed with the
 * app's sparse inputs for this recipe). Split from the cache + stage ({@link appStateFrom}) so a
 * caller can record λz overrides and review decisions through the controller in between
 * (nca-hub-demo-studies PD-3, probe B1) — the pipeline records none.
 */
export function appRunOf(src: DG.DataFrame, m: PipelineMapping, pr: PipelineRules, recipe: AnalysisRecipe):
  {ok: true; run: AppRun} | {ok: false; reason: string} {
  const work = src.clone();
  const blqName = m.blq !== undefined && work.col(m.blq) !== null ? ensureBoolBlq(work, m.blq) : null;
  const prep = prepareAppRun(work, m, blqName);
  if (!prep.ok) return prep;
  const controller = new NcaController(work, prep.contract, prep.scalars, pr.rules,
    {log: new AuditLog(), ...(prep.routeWarnings !== undefined ? {routeWarnings: prep.routeWarnings} : {})});
  const outcome = controller.computeAll();
  if (outcome.status !== 'completed')
    return {ok: false, reason: `validation: ${outcome.errors.map((e) => e.message).join('; ')}`};
  const profiles = buildProfileIndices(work, prep.contract);
  const carried = appendGroupingColumns(outcome.results, work, profiles, prep.contract);
  const observations: ObservationsCtx = {
    df: work, timeName: m.nominalTime ?? m.time, concName: m.concentration, doseName: m.dose ?? null,
    routeCodeName: work.col(DERIVED_ROUTE_CODE_COL) !== null ? DERIVED_ROUTE_CODE_COL : null,
  };
  // The app's sparse inputs for this recipe, on the work frame `prepareAppRun` gave the route codes
  // (AD-5: no third route reader) — the composites the app computes on open.
  const composites = recipe.modes.sparse.on ? sparseCompositeUnits(computeSparseForRecipe(work, {
    subjectName: m.subject, timeName: m.time, concName: m.concentration,
    doseColName: m.dose ?? null, bwColName: m.bodyWeight ?? null,
    units: recipe.units as SparseUnits,
    route: {
      codesCol: work.col(DERIVED_ROUTE_CODE_COL) !== null ? DERIVED_ROUTE_CODE_COL : null,
      scalar: recipe.columns.routeScalar,
    },
    lloq: {colName: recipe.columns.lloqCol, scalar: recipe.columns.lloqScalar},
    blqColName: recipe.columns.blqFlagCol,
    nominalTimeName: recipe.modes.sparse.nominalTimeCol,
    animalIdName: recipe.modes.sparse.animalIdCol,
    declaredTopology: recipe.modes.sparse.samplingType === '(auto)' ? null : recipe.modes.sparse.samplingType,
    acceptance: pr.acceptance, bootstrapIterations: 20,
  }).result) : [];
  return {ok: true, run: {
    controller, outcome, work, observations, composites, viable: denseProfilesViable(work, m.subject),
    groupCandidates: groupByCandidates(carried, outcome.results, work, DERIVED_ROUTE_CODE_COL),
  }};
}

/**
 * The summary cache, stage and confirmed recipe of a run, from its CURRENT state — after any λz
 * override and review decision recorded through `run.controller` (AC-3.2a): the review snapshot is
 * taken here, the `Review` column is written from that same snapshot, and the cache is built from
 * the reportable (Exclude-λz-blanked) results with the patched `perProfile`, exactly as the app's
 * save builds it. With no override and no decision this is today's pipeline output (AC-3.1).
 */
export function appStateFrom(run: AppRun, pr: PipelineRules, recipe: AnalysisRecipe, label: string,
  now: string): AppState {
  const {controller, outcome} = run;
  const snap = controller.reviewSnapshot(pr.acceptance, {viable: run.viable, composites: run.composites});
  writeReviewColumns(outcome.results, snap.units);
  const reportable = reportableResults(outcome.results, controller.excludedKeys()).df;
  const summaryCache = buildSummaryCache(reportable, outcome.perProfile, recipe, label, run.observations);
  const fingerprint = computeRecipeFingerprint(recipe);
  const stage = deriveStage({
    confirmedFingerprint: fingerprint, currentFingerprint: fingerprint, hasResults: true,
    units: snap.units.map((u) => ({need: u.need, decided: u.decision !== null})),
  });
  return {
    summaryCache,
    stageCache: {stage: stage.stage, pending: stage.pending, notApplicable: stage.notApplicable},
    // The app confirms the Setup recipe, whose passthrough selection is always empty (LBC-2).
    confirmed: {recipe: {...recipe, passthroughCols: []}, fingerprint, at: now},
  };
}

/** {@link appRunOf} then {@link appStateFrom}: the app's view of one analysis with no review work. */
function appStateOf(src: DG.DataFrame, m: PipelineMapping, pr: PipelineRules, recipe: AnalysisRecipe,
  label: string, now: string): {ok: true; state: AppState} | {ok: false; reason: string} {
  const r = appRunOf(src, m, pr, recipe);
  return r.ok ? {ok: true, state: appStateFrom(r.run, pr, recipe, label, now)} : r;
}

// ─── Pipeline flags (Q1) ───────────────────────────────────────────────────────

/** Two reads of one dose level. Not `===`: the results `DOSE` passes through a platform column
 *  whose storage can round to single precision (0.3 → 0.30000001…), while the raw source cell is
 *  read exact (measured in the pure tests). Dose levels differ by orders of magnitude more than
 *  this relative tolerance, so two levels can never be merged by it. */
function sameDose(a: number, b: number): boolean {
  return Math.abs(a - b) <= 1e-6 * Math.max(1, Math.abs(b));
}

const parseList = (text: string): string[] => {
  if (text === '') return [];
  const v = JSON.parse(text) as unknown;
  return Array.isArray(v) ? v.map((x) => String(x)) : [];
};

/**
 * The pipeline's per-profile flags of one analysis (AC-8.2), parsed from the `_JSON` columns
 * (never by splitting the text). A dense profile is keyed by `SUBJECT`. A composite is keyed by the
 * app's unit key, `sparse:DOSE=<raw dose>` (`sparse:all` with no dose column), rebuilt from the
 * analysis's own source rows with the adapter's grouping rule (`sparseDoseGroups`) and key builder
 * (`sparseCompositeKey`) — the ones the app keys its composites with. A saved sparse analysis holds
 * one route (F6), so the key never carries a route (architect R2 N-2). Since BUG-19 the app and the
 * pipeline group a dose level identically ("10" and "10.0" are one composite, keyed by the smallest
 * spelling), so a numeric composite matches exactly one key. A composite matching no source row
 * throws: an orphan flag is never written.
 *
 * `pipelineOnly` = the reasons the app's own check does not raise (approved design 2026-10-15):
 * every reason minus the app's λz acceptance reasons for the row (the SAME `evaluateAcceptance`
 * on the row's own fit statistics, with no `acceptanceIgnore`) and the two "not estimated" reasons
 * the app shows as its Not-estimated category.
 */
export function pipelineFlagsOf(results: DG.DataFrame, resultRows: ReadonlyArray<number>,
  src: DG.DataFrame, m: PipelineMapping, acceptance: AcceptanceThresholds): PipelineFlag[] {
  const subj = results.col('SUBJECT')!;
  const verdict = results.col('VERDICT')!;
  const reasons = results.col('REASONS_JSON')!;
  const info = results.col('INFO_FLAGS_JSON')!;
  const doseRes = results.col('DOSE');
  const rawDose = m.dose !== undefined ? src.col(m.dose) : null;
  const num = (name: string, r: number): number => {
    const c = results.col(name);
    return c === null || c.isNone(r) ? NaN : numAt(c, r);
  };
  const doseGroups = rawDose !== null ?
    sparseDoseGroups(rawDose, Array.from({length: rawDose.length}, (_, i) => i)) : [];
  const out: PipelineFlag[] = [];
  for (const r of resultRows) {
    const all = parseList(strAt(reasons, r));
    const appRaised = new Set<string>([REASON_NCA_FAILED, REASON_LZ_NOT_ESTIMABLE, ...evaluateAcceptance({
      spanRatio: num('SPAN_RATIO', r), adjR2: num('LAMBDAZ_ADJ_R2', r),
      nPoints: num('LAMBDAZ_NPOINTS', r), pctExtrap: num('PCT_EXTRAP', r),
    }, acceptance).reasons]);
    const flag = {
      verdict: (strAt(verdict, r) === VERDICT_REVIEW ? 'Review' : 'Accepted') as PipelineFlag['verdict'],
      reasons: all, info: parseList(strAt(info, r)), pipelineOnly: all.filter((x) => !appRaised.has(x)),
    };
    if (strAt(subj, r) !== COMPOSITE_SUBJECT) {
      out.push({profileKey: strAt(subj, r), ...flag});
      continue;
    }
    if (rawDose === null) {
      out.push({profileKey: sparseUnitKey(sparseCompositeKey(null, null)), ...flag});
      continue;
    }
    const d = doseRes === null || doseRes.isNone(r) ? NaN : numAt(doseRes, r);
    // A text dose label has no number in `DOSE` (NaN): it matches every text group, as before.
    const keys = doseGroups.filter((g) => (Number.isNaN(d) ? Number.isNaN(g.value) : sameDose(g.value, d)))
      .map((g) => g.key);
    if (keys.length === 0)
      throw new Error(`ncaSaveToStudy: composite at dose ${d} matches no source row — cannot key its flag`);
    for (const k of keys) out.push({profileKey: sparseUnitKey(sparseCompositeKey(null, k)), ...flag});
  }
  return out;
}

// ─── The pure decision (F3 / F4) ───────────────────────────────────────────────

/**
 * Review 6d: may a re-run update `e` in place? Only when it is the pipeline's entry for this
 * partition and nothing has touched it since the run that wrote it — no λz override, no review
 * ledger (even one reset to `To review`: its history is the scientist's record), and the WHOLE
 * stored recipe (Setup, Group-by, annotations) is what that run wrote AND is the confirmed one (a
 * saved, unconfirmed Setup edit is work too). A hash from another `RECIPE_HASH_VERSION` cannot be
 * compared, so it counts as touched (fail-safe: the run is saved beside).
 */
export function pipelineUntouched(e: AnalysisEntry, partition: string): boolean {
  const o = e.origin;
  if (o === null || o === undefined || o.kind !== 'pipeline' || o.partition !== partition) return false;
  if ((e.lzOverrides ?? []).length > 0 || (e.reviewDecisions ?? []).length > 0) return false;
  if (o.recipeHashVersion !== RECIPE_HASH_VERSION) return false;
  const recipe = recipeOf(e);
  return computeRecipeHash(recipe) === o.recipeHash && e.confirmed?.fingerprint === computeRecipeFingerprint(recipe);
}

/** One analysis a run wants to save (its entry already carries this run's source, cache and origin). */
export interface PlannedAnalysis {
  readonly base: string;
  readonly partition: string;
  readonly sparse: boolean;
  readonly entry: AnalysisEntry;
}

/** What the run decided for one planned analysis. */
export interface PlannedOutcome {
  readonly partition: string;
  readonly analysisId: string;
  readonly label: string;
  readonly status: string;
}

export interface StudyPlan {
  readonly study: NcaStudy;
  readonly outcomes: ReadonlyArray<PlannedOutcome>;
  /** Source refs of entries updated in place — removed once nothing references them (AC-3.6). */
  readonly replacedRefs: ReadonlyArray<string>;
  /** Pipeline entries of another sampling design, left in place (AC-4.2a). */
  readonly superseded: ReadonlyArray<{analysisId: string; label: string; by: string}>;
}

/**
 * The pure save decision (F3 / F4), over the stable read `base` (`undefined` = new study). Every
 * entry of `base` is kept; a planned analysis replaces its partition's untouched pipeline entry in
 * place (same id, the STORED label — a scientist's rename is kept), or is appended beside with a
 * unique label. A pipeline entry of the same partition base but another sampling design (dense ↔
 * sparse) is left in place and reported; the new entry records it in `origin.supersedes`.
 */
export function planStudyUpdate(base: NcaStudy | undefined, planned: ReadonlyArray<PlannedAnalysis>, ctx: {
  studyLabel: string; pipelineKey: {owner: string; studyKey: string}; now: string;
}): StudyPlan {
  const analyses: AnalysisEntry[] = [...(base?.analyses ?? [])];
  const taken = new Set(analyses.map((a) => a.label));
  const outcomes: PlannedOutcome[] = [];
  const replacedRefs: string[] = [];
  const superseded: Array<{analysisId: string; label: string; by: string}> = [];
  const partitions = new Set(planned.map((p) => p.partition));
  const pipelineOf = (a: AnalysisEntry): PipelineOrigin | null => (a.origin?.kind === 'pipeline' ? a.origin : null);

  for (const p of planned) {
    const candidates = analyses.filter((a) => pipelineOf(a)?.partition === p.partition);
    const untouched = candidates.filter((a) => pipelineUntouched(a, p.partition))
      .sort((x, y) => pipelineOf(y)!.runAt.localeCompare(pipelineOf(x)!.runAt));
    // A design flip: this partition base was saved under the other sampling design by an earlier
    // run (dense `A` ↔ sparse `A · <route>`), and no planned analysis of this run updates it.
    // The candidate must be of the OTHER design: a label prefix alone would also match an unrelated
    // dense partition such as `Cmpd · Day1` (another groupBy) and report it superseded.
    const flipped = analyses.find((a) => {
      const o = pipelineOf(a);
      if (o === null || partitions.has(o.partition) || superseded.some((s) => s.analysisId === a.id)) return false;
      if ((a.modes?.sparse?.on === true) === p.sparse) return false;
      return p.sparse ? o.partition === p.base : o.partition.startsWith(`${p.base} · `);
    });
    const origin = pipelineOf(p.entry);
    const entry: AnalysisEntry = flipped !== undefined && origin !== null ?
      {...p.entry, origin: {...origin, supersedes: flipped.id}} : p.entry;
    let outcome: PlannedOutcome;
    if (untouched.length > 0) {
      const old = untouched[0];
      analyses[analyses.findIndex((a) => a.id === old.id)] = {...entry, id: old.id, label: old.label};
      if (old.source.ref !== entry.source.ref) replacedRefs.push(old.source.ref);
      outcome = {partition: p.partition, analysisId: old.id, label: old.label, status: 'updated'};
    } else {
      const label = uniqueLabel(p.partition, taken);
      taken.add(label);
      analyses.push({...entry, label});
      outcome = {
        partition: p.partition, analysisId: entry.id, label,
        status: candidates.length > 0 ?
          `saved beside '${candidates[candidates.length - 1].label}' (has review work)` : 'created',
      };
    }
    outcomes.push(outcome);
    if (flipped !== undefined) superseded.push({analysisId: flipped.id, label: flipped.label, by: outcome.label});
  }
  const study: NcaStudy = base !== undefined ?
    {...base, analyses, updatedAt: ctx.now} :
    {
      id: '', label: ctx.studyLabel, kind: 'study', createdAt: ctx.now, updatedAt: ctx.now,
      access: {canEdit: true, canDelete: true, canShare: true}, analyses, pipelineKey: ctx.pipelineKey,
    };
  return {study, outcomes, replacedRefs, superseded};
}

// ─── Owner / share (Q2) ────────────────────────────────────────────────────────

/** What `share` asks for: the groups to grant View on a NEW study, and whether an Edit list was given
 *  (decision 10: View only until the platform honours row Edit grants — SD-16 reports it). */
interface ShareRequest {
  readonly viewGroupIds: string[];
  readonly editAsked: boolean;
}

/** `share` (AC-3.2): `{"view": [...], "edit": [...]}` of group names or ids. An unknown group or a
 *  malformed value is a configuration error and throws. The runner owns the study already (the
 *  platform grants its author everything), so `share` only adds other groups. */
async function shareRequestOf(shareJson: string): Promise<ShareRequest> {
  if (shareJson.trim() === '') return {viewGroupIds: [], editAsked: false};
  let o: unknown;
  try {
    o = JSON.parse(shareJson);
  } catch (e) {
    throw new Error(`share: invalid JSON — ${e instanceof Error ? e.message : String(e)}`);
  }
  if (o === null || typeof o !== 'object' || Array.isArray(o)) throw new Error('share: expected a JSON object');
  const rec = o as Record<string, unknown>;
  for (const k of Object.keys(rec))
    if (k !== 'view' && k !== 'edit') throw new Error(`share: unknown key '${k}' (expected 'view', 'edit')`);
  const names = (k: 'view' | 'edit'): string[] => {
    const v = rec[k];
    if (v === undefined) return [];
    if (!Array.isArray(v) || v.some((x) => typeof x !== 'string'))
      throw new Error(`share: '${k}' must be an array of group names`);
    return v as string[];
  };
  const groups = await grok.dapi.groups.list();
  const idOf = (n: string): string => {
    const g = groups.find((x) => x.id === n || x.name === n || x.friendlyName === n);
    if (g === undefined) throw new Error(`share: unknown group '${n}'`);
    return g.id;
  };
  // An Edit group is granted View (decision 10) and the run says so (SD-16) — never a silent downgrade.
  const edit = names('edit').map(idOf);
  return {viewGroupIds: [...new Set([...names('view').map(idOf), ...edit])], editAsked: edit.length > 0};
}

/** The runner's id, login and the ids of the groups it belongs to. */
export async function runner(): Promise<{id: string; login: string; groupNames: Set<string>}> {
  const u = await grok.dapi.users.current();
  // PR-11: the groups the caller belongs to (transitively), by id, name and friendly name.
  const groupNames = new Set<string>();
  try {
    for (const g of await grok.dapi.groups.currentUserGroups())
      for (const n of [g.id, g.name, g.friendlyName]) if (n) groupNames.add(n);
  } catch {/* without the list only the caller's own login or id can own a study */}
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  return {id: u.id, login: ((u as any).login as string | undefined) ?? u.name, groupNames};
}

/**
 * PR-11: `owner` names who a pipeline study belongs to. It must be the caller (its login or id) or a
 * group the caller belongs to — otherwise a run could write into, or probe for, another owner's
 * studies. Client-side by nature (the server cannot bind a business key to its caller, PR-19);
 * existing rows stay protected by row security either way.
 */
function checkOwner(owner: string, me: {id: string; login: string; groupNames: Set<string>}): void {
  if (owner === me.login || owner === me.id || me.groupNames.has(owner)) return;
  throw new Error(`Owner '${owner}' is not you or one of your groups`);
}

// ─── Serialisation per study within one client (AC-3.7a) ──────────────────────

const studyChains = new Map<string, Promise<unknown>>();

/** Run `job` after every earlier job of the same (owner, study key) in this client has finished. */
function serialized<T>(key: string, job: () => Promise<T>): Promise<T> {
  const prev = studyChains.get(key) ?? Promise.resolve();
  const next = prev.then(job, job);
  studyChains.set(key, next.catch(() => undefined));
  return next;
}

// ─── The save ──────────────────────────────────────────────────────────────────

/** Attempts to plan on a base that keeps moving before the run reports a conflict (AC-3.4d). */
const MAX_PLAN_ATTEMPTS = 3;

/** One prepared analysis, before its snapshot is written. */
interface PreparedUnit {
  readonly unit: SaveUnit;
  /** The partition's source frame — stored as its d42 snapshot, exactly the frame the run computed. */
  readonly df: DG.DataFrame;
  readonly recipe: AnalysisRecipe;
  readonly state: AppState;
  readonly flags: PipelineFlag[];
  readonly verdicts: string[];
}

/**
 * Save the non-blocked partitions of a pipeline run as NCA Studio analyses.
 * `table` is the observations table given to `ncaCompute`; `results` its output. `share` (JSON
 * `{"view": [...], "edit": [...]}` of group names) grants View on a NEW study only (decision 10); `owner`
 * (default: the runner's id; must be the runner or one of its groups — PR-11) is part of the study
 * identity — set it stably for scheduled runs so a change of service account does not duplicate
 * studies. Returns one row per analysis (or skipped partition) with the app URL and what was done.
 */
export async function saveToStudies(table: DG.DataFrame, results: DG.DataFrame,
  share: string = '', owner: string = '',
  /** Tests only: another repository (an in-memory one, for the save-protocol tests). */
  opts: {repo?: NcaRepository} = {}): Promise<DG.DataFrame> {
  const mappingTag = results.getTag('nca.mapping');
  if (mappingTag === null || mappingTag === undefined || mappingTag === '')
    throw new Error('ncaSaveToStudy: results carry no mapping — pass the output of ncaCompute');
  const m = parseMapping(mappingTag);
  const rulesJson = results.getTag('nca.rules') ?? '';
  const pr = parseRules(rulesJson);
  const keyNames = (JSON.parse(results.getTag('nca.keyColumns') ?? '[]') as string[])
    .filter((k) => results.col(k) !== null);
  const {units, blocked} = saveUnitsOf(table, results, partitionsOf(table, results, keyNames), m);

  const repo = opts.repo ?? getNcaRepository();
  const me = await runner();
  const ownerKey = owner.trim() !== '' ? owner.trim() : me.id;
  checkOwner(ownerKey, me);
  const shareReq = await shareRequestOf(share);
  const runId = newId();
  const now = new Date().toISOString();
  const out: SavedRow[] = [];

  for (const [studyKey, idx] of groupRows(units.map((_, i) => i), (i) => units[i].study)) {
    const studyUnits = idx.map((i) => units[i]);
    out.push(...await serialized(`${ownerKey}\u0001${studyKey}`, () => saveStudy({
      repo, table, results, m, pr, rulesJson, studyKey, units: studyUnits,
      ownerKey, shareReq, user: me.login, runId, now,
    })));
  }
  for (const b of blocked) {
    out.push({study: b.study, analysis: b.label, verdict: VERDICT_BLOCKED, profiles: 0, review: 0,
      studyId: '', analysisId: '', url: '', status: `not saved: ${b.reason}`});
  }

  const df = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('Study', out.map((r) => r.study)),
    DG.Column.fromStrings('Analysis', out.map((r) => r.analysis)),
    DG.Column.fromStrings('Verdict', out.map((r) => r.verdict)),
    DG.Column.fromInt32Array('Profiles', Int32Array.from(out.map((r) => r.profiles))),
    DG.Column.fromInt32Array('Review', Int32Array.from(out.map((r) => r.review))),
    DG.Column.fromStrings('Url', out.map((r) => r.url)),
    DG.Column.fromStrings('Status', out.map((r) => r.status)),
    DG.Column.fromStrings('StudyId', out.map((r) => r.studyId)),
    DG.Column.fromStrings('AnalysisId', out.map((r) => r.analysisId)),
    DG.Column.fromStrings('RunId', out.map((r) => (r.studyId !== '' ? runId : ''))),
  ]);
  df.name = 'Saved to NCA Studio';
  return df;
}

/** Save one study's analyses (see the module header for the protocol). */
async function saveStudy(a: {
  repo: NcaRepository; table: DG.DataFrame; results: DG.DataFrame; m: PipelineMapping;
  pr: PipelineRules; rulesJson: string; studyKey: string; units: SaveUnit[];
  ownerKey: string; shareReq: ShareRequest; user: string; runId: string; now: string;
}): Promise<SavedRow[]> {
  const {repo, m, studyKey, now} = a;
  const rows: SavedRow[] = [];
  const notSaved = (u: SaveUnit, reason: string): SavedRow => ({
    study: studyKey, analysis: u.partition, verdict: VERDICT_BLOCKED, profiles: 0, review: 0,
    studyId: '', analysisId: '', url: '', status: `not saved: ${reason}`,
  });

  // 1. Compute everything per analysis (no storage).
  const prepared: PreparedUnit[] = [];
  for (const u of a.units) {
    const src = sourceFrame(a.table, u.rows, m);
    const sex = sexOf(src, m);
    const recipe = recipeFor(m, a.rulesJson, u.sparse, sex.ok ? sex : null);
    const app = appStateOf(src, m, a.pr, recipe, u.partition, now);
    if (!app.ok) {
      rows.push(notSaved(u, app.reason));
      continue;
    }
    prepared.push({
      unit: u, df: src, recipe, state: app.state,
      flags: pipelineFlagsOf(a.results, u.resultRows, src, m, a.pr.acceptance),
      verdicts: u.resultRows.map((r) => strAt(a.results.col('VERDICT'), r)),
    });
  }
  if (prepared.length === 0) return rows;

  // 2. Our own study by its business key (owner, study key), or a new one with a label no visible
  //    study uses. A study of this key the caller cannot see is not found here; creating it then
  //    hits the key and is refused below (no oracle beyond the caller's own owner namespace, PR-11).
  let base: NcaStudy | undefined;
  try {
    base = (await repo.pipelineBase(a.ownerKey, studyKey)).study;
  } catch (e) {
    return [...rows, ...prepared.map((p) => notSaved(p.unit, e instanceof Error ? e.message : String(e)))];
  }
  let studyLabel = base?.label ?? studyKey;
  let createdNote = '';
  if (base === undefined) {
    const visible = (await repo.listStudies()).filter((s) => s.kind === 'study');
    studyLabel = uniqueLabel(studyKey, new Set([...visible.map((s) => s.label), STANDALONE_ROUTE_SEGMENT]));
    if (studyLabel !== studyKey) createdNote = ` (a study named '${studyKey}' already exists)`;
  }

  // 3. Store every analysis's data as a NEW snapshot (never an existing one).
  const written: SourcePointer[] = [];
  const planned: PlannedAnalysis[] = [];
  try {
    for (const p of prepared) {
      const snap = await repo.writeSnapshot(p.df, now);
      written.push(snap.source);
      const draft = makeAnalysisEntry({
        id: newId(), label: p.unit.partition, source: snap.source, sourceSignature: snap.signature,
        recipe: p.recipe, summaryCache: p.state.summaryCache,
        workflow: {confirmed: p.state.confirmed, lzOverrides: [], reviewDecisions: [], stageCache: p.state.stageCache},
        origin: null, pipelineFlags: p.flags,
      });
      const origin: PipelineOrigin = {
        kind: 'pipeline', runId: a.runId, runAt: now, user: a.user, engineVersion: INSTALLED_ENGINE_VERSION,
        mapping: {...m} as unknown as Record<string, unknown>, rules: a.rulesJson, partition: p.unit.partition,
        recipeHash: computeRecipeHash(recipeOf(draft)), recipeHashVersion: RECIPE_HASH_VERSION,
      };
      planned.push({base: p.unit.base, partition: p.unit.partition, sparse: p.unit.sparse, entry: {...draft, origin}});
    }
  } catch (e) {
    for (const s of written) await repo.deleteSnapshot(s);
    return [...rows, ...prepared.map((p) => notSaved(p.unit, e instanceof Error ? e.message : String(e)))];
  }

  // 4. Plan on the base read and write the plan in one transaction under its versions; a base that
  //    moved meanwhile re-reads and re-plans (≤ 3) — never an array planned on an old base.
  let plan: StudyPlan | null = null;
  let studyId = '';
  // S5: a study of this key the caller cannot write (403), or one the create hit but a re-read still
  // cannot see (another owner's, invisible by row security), is ONE refusal — no oracle beyond it.
  const notYours = `Study '${studyKey}' exists and is not yours`;
  let keyTaken = false;
  try {
    for (let attempt = 0; attempt < MAX_PLAN_ATTEMPTS && plan === null; attempt++) {
      if (attempt > 0) {
        base = (await repo.pipelineBase(a.ownerKey, studyKey)).study;
        if (base === undefined && keyTaken) throw new Error(notYours);
      }
      const candidate = planStudyUpdate(base, planned, {studyLabel, pipelineKey: {owner: a.ownerKey, studyKey}, now});
      try {
        studyId = await repo.savePipelinePlan({base, planned: candidate.study, owner: a.ownerKey, studyKey});
        plan = candidate;
      } catch (e) {
        if (e instanceof DG.DomainForbiddenError) throw new Error(notYours);
        if (!(e instanceof KeepsChangingError)) throw e;
        keyTaken = base === undefined; // the create met the key: our own concurrent run, or not ours
      }
    }
  } catch (e) {
    // The transaction was refused: nothing points at this run's snapshots.
    for (const s of written) await repo.deleteSnapshot(s);
    return [...rows, ...prepared.map((p) => notSaved(p.unit, e instanceof Error ? e.message : String(e)))];
  }
  if (plan === null) {
    for (const s of written) await repo.deleteSnapshot(s);
    return [...rows, ...prepared.map((p) =>
      notSaved(p.unit, `study '${studyLabel}' kept changing during the save — re-run`))];
  }

  // 5. The transaction committed: the snapshots the run replaced are no longer pointed at (AC-3.6).
  for (const ref of plan.replacedRefs) await repo.deleteSnapshot({kind: 'stored', ref});
  // A NEW study gets the requested View groups (decision 10; an Edit list is granted View — SD-16).
  let shareNote = '';
  if (base === undefined) {
    for (const g of a.shareReq.viewGroupIds) await repo.share(studyId, g);
    if (a.shareReq.editAsked) shareNote = ' · share: edit not available yet (View granted)';
  }

  for (let i = 0; i < prepared.length; i++) {
    const p = prepared[i];
    const o = plan.outcomes[i];
    const nReview = p.verdicts.filter((v) => v === VERDICT_REVIEW).length;
    rows.push({
      study: studyKey, analysis: o.label, verdict: nReview > 0 ? VERDICT_REVIEW : 'Accepted',
      profiles: p.verdicts.length, review: nReview, studyId, analysisId: o.analysisId,
      url: `${NCA_APP_ROOT}/${encodeURIComponent(plan.study.label)}/${encodeURIComponent(o.label)}`,
      status: (o.status === 'created' && createdNote !== '' ?
        `created study '${plan.study.label}'${createdNote}` : o.status) + shareNote,
    });
  }
  for (const s of plan.superseded) {
    rows.push({
      study: studyKey, analysis: s.label, verdict: '', profiles: 0, review: 0, studyId, analysisId: s.analysisId,
      url: `${NCA_APP_ROOT}/${encodeURIComponent(plan.study.label)}/${encodeURIComponent(s.label)}`,
      status: `superseded by '${s.by}' (sampling design changed) — left in place`,
    });
  }
  return rows;
}
