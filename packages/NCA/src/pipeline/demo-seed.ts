/**
 * `ncaSeedDemoStudies`: make sure the five Hub demo studies (nca-hub-demo-studies F6) are installed
 * and current — ordinary NCA studies, shared View to `All users`. Any user may call it from the
 * console (user decision 2026-11-04, superseding the admin-only seed): the first caller installs
 * the demo and owns it (owner and analyst on every pre-recorded override and decision). A later call
 * checks each demo study: up to date → nothing is written; missing → installed; out of date (another
 * sci-comp version, another results basis, or bundled data that no longer matches) → re-installed,
 * but only when the caller may delete it (its owner or an administrator); anyone else is told who can
 * refresh it, and nothing is written. Approach 1 was chosen for its simplicity (a rule-12 exception).
 *
 * Every analysis is built the way the app computes it on open (probe B1): `appRunOf` (the SAME
 * `NcaController` compute), then the planted λz overrides and review decisions recorded THROUGH the
 * controller (with the fingerprints of its own `reviewSnapshot`, so nothing auto-resets on open),
 * then `appStateFrom` (cache + stage from the post-decision state). Nothing is written as a literal.
 *
 * A demo study is recognised by its label AND its analysis labels; a same-label study with other
 * analyses is a collision, refused and never deleted. One study at a time: build every analysis (no
 * store), then create → append → share; a failure deletes that study (final), so a re-run never
 * meets a half-built study.
 */

import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {
  DEMO_STUDIES, DEMO_RULES_JSON, DEMO_THRESHOLDS, DEMO_SPARSE_VERDICTS,
  N_ONLY_DECISION, N_ONLY_MIN_ADJ_R2, DEMO_STRUCTURE_COL,
} from '../adapter/demo-studies';
import type {DemoAnalysisSpec, DemoDecision, DemoStudySpec} from '../adapter/demo-studies';
import {
  SPAN_RATIO_POOR, criterionReasons, judgementReasons, OTHER_REASON_CODE, getNcaRepository, numAt,
  passthroughCandidate, sharedFlagCodes, ACCEPT_ACCEPTED, ACCEPT_FLAGGED, ACCEPT_NA, isEntryStale,
  computeSnapshotHash,
} from '../adapter';
import type {
  AnalysisRecipe, NcaRepository, NewAnalysisInput, ReviewReason, ReviewUnit, StudyHeader,
  ReviewCategory, PassthroughColumnSpec,
} from '../adapter';
import {parseRules, type PipelineRules} from './pipeline-compute';
import {detectPipelineMapping, mappingGaps, sexOf, type PipelineMapping} from './pipeline-mapping';
import {appRunOf, appStateFrom, recipeFor, sourceFrame} from './pipeline-save';
import type {AppRun} from './pipeline-save';

/** What the seed needs from its surroundings — injected by the tests (a memory repository). */
export interface DemoSeedDeps {
  readonly repo: NcaRepository;
  /** Load one demo table (`DEMO_STUDIES_DIR/<file>`), as the app loads a bundled table. */
  readonly loadTable: (file: string) => Promise<DG.DataFrame>;
  /** The owner's display name of a study (named when a `replace` cannot delete it). */
  readonly ownerOf: (studyId: string) => Promise<string>;
  /** The group every demo study is shared with (View). */
  readonly shareGroupId: string;
  readonly now: () => string;
  /** Tests only: seed renamed copies (`zz-test …`) so a live run can never touch a real install. */
  readonly specs?: ReadonlyArray<DemoStudySpec>;
}

/** One analysis, built and ready to write. */
export interface BuiltDemoAnalysis {
  readonly spec: DemoAnalysisSpec;
  readonly input: NewAnalysisInput;
  /** The run it was built from (the tests read the controller's final state). */
  readonly run: AppRun;
  readonly recipe: AnalysisRecipe;
  /** The pipeline mapping and rules it was computed with. */
  readonly mapping: PipelineMapping;
  readonly rules: PipelineRules;
  readonly pending: number;
}

/** A demo build that no longer matches its planted design (the data or the engine changed). */
export class DemoDriftError extends Error {
  constructor(message: string) {super(`demo data drifted from the planted verdict: ${message}`);}
}

// ─── Building one analysis (no storage) ─────────────────────────────────────────

/** The pipeline mapping of a demo table — autodetected at `high` confidence with no gap (AC-2.1). */
export function demoMappingOf(df: DG.DataFrame, study: DemoStudySpec): PipelineMapping {
  const det = detectPipelineMapping(df);
  if (det.mapping === null || det.confidence !== 'high')
    throw new DemoDriftError(`${study.label}: mapping not detected with high confidence (${det.issues.join('; ')})`);
  const gaps = mappingGaps(det.mapping);
  if (gaps.length > 0) throw new DemoDriftError(`${study.label}: mapping leaves ${gaps.join(', ')} open`);
  return {
    ...det.mapping,
    ...(study.infusion !== undefined ? {infusionDuration: study.infusion.col} : {}),
    ...(study.sparse ? {design: 'sparse-destructive' as const} : {}),
  };
}

/** The rows of `table` an analysis holds. */
function analysisRows(table: DG.DataFrame, spec: DemoAnalysisSpec): number[] {
  const all = Array.from({length: table.rowCount}, (_, i) => i);
  if (spec.rows === undefined) return all;
  const c = table.col(spec.rows.col);
  if (c === null) throw new DemoDriftError(`no column '${spec.rows.col}'`);
  return all.filter((i) => !c.isNone(i) && c.getString(i) === spec.rows!.value);
}

/** The stored source frame of one demo analysis — ONE builder for the write and the freshness check,
 *  so the snapshot hash compared on a later call is the hash of exactly what was written. */
function demoSourceFrame(table: DG.DataFrame, study: DemoStudySpec, spec: DemoAnalysisSpec,
  m: PipelineMapping): DG.DataFrame {
  const src = sourceFrame(table, analysisRows(table, spec), m);
  src.name = study.frameName;
  // the passthrough semType rides the column; copy it from the detected table (subFrame keeps tags,
  // but say so explicitly — `isStructure` must come from the live detection, never a name match)
  for (const name of study.passthrough) {
    const from = table.col(name);
    const to = src.col(name);
    if (from !== null && to !== null && from.semType !== null) to.semType = from.semType;
  }
  return src;
}

/** The demo recipe: `recipeFor` (the pipeline's, replayed on open to the same fingerprint) with the
 *  pinned thresholds, the grouping, the annotations and — Study 5 only — the infusion fields (PD-5). */
function demoRecipe(m: PipelineMapping, pr: PipelineRules, src: DG.DataFrame, study: DemoStudySpec):
  AnalysisRecipe {
  const sex = sexOf(src, m);
  const base = recipeFor(m, DEMO_RULES_JSON, study.sparse, sex.ok ? sex : null);
  const recipe: AnalysisRecipe = {
    ...base,
    ncaRules: pr.rules,
    acceptanceRSquared: DEMO_THRESHOLDS.acceptanceRSquared,
    groupByCols: [...study.groupBy],
    ...(study.derived !== undefined ? {derived: study.derived} : {}),
    columns: study.infusion !== undefined ?
      {...base.columns, infusionSource: 'value', infusionDurationCol: null,
        infusionDurationScalar: study.infusion.hours} :
      base.columns,
  };
  const subjectCol = src.col(m.subject);
  const passthroughCols: PassthroughColumnSpec[] = study.passthrough.map((name) => {
    const col = src.col(name);
    if (col === null || !passthroughCandidate(col, recipe, subjectCol))
      throw new DemoDriftError(`${study.label}: '${name}' is not an annotation candidate`);
    // Captured ONCE from the live semType after detection (design row "passthrough").
    return {name, isStructure: col.semType === DG.SEMTYPE.MOLECULE};
  });
  return {...recipe, passthroughCols};
}

/** The pinned thresholds are what the recipe carries and what the code still uses (R1/F-5). */
function checkThresholds(pr: PipelineRules): void {
  const lz = pr.rules.lambdaZ;
  const ok = lz.minRSquared === DEMO_THRESHOLDS.minRSquared && lz.minPoints === DEMO_THRESHOLDS.minPoints &&
    pr.rules.extrapWarnPct === DEMO_THRESHOLDS.extrapWarnPct &&
    pr.acceptance.lambdaAcceptanceRSquared === DEMO_THRESHOLDS.acceptanceRSquared &&
    SPAN_RATIO_POOR === DEMO_THRESHOLDS.spanRatio;
  if (!ok) throw new DemoDriftError('the acceptance thresholds differ from the planted ones');
}

/** Rows of the controller's frame for one subject at the given times. */
function rowsAt(run: AppRun, m: PipelineMapping, subject: string, times: ReadonlyArray<number>): number[] {
  const s = run.work.col(m.subject)!;
  const t = run.work.col(m.time)!;
  const out: number[] = [];
  for (let i = 0; i < run.work.rowCount; i++) {
    if (s.getString(i) !== subject) continue;
    const ti = numAt(t, i);
    if (times.some((x) => Math.abs(x - ti) < 1e-9)) out.push(i);
  }
  if (out.length !== times.length) throw new DemoDriftError(`${subject}: override times ${times.join(', ')} not found`);
  return out;
}

/** Record one batch decision through the controller, refusing what the app would refuse. */
function recordBatch(run: AppRun, pr: PipelineRules, d: DemoDecision, where: string): void {
  const snap = run.controller.reviewSnapshot(pr.acceptance, {viable: run.viable, composites: run.composites});
  const units = d.units.map((k) => {
    const u = snap.units.find((x) => x.key === k);
    if (u === undefined) throw new DemoDriftError(`${where} ${d.ref}: no review unit '${k}'`);
    if (u.need !== 'review' || u.decision !== null || !u.allowedActions.includes(d.action))
      throw new DemoDriftError(`${where} ${d.ref}: '${k}' (${u.category}) does not allow ${d.action}`);
    return u;
  });
  const category = units[0].category as ReviewCategory;
  if (units.some((u) => u.category !== category)) throw new DemoDriftError(`${where} ${d.ref}: mixed categories`);
  const judge = judgementReasons(d.action, category);
  const reasons: ReviewReason[] = [
    ...(category === 'flagged' ? criterionReasons(sharedFlagCodes(units), pr.acceptance) : [])
      .map((r) => ({code: r.code, text: r.text})),
    ...d.codes.map((code) => {
      const def = judge.find((r) => r.code === code);
      if (def === undefined) throw new DemoDriftError(`${where} ${d.ref}: reason ${code} not offered for ${category}`);
      return {code: def.code, text: def.text};
    }),
    ...(d.other !== '' ? [{code: OTHER_REASON_CODE, text: d.other}] : []),
  ];
  for (const u of units) run.controller.decide(u.key, d.action, reasons, snap.fingerprints.get(u.key) ?? '');
}

/** The verdict word of a composite unit, as the grid shows it. */
function verdictOf(u: ReviewUnit): string {
  switch (u.category) {
  case 'accepted': return ACCEPT_ACCEPTED;
  case 'flagged': return ACCEPT_FLAGGED;
  case 'no-exposure': return 'No exposure';
  default: return ACCEPT_NA;
  }
}

/**
 * Build one demo analysis exactly as the app would hold it after the planted review: compute,
 * overrides, decisions, cache, stage. Throws {@link DemoDriftError} when anything planted no longer
 * holds (AC-3.2) — the seed never records a decision the app would refuse.
 */
export function buildDemoAnalysis(table: DG.DataFrame, study: DemoStudySpec, spec: DemoAnalysisSpec,
  now: string, opts: {review?: boolean} = {}): BuiltDemoAnalysis {
  const where = `${study.label} / ${spec.label}`;
  const m = demoMappingOf(table, study);
  const pr = parseRules(DEMO_RULES_JSON);
  checkThresholds(pr);
  const src = demoSourceFrame(table, study, spec, m);
  if (study.infusion !== undefined) {
    const c = src.col(study.infusion.col);
    const hours = study.infusion.hours;
    if (c === null || Array.from({length: c.length}, (_, i) => i).some((i) => c.isNone(i) || numAt(c, i) !== hours))
      throw new DemoDriftError(`${where}: '${study.infusion.col}' is not ${hours} h on every row`);
  }
  const recipe = demoRecipe(m, pr, src, study);
  const r = appRunOf(src, m, pr, recipe);
  if (!r.ok) throw new DemoDriftError(`${where}: ${r.reason}`);
  const run = r.run;
  for (const k of study.groupBy) {
    if (!run.groupCandidates.includes(k))
      throw new DemoDriftError(`${where}: '${k}' is not a Group-by candidate (${run.groupCandidates.join(', ')})`);
  }

  if (opts.review === false) {
    // The engine's own verdicts, before any planted review (the tests' "before" reading).
    const state = appStateFrom(run, pr, recipe, spec.label, now);
    return {spec, run, recipe, mapping: m, rules: pr, pending: state.stageCache.pending,
      input: {label: spec.label, df: src, recipe,
        summaryCache: state.summaryCache, origin: null, workflow: {confirmed: state.confirmed, lzOverrides: [],
          reviewDecisions: [], stageCache: state.stageCache}}};
  }
  for (const o of spec.overrides ?? []) {
    run.controller.overrideLambdaZPoints(o.subject,
      rowsAt(run, m, o.subject, o.times).map((rowIdx) => ({rowIdx, action: o.action})));
  }
  for (const d of spec.decisions ?? []) recordBatch(run, pr, d, where);

  if (study.sparse) {
    const pinned = DEMO_SPARSE_VERDICTS[spec.label];
    const snap = run.controller.reviewSnapshot(pr.acceptance, {viable: run.viable, composites: run.composites});
    const composites = snap.units.filter((u) => u.kind === 'composite');
    if (pinned === undefined || composites.length === 0 ||
      composites.some((u) => verdictOf(u) !== pinned.verdict)) {
      throw new DemoDriftError(`${where}: composite verdict ${composites.map(verdictOf).join(', ')} ≠ ` +
        `pinned ${pinned?.verdict ?? '(none)'}`);
    }
    if (pinned.decision !== null)
      recordBatch(run, pr, {...pinned.decision, units: composites.map((u) => u.key)}, where);
  }

  // The unplanned profiles flagged ONLY by N ≤ 3: one bulk Accept as is (user decision 2026-11-04).
  {
    const snap = run.controller.reviewSnapshot(pr.acceptance, {viable: run.viable, composites: run.composites});
    const keep = new Set(spec.leaveToReview ?? []);
    const nOnly = snap.units.filter((u) => u.need === 'review' && u.decision === null && !keep.has(u.key) &&
      u.category === 'flagged' && u.flagCodes.length === 1 && u.flagCodes[0] === 'NPTS_LE');
    const subj = run.outcome.results.col('SUBJECT')!;
    const adj = run.outcome.results.col('LAMBDAZ_ADJ_R2')!;
    for (const u of nOnly) {
      const i = Array.from({length: subj.length}, (_, j) => j).find((j) => subj.getString(j) === u.key);
      const v = i === undefined ? NaN : numAt(adj, i);
      if (!(v >= N_ONLY_MIN_ADJ_R2)) throw new DemoDriftError(`${where}: ${u.key} is N ≤ 3 with adj. R² ${v}`);
    }
    if (nOnly.length > 0)
      recordBatch(run, pr, {ref: 'N≤3', units: nOnly.map((u) => u.key), ...N_ONLY_DECISION}, where);
  }

  const state = appStateFrom(run, pr, recipe, spec.label, now);
  const final = run.controller.reviewSnapshot(pr.acceptance, {viable: run.viable, composites: run.composites});
  const pendingKeys = final.units.filter((u) => u.need === 'review' && u.decision === null).map((u) => u.key);
  if (state.stageCache.pending !== spec.expectedPending ||
    (spec.leaveToReview ?? []).some((k) => !pendingKeys.includes(k))) {
    throw new DemoDriftError(`${where}: To review = ${state.stageCache.pending} (${pendingKeys.join(', ')}), ` +
      `planned ${spec.expectedPending}`);
  }
  return {
    spec, run, recipe, mapping: m, rules: pr, pending: state.stageCache.pending,
    input: {
      label: spec.label, df: src, recipe, summaryCache: state.summaryCache, origin: null,
      workflow: {
        confirmed: state.confirmed, lzOverrides: run.controller.getLzOverrides(),
        reviewDecisions: run.controller.getReviewLedgers(), stageCache: state.stageCache,
      },
    },
  };
}

/**
 * Load a demo table and type its structure column. The column is DECLARED by the demo spec
 * (`DEMO_STRUCTURE_COL` — the data is ours, so it is SMILES), not detected: Chem's detector types a
 * SMILES column asynchronously, and the demo must install without the Chem package (user decision
 * 2026-11-04). With Chem installed the study grid draws the structures; without it the SMILES text
 * shows. `isStructure` is then captured from this semType once (design row "passthrough").
 */
export async function loadDemoTable(deps: Pick<DemoSeedDeps, 'loadTable'>, study: DemoStudySpec):
  Promise<DG.DataFrame> {
  const df = await deps.loadTable(study.file);
  const col = df.col(DEMO_STRUCTURE_COL);
  if (col === null) throw new DemoDriftError(`${study.label}: no '${DEMO_STRUCTURE_COL}' column`);
  col.semType = DG.SEMTYPE.MOLECULE;
  col.meta.units = DG.UNITS.Molecule.SMILES;
  df.name = study.frameName;
  return df;
}

/** Build every analysis of one study (no storage). */
export async function buildDemoStudy(deps: Pick<DemoSeedDeps, 'loadTable' | 'now'>,
  study: DemoStudySpec,
  opts: {review?: boolean} = {}): Promise<BuiltDemoAnalysis[]> {
  const table = await loadDemoTable(deps, study);
  return study.analyses.map((a) => buildDemoAnalysis(table, study, a, deps.now(), opts));
}

// ─── The protocol ───────────────────────────────────────────────────────────────

/** A visible study with a demo label: the demo itself (same analysis labels) or a collision. */
function demoMatch(visible: ReadonlyArray<StudyHeader>, spec: DemoStudySpec):
  {study: StudyHeader; isDemo: boolean} | null {
  const s = visible.find((x) => x.kind === 'study' && x.label === spec.label);
  if (s === undefined) return null;
  const want = spec.analyses.map((a) => a.label).sort().join('\u0001');
  const have = s.analyses.map((a) => a.label).sort().join('\u0001');
  return {study: s, isDemo: want === have};
}

/**
 * Why an installed demo study is out of date, or null when it is current: an analysis computed by
 * another sci-comp version or reporting basis (the ONE staleness predicate, `isEntryStale` — the
 * same that writes `stale (refreshes on open)`), or a stored snapshot that differs from the bundled
 * table (new demo data shipped with the package).
 */
export async function demoStaleness(repo: NcaRepository, studyId: string, spec: DemoStudySpec,
  table: DG.DataFrame): Promise<string | null> {
  const loaded = await repo.loadStudy(studyId);
  const m = demoMappingOf(table, spec);
  for (const a of spec.analyses) {
    const e = loaded.analyses.find((x) => x.label === a.label);
    if (e === undefined) return `analysis '${a.label}' is missing`;
    if (isEntryStale(e)) return `'${a.label}' was computed by another sci-comp version`;
    if (e.sourceSignature.hash !== computeSnapshotHash(demoSourceFrame(table, spec, a, m)))
      return `'${a.label}' holds older demo data`;
  }
  return null;
}

interface ReportRow {study: string; analyses: number; status: string; studyId: string}

/**
 * Make sure every demo study is installed and current (see the module header). Returns one row per
 * study — what was found and what was done; a planned re-install is listed before its delete runs.
 * Re-installing gives the study new ids, so copied links to it stop working.
 */
export async function seedDemoStudies(deps: DemoSeedDeps): Promise<DG.DataFrame> {
  const repo = deps.repo;
  const visible = await repo.listStudies();
  const rows: ReportRow[] = [];
  const todo: Array<{spec: DemoStudySpec; table: DG.DataFrame | null; replaces: StudyHeader | null}> = [];
  // The plan first: what exists, what is current, what is refused, what will be deleted (final).
  for (const spec of deps.specs ?? DEMO_STUDIES) {
    const match = demoMatch(visible, spec);
    if (match === null) {todo.push({spec, table: null, replaces: null}); continue;}
    if (!match.isDemo) {
      rows.push({study: spec.label, analyses: 0, studyId: match.study.id, status:
        `refused: a study named '${spec.label}' exists with other analyses — rename it and run again`});
      continue;
    }
    let table: DG.DataFrame;
    let stale: string | null;
    try {
      table = await loadDemoTable(deps, spec);
      stale = await demoStaleness(repo, match.study.id, spec, table);
    } catch (e) {
      rows.push({study: spec.label, analyses: match.study.analysisCount, studyId: match.study.id,
        status: `not checked: ${e instanceof Error ? e.message : String(e)}`});
      continue;
    }
    if (stale === null) {
      rows.push({study: spec.label, analyses: match.study.analysisCount, studyId: match.study.id,
        status: 'up to date'});
      continue;
    }
    if (!match.study.access.canDelete) {
      rows.push({study: spec.label, analyses: match.study.analysisCount, studyId: match.study.id, status:
        `out of date (${stale}) — only its owner, ${await deps.ownerOf(match.study.id)}, or an administrator ` +
        'can refresh it; nothing was changed'});
      continue;
    }
    rows.push({study: spec.label, analyses: match.study.analysisCount, studyId: match.study.id,
      status: `planned: re-install (${stale})`});
    todo.push({spec, table, replaces: match.study});
  }

  for (const {spec, replaces} of todo) {
    let built: BuiltDemoAnalysis[];
    try {
      built = await buildDemoStudy(deps, spec);
    } catch (e) {
      const msg = e instanceof Error ? e.message : String(e);
      rows.push({study: spec.label, analyses: replaces?.analysisCount ?? 0, studyId: replaces?.id ?? '',
        status: replaces !== null ? `not refreshed (kept): ${msg}` : `not installed: ${msg}`});
      continue;
    }
    if (replaces !== null) {
      await repo.deleteStudy(replaces.id);
      rows.push({study: spec.label, analyses: replaces.analysisCount, studyId: replaces.id, status: 'deleted'});
    }
    let studyId = '';
    try {
      studyId = (await repo.createStudy(spec.label)).id;
      for (const b of built) await repo.appendAnalysis(studyId, b.input);
      await repo.share(studyId, deps.shareGroupId);
      rows.push({study: spec.label, analyses: built.length, studyId, status: 'installed'});
    } catch (e) {
      // No half-built study: the delete is final, so a re-run starts clean (AC-6.6).
      if (studyId !== '') await repo.deleteStudy(studyId).catch(() => {});
      rows.push({study: spec.label, analyses: 0, studyId: '',
        status: `not installed: ${e instanceof Error ? e.message : String(e)}`});
    }
  }

  const df = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('Study', rows.map((r) => r.study)),
    DG.Column.fromInt32Array('Analyses', Int32Array.from(rows.map((r) => r.analyses))),
    DG.Column.fromStrings('Status', rows.map((r) => r.status)),
    DG.Column.fromStrings('StudyId', rows.map((r) => r.studyId)),
  ]);
  df.name = 'NCA demo studies';
  return df;
}

/** The production dependencies: the session repository, the bundled tables. */
export function productionSeedDeps(loadTable: (file: string) => Promise<DG.DataFrame>): DemoSeedDeps {
  return {
    repo: getNcaRepository(),
    loadTable,
    ownerOf: async (id) => {
      try {
        const [e] = await grok.dapi.getEntities([id]);
        // eslint-disable-next-line @typescript-eslint/no-explicit-any
        return ((e as any)?.author?.friendlyName as string | undefined) ?? 'another user';
      } catch {return 'another user';}
    },
    shareGroupId: DG.Group.defaultGroupsIds['All users'],
    now: () => new Date().toISOString(),
  };
}
