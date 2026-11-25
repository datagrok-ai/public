/**
 * The production {@link NcaRepository}: the platform EMS schema `nca` 0.2.0 (`databases/nca/schema.json`,
 * typed client `src/generated/db.ts`) and d42 snapshots in the analysis row's `file` column.
 *
 * Row ↔ domain mapping lives in `ems-nca-rows.ts`; the domain types (`NcaStudy`, `AnalysisEntry`)
 * stay. An analysis is an aggregate: its `analysis` row plus `analysis_result`, `lz_override`,
 * `review_decision` and `pipeline_flag` rows, and (for a pipeline entry) a `pipeline_run` row.
 *
 * Every read and write records, per analysis id, what this session last saw of that aggregate (the
 * {@link Tracked} state): the row `version`, its `source_file`, the content base (recipe, confirmed
 * fingerprint and time, λz overrides, review decisions — PR-1) and, per column and per child row,
 * the stored values. A write diffs its target against that state and sends ONLY what differs, in one
 * transaction guarded by the `analysis` version (decision 14) — so a reopen that recomputes the same
 * numbers writes nothing at all (decision 2). Requests per screen: the lists 2 (studies + analyses),
 * a study 3 in parallel (study + analyses + results) and one curve read after the first render, an
 * open 5 in parallel + the snapshot, a Save 0 reads.
 */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {ncaDb} from '../generated/db';
import type {AnalysisColumn, AnalysisResultColumn} from '../generated/db';
import {INSTALLED_ENGINE_VERSION} from './engine-version';
import {RESULTS_BASIS} from './unit-conversion';
import {recipeOf, DERIVED_BASIS} from './study-schema';
import type {
  NcaStudy, AnalysisEntry, AnalysisRecipe, SourcePointer, SourceSignature, SummaryCache, CompoundCurve,
  LzOverride, ReviewLedger, PipelineFlag, PipelineOrigin,
} from './study-schema';
import {computeRecipeFingerprint} from './fingerprint';
import {
  NCA_STUDIES_CHANGED_EVENT, newId, uniqueLabel, computeSnapshotHash, snapshotSignature, snapshotDiff, snapshotFrame,
  SourceReplacedError, ContentChangedError, KeepsChangingError, SnapshotFidelityError,
  STANDALONE_ROUTE_SEGMENT,
} from './nca-repository';
import type {
  NcaRepository, StudyHeader, AnalysisHeader, LoadedStudy, LoadedAnalysis,
  NewAnalysisInput, SaveAnalysisInput, SummaryWrite, MoveTarget, MoveResult, PipelineBase, PipelineWrite,
  WrittenSnapshot, StudyAccess, WorkflowState,
} from './nca-repository';
import {
  canonical, canonicalColumns, canonicalRow, isoOf, fromFileUrl, entryCols, entryOf, ENTRY_COLUMNS,
  RESULT_ENTRY_COLUMNS, RESULT_COLUMNS, RESULT_COLUMNS_WITHOUT_CURVE, resultRowsOf, tierCols, stageOf, stageCols,
  confirmedCols, recipeSearchCols, originCols, pipelineRunCols, lzRowsOf, lzKey, lzOverridesOf, normLz, LZ_COLUMNS,
  ledgersOf, normLedgers, decisionSteps, decisionKey, storedDecisionOf, flagRowsOf, FLAG_COLUMNS,
} from './ems-nca-rows';
import type {Bag, StoredDecision} from './ems-nca-rows';

/** Where snapshots live (the platform convention `<schema>/<table>/<uuid>/<file>`). */
const SNAPSHOT_ROOT = 'System:DomainFiles/nca/analysis';
/** One page of a list read: the server's JSON cap is 10k, its default 100 (E3). */
const PAGE = 1000;
/** Version races retried before a write gives up (PR-1, PR-14). */
const MAX_ATTEMPTS = 3;
/** The recipe layout this build writes; a row written by a newer one is not opened (SD-15). */
const RECIPE_VERSION = 1;

/** The light analysis columns of the lists (SD-17 / BD-1). */
const LIGHT_ANALYSIS_COLUMNS: AnalysisColumn[] = [
  'study_id', 'label', 'engine_version', 'results_basis', 'has_cache', 'stage', 'stage_pending',
  'stage_not_applicable', 'pipeline_flag_count', 'moved_from_label', 'updated_at',
];

/** Every analysis column the study overview reads (SD-21): all of them — the ledgers and flags
 *  are child tables now, not read here. */
const STUDY_ANALYSIS_COLUMNS: AnalysisColumn[] = ['study_id', 'moved_from_label', 'updated_at',
  'results_computed_at', 'pipeline_run', ...ENTRY_COLUMNS];

/** Ops a Save never sends: the label (rename's), the origin and pipeline provenance (the
 *  pipeline's / a move's) are not the view's to change. */
const NOT_SAVED: ReadonlySet<string> = new Set(['label', 'origin_kind', 'origin_label', 'origin_container_id',
  'pipeline_partition', 'pipeline_recipe_hash', 'pipeline_recipe_hash_version', 'pipeline_supersedes',
  'pipeline_flag_count']);
/** The columns a Save that keeps the snapshot leaves as they are. */
const SNAPSHOT_COLUMNS: ReadonlySet<string> = new Set(['source_kind', 'source_file', 'source_ref', 'source_sig_kind',
  'source_mtime', 'source_hash', 'source_computed_at', 'row_count']);

// ─── Small helpers ─────────────────────────────────────────────────────────────

/** One `=` condition (a list value is `= ANY`; a dotted property filters through a ref). */
const eq = (property: string, value: string | number | boolean | null | string[]): never =>
  ({property, operator: '=', value}) as never;

const nowIso = (): string => new Date().toISOString();

/** A version race: keyed on the typed error (every `/domains` call, transactions included, rethrows
 *  through `domainCall`), never on message text — a label such as "Run 409" must not read as one. */
export function isConflict(e: unknown): boolean {
  return e instanceof DG.DomainVersionConflictError || (e instanceof DG.DomainError && e.code === 'version-conflict');
}

/** A definitive refusal (the write did not land): a server 4xx answer. A timeout (408), a throttle
 *  (429), a transport error (no typed status) or anything else is NOT — the write may have landed,
 *  so its snapshot must stay (SD-5). Exported with {@link isConflict} for the classification test. */
export function isDefinitiveRefusal(e: unknown): boolean {
  if (!(e instanceof DG.DomainError)) return false;
  return e.status >= 400 && e.status < 500 && e.status !== 408 && e.status !== 429;
}

/** Every page of a query. */
async function allPages<T>(read: (offset: number) => Promise<T[]>): Promise<T[]> {
  const out: T[] = [];
  for (let offset = 0; ; offset += PAGE) {
    const page = await read(offset);
    out.push(...page);
    if (page.length < PAGE) return out;
  }
}

/** Rows grouped by a key column. */
function groupBy(rows: readonly Bag[], key: string): Map<string, Bag[]> {
  const out = new Map<string, Bag[]>();
  for (const r of rows) out.set(r[key], [...(out.get(r[key]) ?? []), r]);
  return out;
}

/** The content a Save must not overwrite when another session changed it (PR-1): the recipe, the
 *  confirmed fingerprint and time (the snapshot is derived from the recipe they hash), the λz
 *  overrides and the review decisions — in the normalised form a stored set reads back as. */
function baseOf(recipe: unknown, confirmedFingerprint: string | null, confirmedAt: unknown,
  lz: ReadonlyArray<LzOverride>, review: ReadonlyArray<ReviewLedger>): string {
  return canonical({recipe: recipe ?? null, cf: confirmedFingerprint, ca: isoOf(confirmedAt), lz: normLz(lz),
    review: normLedgers(review)});
}

const baseOfEntry = (e: AnalysisEntry): string => baseOf(recipeOf(e), e.confirmed?.fingerprint ?? null,
  e.confirmed?.at ?? null, e.lzOverrides ?? [], e.reviewDecisions ?? []);

// ─── Row ↔ domain (the parts not in ems-nca-rows) ─────────────────────────────────

function accessOf(r: Bag): StudyAccess {
  return {canEdit: r['~can_edit'] === true, canDelete: r['~can_delete'] === true, canShare: r['~can_share'] === true};
}

function headerOf(r: Bag): AnalysisHeader {
  return {
    id: r.id, studyId: r.study_id, label: r.label, engineVersion: r.engine_version ?? null,
    resultsBasis: r.results_basis ?? null, hasCache: r.has_cache === true, stageCache: stageOf(r),
    pipelineFlagCount: r.pipeline_flag_count ?? 0,
    movedFromLabel: r.moved_from_label ?? null, updatedAt: isoOf(r.updated_at ?? r.updated_on) ?? '',
  };
}

/** The domain entry a write input describes (source = the written snapshot). */
function entryFromInput(id: string, input: NewAnalysisInput & {label: string}, snap: WrittenSnapshot): AnalysisEntry {
  const w: WorkflowState = input.workflow;
  return {
    ...input.recipe,
    id, label: input.label, source: snap.source, sourceSignature: snap.signature,
    summaryCache: input.summaryCache,
    engineVersion: input.summaryCache === null ? null : INSTALLED_ENGINE_VERSION,
    resultsBasis: input.summaryCache === null ? null : RESULTS_BASIS,
    derivedBasis: input.summaryCache === null ? null : DERIVED_BASIS,
    confirmed: w.confirmed, lzOverrides: w.lzOverrides, reviewDecisions: w.reviewDecisions, stageCache: w.stageCache,
    origin: input.origin ?? null,
    ...(input.pipelineFlags !== undefined ? {pipelineFlags: input.pipelineFlags} : {}),
  };
}

// ─── The tracked aggregate and the write planner ──────────────────────────────────

/** A stored child row as a later write diffs against it (`canon` null = content not read). */
interface StoredRow {id: string; canon: string | null}
interface StoredLz {id: string; key: string; canon: string}

/** The child rows of an analysis as last read or written. `results` / `flags` null = not read. */
interface ChildState {
  results: StoredRow[] | null;
  lz: StoredLz[];
  decisions: StoredDecision[];
  flags: StoredRow[] | null;
}

/** What this session last read or wrote of one analysis — the CAS token and diff base of its next write. */
interface Tracked extends ChildState {
  version: number;
  studyId: string;
  label: string;
  sourceFile: string;
  hash: string | null;
  recipeFingerprint: string;
  base: string;
  /** Canonical stored value per {@link ENTRY_COLUMNS} column. */
  cols: Record<string, string>;
}

function childStateOf(results: readonly Bag[] | null, lz: readonly Bag[], decisions: readonly Bag[],
  flags: readonly Bag[] | null): ChildState {
  return {
    results: results === null ? null : [...results].sort((x, y) => (x.ord ?? 0) - (y.ord ?? 0))
      .map((r) => ({id: r.id, canon: canonicalRow(r, RESULT_COLUMNS)})),
    lz: lz.map((r) => ({id: r.id, key: lzKey(r.profile_key, r.time_key, r.ordinal ?? 0),
      canon: canonicalRow(r, LZ_COLUMNS)})),
    decisions: decisions.map(storedDecisionOf),
    flags: flags === null ? null : [...flags].sort((x, y) => (x.profile_key < y.profile_key ? -1 : 1))
      .map((r) => ({id: r.id, canon: canonicalRow(r, FLAG_COLUMNS)})),
  };
}

function trackRead(a: Bag, results: readonly Bag[] | null, lz: readonly Bag[], decisions: readonly Bag[],
  flags: readonly Bag[] | null): Tracked {
  return {
    version: a.version, studyId: a.study_id, label: a.label, sourceFile: fromFileUrl(a.source_file),
    hash: a.source_hash ?? null, recipeFingerprint: computeRecipeFingerprint(a.recipe as AnalysisRecipe),
    base: baseOf(a.recipe, a.confirmed_fingerprint ?? null, a.confirmed_at, lzOverridesOf(lz), ledgersOf(decisions)),
    cols: canonicalColumns(a, ENTRY_COLUMNS),
    ...childStateOf(results, lz, decisions, flags),
  };
}

/** The columns of `cols` whose value differs from the tracked one (all of them without a tracker). */
function changedCols(cols: Bag, t: Tracked | null): Bag {
  const keys = Object.keys(cols);
  if (t === null) return {...cols};
  const now = canonicalColumns(cols, keys);
  return Object.fromEntries(keys.filter((k) => now[k] !== t.cols[k]).map((k) => [k, cols[k]]));
}

/** What the child rows of a write are built from (an entry, or a background summary write). */
interface ChildInput {
  readonly summaryCache: SummaryCache | null;
  readonly lzOverrides?: ReadonlyArray<LzOverride>;
  readonly reviewDecisions?: ReadonlyArray<ReviewLedger>;
  readonly pipelineFlags?: ReadonlyArray<PipelineFlag>;
  readonly engineVersion?: string | null;
  readonly resultsBasis?: string | null;
}

/** The child ops of one write, and how to read the child state they leave once committed. */
interface ChildPlan {
  readonly ops: Bag[];
  readonly resultsChanged: boolean;
  /** The child state after the commit, given the transaction results of {@link ops} (same order). */
  finish(res: readonly Bag[], studyId: string): ChildState;
}

/**
 * The child-row ops that bring `prev` to `x`, for the parts `what` names (the others stay as they
 * are): result rows overwritten in place by `ord` only where they differ, extras deleted, new ones
 * inserted; λz rows by their business key (update, insert, delete); review decisions append-only
 * ({@link decisionSteps}: inserts + supersede updates, never a delete); pipeline flags replaced as a
 * set when they differ. `ref` / `studyRef` are the analysis and study ids, or `$ref`s of the same
 * transaction.
 */
function planChildren(ref: string, studyRef: string, x: ChildInput, prev: ChildState,
  what: {results?: boolean; lz?: boolean; decisions?: boolean; flags?: boolean}): ChildPlan {
  const ops: Bag[] = [];
  const after: Array<(r: Bag) => void> = [];
  const push = (op: Bag, f: (r: Bag) => void = () => {}): void => {ops.push(op); after.push(f);};

  let results: Array<{id: string; row: Bag}> | null = null;
  if (what.results) {
    if (prev.results === null) throw new Error('planChildren: the result rows were not read');
    const old = prev.results;
    const rows = x.summaryCache === null ? [] : resultRowsOf(ref, studyRef, x.summaryCache,
      {engineVersion: x.engineVersion ?? null, resultsBasis: x.resultsBasis ?? null});
    results = rows.map((row, i) => {
      const slot = {id: old[i]?.id ?? '', row};
      if (old[i] === undefined) push({op: 'insert', table: 'analysis_result', values: row}, (r) => {slot.id = r.id;});
      else if (old[i].canon !== canonicalRow(row, RESULT_COLUMNS)) {
        const {analysis_id: _a, ...values} = row;
        push({op: 'update', table: 'analysis_result', id: old[i].id, values});
      }
      return slot;
    });
    for (const o of old.slice(rows.length)) push({op: 'delete', table: 'analysis_result', id: o.id});
  }
  const resultsChanged = ops.length > 0;

  let lz: StoredLz[] | null = null;
  if (what.lz) {
    const old = new Map(prev.lz.map((o) => [o.key, o]));
    const want = new Map(lzRowsOf(ref, x.lzOverrides ?? [])
      .map((row) => [lzKey(row.profile_key, row.time_key, row.ordinal), row]));
    lz = [...want].map(([key, row]) => {
      const canon = canonicalRow(row, LZ_COLUMNS);
      const o = old.get(key);
      old.delete(key);
      const slot: StoredLz = {id: o?.id ?? '', key, canon};
      if (o === undefined) push({op: 'insert', table: 'lz_override', values: row}, (r) => {slot.id = r.id;});
      else if (o.canon !== canon) {
        push({op: 'update', table: 'lz_override', id: o.id, values: {conc_key: row.conc_key, action: row.action,
          set_by: row.set_by, set_by_label: row.set_by_label, set_at: row.set_at}});
      }
      return slot;
    });
    for (const o of old.values()) push({op: 'delete', table: 'lz_override', id: o.id});
  }

  let decisions: StoredDecision[] | null = null;
  if (what.decisions) {
    const {steps, unmatched} = decisionSteps(prev.decisions, ref, studyRef, x.reviewDecisions ?? []);
    decisions = [...unmatched];
    for (const s of steps) {
      const slot: StoredDecision = {id: s.match?.id ?? '', key: decisionKey(s.row),
        supersededBy: s.row.superseded_by ?? null, supersededAt: s.row.superseded_at ?? null};
      decisions.push(slot);
      if (s.op !== null) push(s.op, (r) => {if (s.match === null) (slot as {id: string}).id = r.id;});
    }
  }

  let flags: Array<{id: string; canon: string}> | null = null;
  if (what.flags && x.pipelineFlags !== undefined) {
    const rows = flagRowsOf(ref, x.pipelineFlags)
      .sort((a, b) => (a.profile_key < b.profile_key ? -1 : 1));
    const canons = rows.map((r) => canonicalRow(r, FLAG_COLUMNS));
    const old = prev.flags ?? [];
    if (canonical(canons) === canonical(old.map((o) => o.canon))) flags = old.map((o) => ({id: o.id, canon: o.canon!}));
    else {
      for (const o of old) push({op: 'delete', table: 'pipeline_flag', id: o.id});
      flags = rows.map((row, i) => {
        const slot = {id: '', canon: canons[i]};
        push({op: 'insert', table: 'pipeline_flag', values: row}, (r) => {slot.id = r.id;});
        return slot;
      });
    }
  }

  return {
    ops, resultsChanged,
    finish(res, studyId) {
      after.forEach((f, i) => f(res[i] ?? {}));
      return {
        // The canonical form with the REAL study id: a `$ref` must not read as a later difference.
        results: results === null ? prev.results :
          results.map((s) => ({id: s.id, canon: canonicalRow({...s.row, study_id: studyId}, RESULT_COLUMNS)})),
        lz: lz ?? prev.lz,
        decisions: decisions ?? prev.decisions,
        flags: flags ?? prev.flags,
      };
    },
  };
}

/** The tracker of an aggregate after a committed write: the stored columns merged with the written
 *  ones, the child state the plan left. */
function trackWrite(prev: Tracked | null, a: {version: number; studyId: string; label: string; sourceFile: string;
  hash: string | null; recipeFingerprint: string; base: string; written: Bag; children: ChildState}): Tracked {
  const written = canonicalColumns(a.written, Object.keys(a.written).filter((k) =>
    (ENTRY_COLUMNS as readonly string[]).includes(k)));
  return {
    version: a.version, studyId: a.studyId, label: a.label, sourceFile: a.sourceFile, hash: a.hash,
    recipeFingerprint: a.recipeFingerprint, base: a.base, cols: {...(prev?.cols ?? {}), ...written}, ...a.children,
  };
}

const EMPTY_CHILDREN: ChildState = {results: [], lz: [], decisions: [], flags: []};

// ─── The repository ───────────────────────────────────────────────────────────

export class EmsNcaRepository implements NcaRepository {
  private readonly tracked = new Map<string, Tracked>();
  /** What each pipeline base read saw, keyed by the base object it returned: a plan written on that
   *  base is checked against THESE versions and diffed against THESE rows, so a plan on a stale base
   *  always re-plans (SD-20). Keyed by the object, not the study id — two overlapping runs in one
   *  session each keep their own read, and a later read cannot refresh an earlier base into passing. */
  private readonly pipelineReads = new WeakMap<NcaStudy, {study: number; analyses: Map<string, Tracked>}>();

  // ── reads ──────────────────────────────────────────────────────────────────

  /** The list read in flight: the Hub's lists and the Browse tree refresh on the same change event,
   *  and share one read instead of each issuing their own (BD-12). Never a stale result — only
   *  callers that arrive while the read is running share it. */
  private listing: Promise<StudyHeader[]> | null = null;

  listStudies(): Promise<StudyHeader[]> {
    if (this.listing === null) {
      const read = this.readStudies().finally(() => {if (this.listing === read) this.listing = null;});
      this.listing = read;
    }
    return this.listing;
  }

  /** After a structural write: a list read started before it must not be shared any more, then
   *  the Hub lists and the tree refresh (NCA_STUDIES_CHANGED_EVENT). */
  private changed(): void {
    this.listing = null;
    grok.events.fireCustomEvent(NCA_STUDIES_CHANGED_EVENT, {});
  }

  private async readStudies(): Promise<StudyHeader[]> {
    const [studies, analyses] = await Promise.all([
      allPages((offset) => ncaDb.studies.query({columns: ['label', 'kind', 'updated_at'], withAccess: true,
        sort: 'label', limit: PAGE, offset}) as Promise<Bag[]>),
      allPages((offset) => ncaDb.analyses.query({columns: LIGHT_ANALYSIS_COLUMNS, sort: 'label',
        limit: PAGE, offset}) as Promise<Bag[]>),
    ]);
    const byStudy = new Map<string, AnalysisHeader[]>();
    for (const a of analyses) {
      const h = headerOf(a);
      byStudy.set(h.studyId, [...(byStudy.get(h.studyId) ?? []), h]);
    }
    return studies.map((s) => {
      const list = byStudy.get(s.id) ?? [];
      return {id: s.id, label: s.label, kind: s.kind, updatedAt: isoOf(s.updated_at ?? s.updated_on) ?? '',
        analysisCount: list.length, access: accessOf(s), analyses: list};
    });
  }

  /** Three parallel reads: the study, its analyses, and every result row of the study by the
   *  denormalised `study_id` (decision 9) without its curve. The curves come later
   *  ({@link loadCurves}, listed in `pendingCurves`). Nothing is tracked: no write starts here. */
  async loadStudy(studyId: string): Promise<LoadedStudy> {
    const [study, rows, results] = await Promise.all([
      ncaDb.studies.get(studyId, {withAccess: true}) as Promise<Bag | null>,
      allPages((offset) => ncaDb.analyses.query({filter: eq('study_id', studyId),
        columns: STUDY_ANALYSIS_COLUMNS, expand: ['pipeline_run'], sort: 'label', limit: PAGE,
        offset}) as Promise<Bag[]>),
      allPages((offset) => ncaDb.analysisResults.query({filter: eq('study_id', studyId),
        columns: RESULT_COLUMNS_WITHOUT_CURVE as AnalysisResultColumn[], sort: 'analysis_id,ord', limit: PAGE,
        offset}) as Promise<Bag[]>),
    ]);
    if (study === null || study === undefined) throw new Error(`Study '${studyId}' is not available.`);
    const byAnalysis = groupBy(results, 'analysis_id');
    const pendingCurves = new Map<string, ReadonlyArray<string | null>>();
    const analyses = rows.map((r) => {
      const own = [...(byAnalysis.get(r.id) ?? [])].sort((x, y) => (x.ord ?? 0) - (y.ord ?? 0));
      if (r.has_cache === true && own.some((x) => x.has_curve === true))
        pendingCurves.set(r.id, own.map((x) => (x.has_curve === true ? x.id as string : null)));
      return entryOf(r, own);
    });
    return {...studyOf(study, analyses), ...(pendingCurves.size > 0 ? {pendingCurves} : {})};
  }

  async loadCurves(ids: readonly string[]): Promise<Map<string, CompoundCurve | null>> {
    const out = new Map<string, CompoundCurve | null>();
    for (let i = 0; i < ids.length; i += PAGE) {
      const rows = await ncaDb.analysisResults.query({filter: eq('id', ids.slice(i, i + PAGE)), columns: ['curve'],
        limit: PAGE}) as Bag[];
      for (const r of rows) out.set(r.id, (r.curve ?? null) as CompoundCurve | null);
    }
    return out;
  }

  /** The five parallel reads of one aggregate: the analysis row (with its study and run), its λz
   *  overrides, review decisions, pipeline flags and result rows (full — a later write diffs them). */
  private async readAggregate(analysisId: string): Promise<{a: Bag; results: Bag[]; lz: Bag[]; decisions: Bag[];
    flags: Bag[]}> {
    const byAnalysis = (q: (offset: number) => Promise<Bag[]>): Promise<Bag[]> => allPages(q);
    const [rows, lz, decisions, flags, results] = await Promise.all([
      ncaDb.analyses.query({filter: eq('id', analysisId), expand: ['study_id', 'pipeline_run'], withAccess: true,
        limit: 1}) as Promise<Bag[]>,
      byAnalysis((offset) => ncaDb.lzOverrides.query({filter: eq('analysis_id', analysisId), limit: PAGE,
        offset}) as Promise<Bag[]>),
      byAnalysis((offset) => ncaDb.reviewDecisions.query({filter: eq('analysis_id', analysisId), limit: PAGE,
        offset}) as Promise<Bag[]>),
      byAnalysis((offset) => ncaDb.pipelineFlags.query({filter: eq('analysis_id', analysisId), limit: PAGE,
        offset}) as Promise<Bag[]>),
      byAnalysis((offset) => ncaDb.analysisResults.query({filter: eq('analysis_id', analysisId), sort: 'ord',
        limit: PAGE, offset}) as Promise<Bag[]>),
    ]);
    if (rows.length === 0) throw new Error(`Analysis '${analysisId}' is not available.`);
    return {a: rows[0], results, lz, decisions, flags};
  }

  async loadAnalysis(analysisId: string): Promise<LoadedAnalysis> {
    let g = await this.readAggregate(analysisId);
    if ((g.a.recipe_version ?? RECIPE_VERSION) > RECIPE_VERSION)
      throw new Error('Saved by a newer NCA Studio — update the app.');
    let df: DG.DataFrame;
    try {
      df = await this.readSnapshotFile(fromFileUrl(g.a.source_file));
    } catch {
      // PR-8: the snapshot was replaced between the row read and the file read — re-read once.
      g = await this.readAggregate(analysisId);
      try {
        df = await this.readSnapshotFile(fromFileUrl(g.a.source_file));
      } catch {
        throw new Error(`The data of '${g.a.label}' could not be read.`);
      }
    }
    const {a} = g;
    this.tracked.set(a.id, trackRead(a, g.results, g.lz, g.decisions, g.flags));
    const study: StudyHeader = {
      id: a.study_id, label: a['study_id.label'], kind: a['study_id.kind'],
      updatedAt: isoOf(a['study_id.updated_at']) ?? '', analysisCount: 0,
      access: {...accessOf(a), canShare: false}, analyses: [],
    };
    return {study, entry: entryOf(a, g.results, {lz: g.lz, decisions: g.decisions, flags: g.flags}), df};
  }

  // ── studies ────────────────────────────────────────────────────────────────

  async createStudy(label: string): Promise<StudyHeader> {
    const now = nowIso();
    const [s] = await ncaDb.studies.insert({label, kind: 'study', updated_at: now});
    this.changed();
    return {id: s.id, label, kind: 'study', updatedAt: now, analysisCount: 0,
      access: {canEdit: true, canDelete: true, canShare: true}, analyses: []};
  }

  async renameStudy(studyId: string, label: string): Promise<void> {
    await ncaDb.studies.update(studyId, {label, updated_at: nowIso()} as never);
    this.changed();
  }

  async deleteStudy(studyId: string): Promise<void> {
    await this.deleteContainer(studyId);
    this.changed();
  }

  /** Delete a study (the cascade takes its analyses, their child rows and its pipeline runs) under
   *  its version — an append in between moves the version, so the analyses are re-read and the new
   *  one's snapshot is not left behind — then delete every snapshot the deleted analyses pointed at
   *  (delete is final, BD-16). */
  private async deleteContainer(studyId: string): Promise<void> {
    for (let attempt = 0; ; attempt++) {
      const study = await ncaDb.studies.get(studyId) as Bag | null;
      if (study === null || study === undefined) return; // already gone
      const rows = await allPages((offset) => ncaDb.analyses.query({filter: eq('study_id', studyId),
        columns: ['source_file'], limit: PAGE, offset}) as Promise<Bag[]>);
      try {
        await grok.dapi.domains.transaction('nca', [
          {op: 'update', table: 'study', id: studyId, expectedVersion: study.version, values: {updated_at: nowIso()}},
          {op: 'delete', table: 'study', id: studyId},
        ] as never);
      } catch (e) {
        if (!isConflict(e) || attempt + 1 >= MAX_ATTEMPTS) throw e;
        continue;
      }
      for (const r of rows) {
        this.tracked.delete(r.id);
        await this.deleteSnapshot({kind: 'stored', ref: fromFileUrl(r.source_file)});
      }
      return;
    }
  }

  async share(studyId: string, groupId: string): Promise<void> {
    const s = await ncaDb.studies.get(studyId) as Bag | null;
    if (s === null || s === undefined) throw new Error('The study is not available.');
    if (s.kind !== 'study') throw new Error('Only a study can be shared; move the analysis into a study first.');
    const [entity] = await grok.dapi.getEntities([studyId]);
    const group = await grok.dapi.groups.find(groupId);
    // View only until the platform honours row Edit grants (decision 10).
    await grok.dapi.permissions.grant(entity, group, false);
  }

  // ── analyses ───────────────────────────────────────────────────────────────

  /**
   * The ops that insert a whole aggregate (`entry` under `studyRef`), appended to `ops`; returns how
   * to track it once committed (`res` = the results of the whole transaction).
   */
  private insertOps(ops: Bag[], entry: AnalysisEntry, studyRef: string, now: string, extra: Bag = {}):
    (res: readonly Bag[], studyId: string) => Tracked {
    const ref = `a${ops.length}`;
    const cols = entryCols(entry, entry.sourceSignature.hash, RECIPE_VERSION);
    const values = {id: entry.id, ...cols, study_id: studyRef, updated_at: now,
      ...(entry.summaryCache !== null ? {results_computed_at: now} : {}), ...extra};
    const at = ops.length;
    ops.push({op: 'insert', table: 'analysis', ref, values});
    const plan = planChildren(`$${ref}`, studyRef, entry, EMPTY_CHILDREN,
      {results: true, lz: true, decisions: true, flags: true});
    const from = ops.length;
    ops.push(...plan.ops);
    return (res, studyId) => trackWrite(null, {
      version: res[at].version, studyId, label: entry.label, sourceFile: entry.source.ref,
      hash: entry.sourceSignature.hash, recipeFingerprint: computeRecipeFingerprint(recipeOf(entry)),
      base: baseOfEntry(entry), written: cols, children: plan.finish(res.slice(from, from + plan.ops.length), studyId),
    });
  }

  async createStandalone(input: NewAnalysisInput): Promise<{studyId: string; analysisId: string; label: string;
    source: SourcePointer}> {
    const taken = new Set((await allPages((offset) => ncaDb.studies.query({filter: eq('kind', 'standalone'),
      columns: ['label'], sort: 'label', limit: PAGE, offset}) as Promise<Bag[]>)).map((r) => r.label as string));
    const label = uniqueLabel(input.label.trim() || 'Analysis', taken);
    const snap = await this.writeSnapshot(input.df, nowIso());
    const id = input.id ?? newId();
    const entry = entryFromInput(id, {...input, label}, snap);
    const now = nowIso();
    try {
      const ops: Bag[] = [{op: 'insert', table: 'study', ref: 's',
        values: {label, kind: 'standalone', updated_at: now}}];
      const track = this.insertOps(ops, entry, '$s', now);
      const res = await grok.dapi.domains.transaction('nca', ops as never) as Bag[];
      const studyId = res[0].id as string;
      this.tracked.set(id, track(res, studyId));
      this.changed();
      return {studyId, analysisId: id, label, source: snap.source};
    } catch (e) {
      if (isDefinitiveRefusal(e)) await this.deleteSnapshot(snap.source);
      throw e;
    }
  }

  async appendAnalysis(studyId: string, input: NewAnalysisInput): Promise<{analysisId: string; label: string;
    source: SourcePointer}> {
    const taken = await this.analysisLabels(studyId);
    const label = uniqueLabel(input.label.trim() || 'Analysis', taken);
    const snap = await this.writeSnapshot(input.df, nowIso());
    const id = input.id ?? newId();
    const entry = entryFromInput(id, {...input, label}, snap);
    const now = nowIso();
    try {
      const ops: Bag[] = [{op: 'update', table: 'study', id: studyId, values: {updated_at: now}}];
      const track = this.insertOps(ops, entry, studyId, now);
      const res = await grok.dapi.domains.transaction('nca', ops as never) as Bag[];
      this.tracked.set(id, track(res, studyId));
      this.changed();
      return {analysisId: id, label, source: snap.source};
    } catch (e) {
      if (isDefinitiveRefusal(e)) await this.deleteSnapshot(snap.source);
      throw e;
    }
  }

  /** Every analysis label of a study, all pages (SD-6: a label check must not stop at one page). */
  private async analysisLabels(studyId: string): Promise<Set<string>> {
    return new Set((await allPages((offset) => ncaDb.analyses.query({filter: eq('study_id', studyId),
      columns: ['label'], sort: 'label', limit: PAGE, offset}) as Promise<Bag[]>)).map((r) => r.label as string));
  }

  private async trackedOf(analysisId: string): Promise<Tracked> {
    const t = this.tracked.get(analysisId);
    if (t !== undefined) return t;
    return this.refresh(analysisId);
  }

  /** Re-read one aggregate into the tracker. */
  private async refresh(analysisId: string): Promise<Tracked> {
    let g;
    try {
      g = await this.readAggregate(analysisId);
    } catch {
      throw new Error('The analysis no longer exists.');
    }
    const t = trackRead(g.a, g.results, g.lz, g.decisions, g.flags);
    this.tracked.set(analysisId, t);
    return t;
  }

  async saveAnalysis(analysisId: string, expectedRef: string, input: SaveAnalysisInput): Promise<SourcePointer> {
    let t = await this.trackedOf(analysisId);
    if (t.sourceFile !== expectedRef) throw new SourceReplacedError(t.sourceFile, expectedRef);
    // SD-19 / PR-13: a new snapshot only when the data changed (full content hash).
    const hash = computeSnapshotHash(input.df);
    const snap: WrittenSnapshot | null = hash === t.hash ? null : await this.writeSnapshot(input.df, nowIso());
    const source: SourcePointer = snap?.source ?? {kind: 'stored', ref: t.sourceFile};
    const signature: SourceSignature = snap?.signature ??
      {kind: 'content-hash', rowCount: input.df.rowCount, mtime: null, hash: t.hash, computedAt: nowIso()};
    try {
      for (let attempt = 0; ; attempt++) {
        const now = nowIso();
        const entry: AnalysisEntry = {
          ...input.recipe, id: analysisId, label: t.label, source, sourceSignature: signature,
          summaryCache: input.summaryCache,
          engineVersion: input.summaryCache === null ? null : INSTALLED_ENGINE_VERSION,
          resultsBasis: input.summaryCache === null ? null : RESULTS_BASIS,
          derivedBasis: input.summaryCache === null ? null : DERIVED_BASIS,
          confirmed: input.workflow.confirmed, lzOverrides: input.workflow.lzOverrides,
          reviewDecisions: input.workflow.reviewDecisions, stageCache: input.workflow.stageCache,
        };
        const all = entryCols(entry, snap?.signature.hash ?? t.hash, RECIPE_VERSION);
        const owned = Object.fromEntries(Object.entries(all).filter(([k]) =>
          !NOT_SAVED.has(k) && (snap !== null || !SNAPSHOT_COLUMNS.has(k))));
        const values = changedCols(owned, t);
        const children = planChildren(analysisId, t.studyId, entry, t, {results: true, lz: true, decisions: true});
        const resultsChanged = children.resultsChanged || Object.keys(values).some((k) => RESULT_ENTRY_COLUMNS.has(k));
        const ops: Bag[] = [
          {op: 'update', table: 'analysis', id: analysisId, expectedVersion: t.version,
            values: {...values, updated_at: now, ...(resultsChanged ? {results_computed_at: now} : {})}},
          {op: 'update', table: 'study', id: t.studyId, values: {updated_at: now}},
          ...children.ops,
        ];
        try {
          const res = await grok.dapi.domains.transaction('nca', ops as never) as Bag[];
          const prev = t.sourceFile;
          this.tracked.set(analysisId, trackWrite(t, {
            version: res[0].version, studyId: t.studyId, label: t.label, sourceFile: source.ref,
            hash: snap?.signature.hash ?? t.hash, recipeFingerprint: computeRecipeFingerprint(input.recipe),
            base: baseOfEntry(entry), written: owned, children: children.finish(res.slice(2), t.studyId),
          }));
          if (snap !== null && prev !== '') await this.deleteSnapshot({kind: 'stored', ref: prev});
          this.changed();
          return source;
        } catch (e) {
          if (!isConflict(e)) throw e;
          const base = t.base;
          t = await this.refresh(analysisId);
          if (t.sourceFile !== expectedRef) throw new SourceReplacedError(t.sourceFile, expectedRef);
          if (t.base !== base) throw new ContentChangedError(t.label, null);
          if (attempt + 1 >= MAX_ATTEMPTS) throw new KeepsChangingError(t.label);
        }
      }
    } catch (e) {
      if (snap !== null && (e instanceof SourceReplacedError || e instanceof ContentChangedError ||
        isDefinitiveRefusal(e))) await this.deleteSnapshot(snap.source);
      throw e;
    }
  }

  /**
   * Background results refresh of a clean view (SD-18, PR-2), diff-gated against the tracker: only
   * the columns and child rows that differ are written, and nothing at all when nothing differs (a
   * plain reopen — decision 2). Otherwise ONE transaction guarded by the `analysis` version
   * (decision 14); `results_computed_at` moves only when the results changed. A version race drops
   * the write (the next clean compute writes again); any other failure is the caller's to show.
   */
  async writeSummary(analysisId: string, expectedRef: string, write: SummaryWrite): Promise<void> {
    const t = this.tracked.get(analysisId);
    if (t === undefined || t.sourceFile !== expectedRef) return; // replaced / never read here
    if (t.hash !== null && write.frameHash !== t.hash) return; // PR-2: computed from unsaved data
    const w = write.workflow;
    const now = nowIso();
    // A freshly computed summary carries the installed engine and the current basis (AD-8).
    const x: ChildInput = {summaryCache: write.summaryCache, lzOverrides: w.lzOverrides,
      reviewDecisions: w.reviewDecisions, engineVersion: INSTALLED_ENGINE_VERSION, resultsBasis: RESULTS_BASIS};
    const cols: Bag = {engine_version: INSTALLED_ENGINE_VERSION, results_basis: RESULTS_BASIS,
      derived_basis: DERIVED_BASIS,
      ...tierCols(write.summaryCache), results_snapshot_hash: t.hash,
      ...confirmedCols(w.confirmed, t.recipeFingerprint), ...stageCols(w.stageCache)};
    const values = changedCols(cols, t);
    const children = planChildren(analysisId, t.studyId, x, t, {results: true, lz: true, decisions: true});
    const resultsChanged = children.resultsChanged || Object.keys(values).some((k) => RESULT_ENTRY_COLUMNS.has(k));
    if (Object.keys(values).length === 0 && children.ops.length === 0) return;
    const guard: Bag = {...values, ...(resultsChanged ? {results_computed_at: now} : {})};
    // A write of child rows alone still moves the analysis version (decision 14).
    if (Object.keys(guard).length === 0) guard.updated_at = now;
    const ops: Bag[] = [{op: 'update', table: 'analysis', id: analysisId, expectedVersion: t.version, values: guard},
      ...children.ops];
    let res: Bag[];
    try {
      res = await grok.dapi.domains.transaction('nca', ops as never) as Bag[];
    } catch (e) {
      if (isConflict(e)) return; // SD-18: dropped; the next clean compute writes again
      throw e;
    }
    // PR-12: this session's own write moves the base, so its next Save is not mistaken for a foreign one.
    const base = JSON.parse(t.base) as Bag;
    this.tracked.set(analysisId, trackWrite(t, {
      version: res[0].version, studyId: t.studyId, label: t.label, sourceFile: t.sourceFile, hash: t.hash,
      recipeFingerprint: t.recipeFingerprint,
      base: baseOf(base.recipe, w.confirmed?.fingerprint ?? null, w.confirmed?.at ?? null, w.lzOverrides,
        w.reviewDecisions),
      written: cols, children: children.finish(res.slice(1), t.studyId),
    }));
  }

  async renameAnalysis(analysisId: string, label: string): Promise<void> {
    const rows = await ncaDb.analyses.query({filter: eq('id', analysisId), expand: ['study_id'],
      columns: ['study_id', 'label'], limit: 1}) as Bag[];
    if (rows.length === 0) throw new Error('The analysis no longer exists.');
    const r = rows[0];
    const now = nowIso();
    const standalone = r['study_id.kind'] === 'standalone';
    await grok.dapi.domains.transaction('nca', [
      {op: 'update', table: 'analysis', id: analysisId, values: {label, updated_at: now}},
      {op: 'update', table: 'study', id: r.study_id, values: standalone ? {label, updated_at: now} : {updated_at: now}},
    ] as never);
    // Keep the open's version and content base (PR-1): the next Save loses one version race on this
    // rename, re-reads, and still refuses if another session changed the content in between. Dropping
    // the tracker would take the re-read row as the base and overwrite that change.
    const t = this.tracked.get(analysisId);
    if (t !== undefined) this.tracked.set(analysisId, {...t, label, cols: {...t.cols, label: canonical(label)}});
    this.changed();
  }

  async deleteAnalysis(analysisId: string): Promise<void> {
    const rows = await ncaDb.analyses.query({filter: eq('id', analysisId), expand: ['study_id'],
      columns: ['study_id', 'source_file'], limit: 1}) as Bag[];
    if (rows.length === 0) return;
    const r = rows[0];
    if (r['study_id.kind'] === 'standalone') await this.deleteContainer(r.study_id);
    else {
      // The cascade takes the analysis's result, override, decision and flag rows.
      await grok.dapi.domains.transaction('nca', [
        {op: 'delete', table: 'analysis', id: analysisId},
        {op: 'update', table: 'study', id: r.study_id, values: {updated_at: nowIso()}},
      ] as never);
      await this.deleteSnapshot({kind: 'stored', ref: fromFileUrl(r.source_file)});
    }
    this.tracked.delete(analysisId);
    this.changed();
  }

  /**
   * A re-parent, not a copy (architecture § "Storage on EMS"): the same analysis id moves into the
   * study, with its child rows — the denormalised `study_id` of its result and decision rows is
   * updated in the same transaction (`lz_override` has none), and every `set_*` / `decided_*` stays.
   */
  async moveToStudy(standaloneId: string, target: MoveTarget, opts: {
    label?: string; passthroughCols?: MoveResult['passthroughCols']; summaryCache?: SummaryCache;
  } = {}): Promise<MoveResult> {
    const [rows, resultRows, decisionRows] = await Promise.all([
      ncaDb.analyses.query({filter: eq('study_id', standaloneId), expand: ['study_id'], limit: 2}) as Promise<Bag[]>,
      // Whole rows when a fresh summary comes along: it is compared with them (decision 2).
      allPages((offset) => ncaDb.analysisResults.query({filter: eq('study_id', standaloneId),
        columns: opts.summaryCache !== undefined ? ['analysis_id', ...RESULT_COLUMNS] as AnalysisResultColumn[] :
          ['ord', 'analysis_id'],
        sort: 'ord', limit: PAGE, offset}) as Promise<Bag[]>),
      allPages((offset) => ncaDb.reviewDecisions.query({filter: eq('study_id', standaloneId), columns: ['unit_key'],
        limit: PAGE, offset}) as Promise<Bag[]>),
    ]);
    const r = rows[0];
    if (r === undefined || r['study_id.kind'] !== 'standalone')
      throw new Error('The analysis is no longer in Ungrouped.');
    const containerLabel = r['study_id.label'] as string;
    const stored = r.recipe as AnalysisRecipe;
    const recipe: AnalysisRecipe = {...stored, passthroughCols: opts.passthroughCols ?? stored.passthroughCols};
    // AC-S1.7.2: a study entry always carries results.
    if (opts.summaryCache === undefined && r.has_cache !== true)
      throw new Error('Confirm the setup first: an analysis without results cannot join a study.');
    const now = nowIso();
    let targetId: string;
    let targetLabel: string;
    const ops: Bag[] = [];
    let taken = new Set<string>();
    if (target.kind === 'existing') {
      // The target's labels + version (the guard below).
      const [t, labels] = await Promise.all([ncaDb.studies.get(target.studyId) as Promise<Bag | null>,
        this.analysisLabels(target.studyId)]);
      if (t === null || t === undefined) throw new Error('The target study is not available.');
      targetId = t.id;
      targetLabel = t.label;
      taken = labels;
      // The guard (AC-3.5): this touch carries the target's version — a View-only caller gets 403.
      ops.push({op: 'update', table: 'study', id: targetId, expectedVersion: t.version, values: {updated_at: now}});
    } else {
      const studies = new Set((await allPages((offset) => ncaDb.studies.query({filter: eq('kind', 'study'),
        columns: ['label'], sort: 'label', limit: PAGE, offset}) as Promise<Bag[]>)).map((x) => x.label as string));
      targetLabel = uniqueLabel(target.label.trim() || 'New study', new Set([...studies, STANDALONE_ROUTE_SEGMENT]));
      targetId = '$t';
      ops.push({op: 'insert', table: 'study', ref: 't', values: {label: targetLabel, kind: 'study', updated_at: now}});
    }
    const label = uniqueLabel((opts.label ?? r.label).trim() || r.label, taken);
    // A fresh summary (a move straight after a compute) rewrites the result rows in place; the
    // compute time moves only if the results differ from the stored ones (decision 2).
    const own = resultRows.filter((x) => x.analysis_id === r.id);
    let fresh: Bag = {};
    if (opts.summaryCache !== undefined) {
      const stamps = {engineVersion: INSTALLED_ENGINE_VERSION, resultsBasis: RESULTS_BASIS};
      const cols: Bag = {...tierCols(opts.summaryCache), engine_version: stamps.engineVersion,
        results_basis: stamps.resultsBasis, derived_basis: DERIVED_BASIS, results_snapshot_hash: r.source_hash ?? null};
      const before = canonicalColumns(r, Object.keys(cols));
      const after = canonicalColumns(cols, Object.keys(cols));
      const same = (x: Bag): string => canonicalRow({...x, study_id: ''}, RESULT_COLUMNS);
      const rowsNew = resultRowsOf(r.id, '', opts.summaryCache, stamps);
      const changed = Object.keys(cols).some((k) => before[k] !== after[k]) || rowsNew.length !== own.length ||
        rowsNew.some((x, i) => same(x) !== same(own[i]));
      fresh = {...cols, ...(changed ? {results_computed_at: now} : {})};
    }
    // Re-parent first, then delete the emptied container (S0 AC-0.4: no cascade onto the analysis).
    ops.push({op: 'update', table: 'analysis', id: r.id, values: {
      study_id: targetId, label, recipe, ...recipeSearchCols(recipe), moved_from_label: containerLabel,
      updated_at: now, ...originCols({kind: 'standalone', label: containerLabel, containerId: standaloneId}), ...fresh,
    }});
    if (opts.summaryCache !== undefined) {
      const plan = planChildren(r.id, targetId, {summaryCache: opts.summaryCache,
        engineVersion: INSTALLED_ENGINE_VERSION, resultsBasis: RESULTS_BASIS},
      {...EMPTY_CHILDREN, results: own.map((x) => ({id: x.id, canon: null}))}, {results: true});
      ops.push(...plan.ops);
    } else
      for (const x of own) ops.push({op: 'update', table: 'analysis_result', id: x.id, values: {study_id: targetId}});
    for (const d of decisionRows)
      ops.push({op: 'update', table: 'review_decision', id: d.id, values: {study_id: targetId}});
    ops.push(
      {op: 'update', table: 'study', id: standaloneId, values: {moved_to: targetId}},
      {op: 'delete', table: 'study', id: standaloneId},
    );
    const res = await grok.dapi.domains.transaction('nca', ops as never) as Bag[];
    if (target.kind === 'new') targetId = res[0].id as string;
    this.tracked.delete(r.id);
    this.changed();
    return {studyId: targetId, studyLabel: targetLabel, analysisId: r.id, label,
      source: {kind: 'stored', ref: fromFileUrl(r.source_file)}, passthroughCols: recipe.passthroughCols};
  }

  // ── pipeline ───────────────────────────────────────────────────────────────

  /** The pipeline's own study by business key, then its whole aggregates in five parallel reads:
   *  the child tables filter through `analysis_id.study_id` (or their own `study_id`), so none waits
   *  for the analyses (FIELD-10's "untouched" test needs the overrides and decisions). */
  async pipelineBase(owner: string, studyKey: string): Promise<PipelineBase> {
    const s = await ncaDb.studies.getByKey({pipeline_owner: owner, pipeline_study_key: studyKey}) as Bag | null;
    if (s === null || s === undefined) return {study: undefined};
    const pages = (q: (offset: number) => Promise<Bag[]>): Promise<Bag[]> => allPages(q);
    const [rows, results, lz, decisions, flags] = await Promise.all([
      pages((offset) => ncaDb.analyses.query({filter: eq('study_id', s.id), expand: ['pipeline_run'], limit: PAGE,
        offset}) as Promise<Bag[]>),
      pages((offset) => ncaDb.analysisResults.query({filter: eq('study_id', s.id), sort: 'analysis_id,ord',
        limit: PAGE, offset}) as Promise<Bag[]>),
      pages((offset) => ncaDb.lzOverrides.query({filter: eq('analysis_id.study_id', s.id), limit: PAGE,
        offset}) as Promise<Bag[]>),
      pages((offset) => ncaDb.reviewDecisions.query({filter: eq('study_id', s.id), limit: PAGE,
        offset}) as Promise<Bag[]>),
      pages((offset) => ncaDb.pipelineFlags.query({filter: eq('analysis_id.study_id', s.id), limit: PAGE,
        offset}) as Promise<Bag[]>),
    ]);
    const [byResults, byLz, byDecisions, byFlags] =
      [results, lz, decisions, flags].map((x) => groupBy(x, 'analysis_id'));
    const tracked = new Map<string, Tracked>();
    const analyses = rows.map((r) => {
      const c = {results: byResults.get(r.id) ?? [], lz: byLz.get(r.id) ?? [], decisions: byDecisions.get(r.id) ?? [],
        flags: byFlags.get(r.id) ?? []};
      tracked.set(r.id, trackRead(r, c.results, c.lz, c.decisions, c.flags));
      return entryOf(r, c.results, c);
    });
    const study: NcaStudy = {...studyOf({...s, '~can_edit': true, '~can_delete': true}, analyses),
      pipelineKey: {owner, studyKey}};
    this.pipelineReads.set(study, {study: s.version, analyses: tracked});
    return {study};
  }

  async savePipelinePlan(write: PipelineWrite): Promise<string> {
    const now = nowIso();
    const baseIds = new Set((write.base?.analyses ?? []).map((a) => a.id));
    const seen = write.base === undefined ? undefined : this.pipelineReads.get(write.base);
    if (write.base !== undefined && seen === undefined)
      throw new Error('savePipelinePlan: the base was not read with pipelineBase');
    const ops: Bag[] = [];
    let studyRef: string;
    if (write.base === undefined) {
      ops.push({op: 'insert', table: 'study', ref: 's', values: {label: write.planned.label, kind: 'study',
        pipeline_owner: write.owner, pipeline_study_key: write.studyKey, updated_at: now}, onDuplicate: 'error'});
      studyRef = '$s';
    } else {
      studyRef = write.base.id;
      // The study version of the base read: any other write to the study since → re-plan.
      ops.push({op: 'update', table: 'study', id: studyRef, expectedVersion: seen!.study, values: {updated_at: now}});
    }
    // One pipeline_run row per (run, study), stated once for every analysis of that run written here
    // (decision 12). An insert merges into an existing row by its business key, and its `$ref` then
    // names that row (probe "a merged insert …").
    const runRefs = new Map<string, string>();
    const runRef = (o: PipelineOrigin): string => {
      let ref = runRefs.get(o.runId);
      if (ref === undefined) {
        ref = `r${runRefs.size + 1}`;
        runRefs.set(o.runId, ref);
        ops.push({op: 'insert', table: 'pipeline_run', ref, values: pipelineRunCols(o, studyRef)});
      }
      return `$${ref}`;
    };
    for (const e of write.planned.analyses) {
      const run = e.origin?.kind === 'pipeline' ? {pipeline_run: runRef(e.origin)} : {pipeline_run: null};
      if (baseIds.has(e.id)) {
        const before = write.base!.analyses.find((a) => a.id === e.id)!;
        if (canonical(before) === canonical(e)) continue; // untouched by this run
        const t = seen!.analyses.get(e.id)!;
        const cols = entryCols(e, e.sourceSignature.hash, RECIPE_VERSION);
        const values = changedCols(cols, t);
        const children = planChildren(e.id, studyRef, e, t, {results: true, lz: true, decisions: true, flags: true});
        const resultsChanged = children.resultsChanged ||
          Object.keys(values).some((k) => RESULT_ENTRY_COLUMNS.has(k));
        ops.push({op: 'update', table: 'analysis', id: e.id, expectedVersion: t.version, values: {...values,
          ...run, study_id: studyRef, updated_at: now, ...(resultsChanged ? {results_computed_at: now} : {})}});
        ops.push(...children.ops);
      } else
        this.insertOps(ops, e, studyRef, now, run);
    }
    try {
      const res = await grok.dapi.domains.transaction('nca', ops as never) as Bag[];
      for (const e of write.planned.analyses) this.tracked.delete(e.id); // versions moved
      this.changed();
      return write.base === undefined ? res[0].id as string : write.base.id;
    } catch (e) {
      // A race on the base, or another run inserting the same pipeline key first (business-key duplicate).
      // The raw unique violation (23505) is how S0 saw a key race surface; it only routes to a re-plan.
      if (isConflict(e) || (e instanceof DG.DomainValidationError && e.isDuplicate) ||
        /23505/.test(e instanceof Error ? e.message : String(e)))
        throw new KeepsChangingError(write.planned.label);
      throw e;
    }
  }

  // ── snapshots ──────────────────────────────────────────────────────────────

  async writeSnapshot(frame: DG.DataFrame, computedAt: string): Promise<WrittenSnapshot> {
    const df = snapshotFrame(frame); // PR-13: the source columns only
    const bytes = df.toByteArray();
    // BD-3: decode our own bytes before uploading; refuse rather than store a different table.
    const diff = snapshotDiff(df, DG.DataFrame.fromByteArray(bytes));
    if (diff !== null) throw new SnapshotFidelityError(diff);
    const path = `${SNAPSHOT_ROOT}/${newId()}/source.d42`;
    await grok.dapi.files.write(path, bytes);
    return {source: {kind: 'stored', ref: path}, signature: snapshotSignature(df, computedAt)};
  }

  async deleteSnapshot(source: SourcePointer): Promise<void> {
    if (source.kind !== 'stored' || !source.ref.startsWith(SNAPSHOT_ROOT)) return;
    try {
      await grok.dapi.files.delete(source.ref.slice(0, source.ref.lastIndexOf('/')));
    } catch {/* an orphan is harmless */}
  }

  private async readSnapshotFile(path: string): Promise<DG.DataFrame> {
    return DG.DataFrame.fromByteArray(await grok.dapi.files.readAsBytes(path));
  }

  async readSource(source: SourcePointer): Promise<DG.DataFrame> {
    switch (source.kind) {
    case 'stored':
      return this.readSnapshotFile(source.ref);
    case 'tableId': {
      const df = grok.shell.tables.find((t) => t.id === source.ref || t.name === source.ref);
      if (df === undefined) throw new Error(`table '${source.ref}' is not open`);
      return df;
    }
    default:
      throw new Error(`${source.kind} sources are not supported in this version`);
    }
  }
}

/** The domain study of a study row and its entries. */
function studyOf(s: Bag, analyses: AnalysisEntry[]): NcaStudy {
  return {
    id: s.id, label: s.label, kind: s.kind, createdAt: isoOf(s.created_on) ?? '',
    updatedAt: isoOf(s.updated_at ?? s.updated_on) ?? '', access: accessOf(s), analyses,
    ...(s.pipeline_owner ?
      {pipelineKey: {owner: s.pipeline_owner as string, studyKey: s.pipeline_study_key as string}} : {}),
  };
}
