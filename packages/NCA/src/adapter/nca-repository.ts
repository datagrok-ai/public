/**
 * The NCA storage contract. Every read and write of studies,
 * analyses and their source snapshots goes through ONE {@link NcaRepository}; the views and the
 * headless pipeline never touch the platform storage APIs themselves.
 *
 * Two implementations: `EmsNcaRepository` (the platform EMS schema `nca` + d42 snapshots in a
 * `file` column, the production store) and an in-memory one for tests (decision 7). One contract
 * suite runs against both (AC-3.8), so the fake cannot drift from the EMS behaviour.
 *
 * **The caller names the container** (guardrails, BUG-16): every method takes the study /
 * analysis id it acts on. The repository holds no "current study". What it DOES hold, per
 * analysis id, is the last `version`, `source_file` and content base it read or wrote — so every
 * write response rebinds the next write (PR-12) without threading a mutable handle through
 * `view.temp`.
 */

import * as DG from 'datagrok-api/dg';

import type {
  NcaStudy, AnalysisEntry, AnalysisRecipe, SummaryCache, StageCache, ConfirmedRecipe, LzOverride,
  ReviewLedger, AnalysisOrigin, PipelineFlag, SourcePointer, SourceSignature, PassthroughColumnSpec,
  CompoundCurve,
} from './study-schema';
import {INSTALLED_ENGINE_VERSION} from './engine-version';
import {RESULTS_BASIS} from './unit-conversion';
import {DERIVED_BASIS} from './study-schema';
import {isOccasionName} from './occasion';

// ─── Shared constants (moved from the deleted StudyCoordinator) ─────────────────

/** The URL segment of the standalone-analyses route (`/apps/Nca/analyses/<label>`). A RESERVED
 *  study label: no study may be named so, or its overview link would be shadowed. */
export const STANDALONE_ROUTE_SEGMENT = 'analyses';

/** The app's URL root (`/apps/<package name>`). Here, not in the router, so the headless
 *  pipeline can build its output links without importing the view layer. */
export const NCA_APP_ROOT = '/apps/Nca';

/** Fired after every structural change (create / rename / delete / add / move) so the
 *  Hub lists and the Browse tree rebuild. NOT fired by a background summary refresh. */
export const NCA_STUDIES_CHANGED_EVENT = 'nca-studies-changed';

/** A new id (`crypto.randomUUID` when available). */
export function newId(): string {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const c = (globalThis as any).crypto;
  if (c?.randomUUID) return c.randomUUID() as string;
  let s = '';
  for (let i = 0; i < 32; i++) s += Math.floor(Math.random() * 16).toString(16);
  return s;
}

/** `desired`, or `desired-2`, `desired-3`, … — the first one not in `taken`. Label uniqueness is
 *  an APP rule over the caller's visible labels in one namespace (SD-11): studies among studies,
 *  standalone analyses among standalone analyses, analyses within their study. */
export function uniqueLabel(desired: string, taken: ReadonlySet<string>): string {
  if (!taken.has(desired)) return desired;
  for (let n = 2; ; n++) {
    const candidate = `${desired}-${n}`;
    if (!taken.has(candidate)) return candidate;
  }
}

/** What the staleness predicate reads of an analysis (a header or an entry). */
export interface StalenessStamp {
  readonly engineVersion?: string | null;
  readonly resultsBasis?: string | null;
  readonly hasCache: boolean;
  /** tk-study-day-stratification F6 — read by {@link isDerivedBasisStale} with the recipe parts. */
  readonly derivedBasis?: string | null;
}

/**
 * Was the cached summary computed under another reporting basis of CL / Vz / Vss
 * (nca-cl-bw-normalised-units AD-8)? A missing stamp predates the per-kg basis. Half of
 * {@link isStale}; read on its own only by the study grid, which blanks the basis-dependent cells
 * of such a row (an old absolute mean must not sit under a mL/min/kg header — architect B3).
 */
export function isBasisStale(a: StalenessStamp): boolean {
  if (!a.hasCache) return false;
  return a.resultsBasis == null || a.resultsBasis !== RESULTS_BASIS;
}

/**
 * Is an analysis's cached summary STALE? The ONE staleness predicate every reader of the cache
 * uses, so the cue cannot drift between surfaces: computed by a different sci-comp version than the installed
 * one, or by an unknown one — OR under another reporting basis ({@link isBasisStale}). An
 * analysis with no cached numbers is never stale — there is nothing to predate a release.
 */
export function isStale(a: StalenessStamp, installed: string = INSTALLED_ENGINE_VERSION): boolean {
  if (!a.hasCache) return false;
  return a.engineVersion == null || a.engineVersion !== installed || isBasisStale(a) || isDerivedBasisStale(a);
}

/** The parts of an entry {@link isDerivedBasisStale} reads (all optional on a bare header). */
type DerivedStampSource = Partial<Pick<AnalysisEntry, 'derivedBasis' | 'summaryCache' | 'groupByCols' |
  'passthroughCols' | 'derived' | 'columns'>>;

/**
 * tk-study-day-stratification F6 (SF-STALE, narrowed by the user 2026-11-10): was a cached %F
 * computed before %F was stratified per occasion, on an analysis where that moves it? True iff
 * the cache holds a %F value, the analysis is multi-occasion and `derivedBasis` is not
 * `DERIVED_BASIS`. A fresh entry always carries the current stamp, so only legacy entries are
 * tested, and "multi-occasion" for them = the ONE occasion-name matcher (`isOccasionName`)
 * matches a column the recipe names. Accepted residuals (cue-only): a legacy repeat-dose
 * analysis whose occasion column the recipe never named is not marked; a legacy analysis grouped
 * by `Period` or by a single-valued `Day` reads stale although its %F does not move.
 */
export function isDerivedBasisStale(a: StalenessStamp & DerivedStampSource): boolean {
  if (!a.hasCache || a.derivedBasis === DERIVED_BASIS) return false;
  const tiers = a.summaryCache?.compoundTier ?? [];
  if (!tiers.some((t) => t.fPct !== null && Number.isFinite(t.fPct.mean))) return false;
  const c = a.columns;
  const named = [
    ...(a.groupByCols ?? []), ...(a.passthroughCols ?? []).map((p) => p.name),
    ...(a.derived?.occasionCol !== undefined ? [a.derived.occasionCol] : []),
    ...(c !== undefined ? Object.values(c).filter((v): v is string => typeof v === 'string') : []),
  ];
  return named.some((n) => isOccasionName(n));
}

/** {@link isStale} for a full entry (its cache present or not). */
export function isEntryStale(e: Pick<AnalysisEntry, 'engineVersion' | 'resultsBasis' | 'summaryCache'> &
  DerivedStampSource, installed: string = INSTALLED_ENGINE_VERSION): boolean {
  return isStale({...e, hasCache: e.summaryCache !== null}, installed);
}

// ─── Snapshot hash (PR-3 / PR-13) ───────────────────────────────────────────────

/** FNV-1a 32-bit over bytes with a seed, so four passes give 128 bits. */
function fnvBytes(bytes: Uint8Array, seed: number): number {
  let h = (0x811c9dc5 ^ seed) >>> 0;
  for (let i = 0; i < bytes.length; i++) {
    h ^= bytes[i];
    h = Math.imul(h, 0x01000193);
  }
  return h >>> 0;
}

/** A column the app derives on the source frame during a compute (`~RouteCode`, `~BLQ`,
 *  `~InfusionDuration` — every derived column carries the `~` prefix). Not data: it is re-derived
 *  on every open. */
const isDerivedColumn = (c: DG.Column): boolean => c.name.startsWith('~');

/**
 * The frame a snapshot stores (PR-13): the source columns only, never the app-derived ones. The
 * frame itself when it has none; otherwise a clone without them.
 */
export function snapshotFrame(df: DG.DataFrame): DG.DataFrame {
  const cols = df.columns.toList();
  if (!cols.some(isDerivedColumn)) return df;
  const out = df.clone(null, cols.filter((c) => !isDerivedColumn(c)).map((c) => c.name), false, true);
  out.name = df.name;
  return out;
}

/**
 * The content hash of a source snapshot (PR-3 / PR-13): every row of every SOURCE column, in order
 * (exactly what {@link snapshotFrame} stores — a derived `~` column added by a compute must not make
 * the frame read as edited). The canonical bytes are the column count, then per column its UTF-8
 * name, its DG type and its raw
 * values — a float as its 8 IEEE bytes, an int as int32, a missing cell as a distinct marker byte,
 * anything else as its length-prefixed UTF-8 `getString`. 128 bits from four seeded FNV-1a passes:
 * change detection, not security, and synchronous on every install (`crypto.subtle` is absent on
 * plain-HTTP servers).
 */
export function computeSnapshotHash(df: DG.DataFrame): string {
  const chunks: Uint8Array[] = [];
  const enc = new TextEncoder();
  const u32 = (n: number): Uint8Array => {
    const b = new Uint8Array(4);
    new DataView(b.buffer).setUint32(0, n >>> 0);
    return b;
  };
  const str = (s: string): void => {
    const b = enc.encode(s);
    chunks.push(u32(b.length), b);
  };
  const source = df.columns.toList().filter((c) => !isDerivedColumn(c));
  chunks.push(u32(source.length), u32(df.rowCount));
  const NULL = Uint8Array.of(0);
  const SET = Uint8Array.of(1);
  for (const col of source) {
    str(col.name);
    str(col.type);
    const isFloat = col.type === DG.COLUMN_TYPE.FLOAT || col.type === DG.COLUMN_TYPE.QNUM;
    const isInt = col.type === DG.COLUMN_TYPE.INT;
    for (let i = 0; i < df.rowCount; i++) {
      if (col.isNone(i)) {chunks.push(NULL); continue;}
      chunks.push(SET);
      if (isFloat) {
        const b = new Uint8Array(8);
        new DataView(b.buffer).setFloat64(0, Number(col.get(i)));
        chunks.push(b);
      } else if (isInt) {
        const b = new Uint8Array(4);
        new DataView(b.buffer).setInt32(0, Number(col.get(i)));
        chunks.push(b);
      } else
        str(col.getString(i));
    }
  }
  let total = 0;
  for (const c of chunks) total += c.length;
  const all = new Uint8Array(total);
  let off = 0;
  for (const c of chunks) {all.set(c, off); off += c.length;}
  return [1, 2, 3, 4].map((s) => fnvBytes(all, s).toString(16).padStart(8, '0')).join('');
}

/** The source signature of a stored snapshot: row count + {@link computeSnapshotHash}. Written
 *  ONLY together with the snapshot it describes (SD-19) — never re-stamped by a cache refresh. */
export function snapshotSignature(df: DG.DataFrame, computedAt: string): SourceSignature {
  return {kind: 'content-hash', rowCount: df.rowCount, mtime: null, hash: computeSnapshotHash(df), computedAt};
}

/**
 * The first difference between two frames as data, or null: row count, column order, names,
 * types, semTypes, data tags and every value (floats by `Object.is`, so -0 and NaN count).
 * `.`-prefixed tags (the platform's per-session bookkeeping, e.g. `.id`) and null-valued tags are
 * not data. Used by the snapshot writer's verify-on-write (BD-3) and the S0 probes.
 */
export function snapshotDiff(a: DG.DataFrame, b: DG.DataFrame): string | null {
  if (a.rowCount !== b.rowCount) return `rows ${a.rowCount} vs ${b.rowCount}`;
  const ac = a.columns.toList();
  const bc = b.columns.toList();
  if (ac.length !== bc.length) return `columns ${ac.length} vs ${bc.length}`;
  const dataTags = (c: DG.Column): string => JSON.stringify(Object.fromEntries(
    Object.entries(c.tags ?? {}).filter(([k, v]) => !k.startsWith('.') && v != null).sort()));
  for (let j = 0; j < ac.length; j++) {
    const x = ac[j]; const y = bc[j];
    if (x.name !== y.name) return `column ${j}: '${x.name}' vs '${y.name}'`;
    if (x.type !== y.type) return `'${x.name}': type ${x.type} vs ${y.type}`;
    if ((x.semType ?? '') !== (y.semType ?? '')) return `'${x.name}': semType ${x.semType} vs ${y.semType}`;
    if (dataTags(x) !== dataTags(y)) return `'${x.name}': tags ${dataTags(x)} vs ${dataTags(y)}`;
    for (let i = 0; i < a.rowCount; i++) {
      if (x.isNone(i) !== y.isNone(i)) {
        const shown = (c: DG.Column): string => c.isNone(i) ? 'missing' : c.getString(i);
        return `'${x.name}' row ${i + 1}: ${shown(x)} vs ${shown(y)}`;
      }
      if (x.isNone(i)) continue;
      const xv = x.get(i); const yv = y.get(i);
      if (typeof xv === 'number' && typeof yv === 'number') {
        if (!Object.is(xv, yv)) return `'${x.name}' row ${i + 1}: ${xv} vs ${yv}`;
      } else if (x.getString(i) !== y.getString(i))
        return `'${x.name}' row ${i + 1}: ${x.getString(i)} vs ${y.getString(i)}`;
    }
  }
  return null;
}

/** Thrown by a snapshot write when the platform's d42 encoding would not give the frame back
 *  exactly (e.g. a null in a bigint column). Nothing is written. */
export class SnapshotFidelityError extends Error {
  constructor(public readonly detail: string) {
    super(`This table cannot be stored exactly (${detail}). Nothing was saved.`);
    this.name = 'SnapshotFidelityError';
  }
}

// ─── Errors ──────────────────────────────────────────────────────────────────────

/** The stored analysis no longer points at the snapshot the caller opened: its data was replaced
 *  (a pipeline re-run, another session). Nothing was written. The caller saves beside (study) or
 *  refuses (standalone), so overrides never land on different points. */
export class SourceReplacedError extends Error {
  constructor(public readonly storedRef: string, public readonly expectedRef: string) {
    super(`The analysis data was replaced while it was open (stored '${storedRef}', open '${expectedRef}').`);
    this.name = 'SourceReplacedError';
  }
}

/** Another session saved CONTENT of this analysis (recipe, confirmed recipe, λz overrides or
 *  review decisions) since this view loaded it (PR-1). Nothing was written; the caller takes the
 *  same path as {@link SourceReplacedError}. Never merged: a merged override set would not match
 *  the numbers computed in this view. */
export class ContentChangedError extends Error {
  constructor(public readonly label: string, public readonly by: string | null) {
    super(`'${label}' was changed in another session${by !== null && by !== '' ? ` (by ${by})` : ''}.`);
    this.name = 'ContentChangedError';
  }
}

/** A write kept losing a version race (≤ 3 attempts). */
export class KeepsChangingError extends Error {
  constructor(label: string) {
    super(`'${label}' kept changing during the save — try again.`);
    this.name = 'KeepsChangingError';
  }
}

// ─── Read models ─────────────────────────────────────────────────────────────────

/** What the caller may do with a study (server-computed `~can_*`). Sharing is View-only until
 *  the platform honours row Edit grants (decision 10). */
export interface StudyAccess {
  readonly canEdit: boolean;
  readonly canDelete: boolean;
  readonly canShare: boolean;
}

/** One analysis as the Hub, the tree and the router need it (light columns only, SD-17). */
export interface AnalysisHeader {
  readonly id: string;
  readonly studyId: string;
  readonly label: string;
  readonly engineVersion: string | null;
  /** The reporting-basis stamp of the cache (`RESULTS_BASIS`); null = written before it. */
  readonly resultsBasis: string | null;
  readonly hasCache: boolean;
  readonly stageCache: StageCache | null;
  readonly pipelineFlagCount: number;
  readonly movedFromLabel: string | null;
  readonly updatedAt: string;
}

/** One study (or standalone container) as the lists need it. */
export interface StudyHeader {
  readonly id: string;
  readonly label: string;
  readonly kind: 'study' | 'standalone';
  readonly updatedAt: string;
  readonly analysisCount: number;
  readonly access: StudyAccess;
  readonly analyses: ReadonlyArray<AnalysisHeader>;
}

/** A study loaded for its overview: every analysis column except the ledgers (SD-21). The
 *  entries' `lzOverrides` / `reviewDecisions` are `[]` here — only {@link NcaRepository.loadAnalysis}
 *  reads them — and their pipeline flags are only counted (`pipelineFlagCount`).
 *
 * Curves may come later (ems-modification step 5): `pendingCurves` maps an analysis id to, per
 * tier row, the id to pass to {@link NcaRepository.loadCurves} for a curve not loaded yet (null where
 * the row has none). Such a row's `curve` is `null` until {@link withCurves} fills it — `curve: null`
 * on a row NOT listed here still means "no curve". Absent = every curve is already in. */
export interface LoadedStudy extends NcaStudy {
  readonly pendingCurves?: ReadonlyMap<string, ReadonlyArray<string | null>>;
}

/** Is the curve of tier row `tierIndex` of analysis `analysisId` still to be loaded? */
export function isCurvePending(study: LoadedStudy, analysisId: string, tierIndex: number): boolean {
  return (study.pendingCurves?.get(analysisId)?.[tierIndex] ?? null) !== null;
}

/** Every curve id {@link LoadedStudy.pendingCurves} lists. */
export function pendingCurveIds(study: LoadedStudy): string[] {
  return [...(study.pendingCurves?.values() ?? [])].flat().filter((id): id is string => id !== null);
}

/** The study with the loaded curves filled in. Nothing stays pending: an id the load did not
 *  return (the row was deleted or recomputed meanwhile) reads as no curve, so a caller that loads
 *  until nothing is pending stops; the next load of the study reads that row afresh. */
export function withCurves(study: LoadedStudy, curves: ReadonlyMap<string, CompoundCurve | null>): LoadedStudy {
  const analyses = study.analyses.map((a) => {
    const ids = study.pendingCurves?.get(a.id);
    if (ids === undefined || a.summaryCache === null) return a;
    return {...a, summaryCache: {...a.summaryCache, compoundTier: a.summaryCache.compoundTier.map((r, i) => {
      const id = ids[i] ?? null;
      return id !== null ? {...r, curve: curves.get(id) ?? null} : r;
    })}};
  });
  return {...study, analyses, pendingCurves: undefined};
}

/** One analysis loaded for an open: the full entry, its study and its source frame. */
export interface LoadedAnalysis {
  readonly study: StudyHeader;
  readonly entry: AnalysisEntry;
  /** The snapshot frame (a fresh copy — the caller may mutate it). */
  readonly df: DG.DataFrame;
}

// ─── Write inputs ────────────────────────────────────────────────────────────────

/** The workflow state that travels with a cache (v15). */
export interface WorkflowState {
  readonly confirmed: ConfirmedRecipe | null;
  readonly lzOverrides: ReadonlyArray<LzOverride>;
  readonly reviewDecisions: ReadonlyArray<ReviewLedger>;
  readonly stageCache: StageCache | null;
}

/** A new analysis to write (Add to study, first save, save beside, pipeline). */
export interface NewAnalysisInput {
  readonly id?: string;
  readonly label: string;
  /** The source frame; snapshotted as d42. */
  readonly df: DG.DataFrame;
  readonly recipe: AnalysisRecipe;
  readonly summaryCache: SummaryCache | null;
  readonly workflow: WorkflowState;
  readonly origin: AnalysisOrigin | null;
  readonly pipelineFlags?: ReadonlyArray<PipelineFlag>;
}

/** An explicit Save of an open analysis (PR-1). */
export interface SaveAnalysisInput {
  readonly df: DG.DataFrame;
  readonly recipe: AnalysisRecipe;
  readonly summaryCache: SummaryCache | null;
  readonly workflow: WorkflowState;
}

/** A clean compute's background cache refresh (SD-18). */
export interface SummaryWrite {
  readonly summaryCache: SummaryCache;
  readonly workflow: WorkflowState;
  /** {@link computeSnapshotHash} of the frame the compute ran on (PR-2): skipped when it differs
   *  from the stored snapshot's hash — the view computed from unsaved data. */
  readonly frameHash: string;
}

/** Where a Move puts a standalone analysis. */
export type MoveTarget =
  | {readonly kind: 'existing'; readonly studyId: string}
  | {readonly kind: 'new'; readonly label: string};

/** What a Move did. */
export interface MoveResult {
  readonly studyId: string;
  readonly studyLabel: string;
  readonly analysisId: string;
  readonly label: string;
  readonly source: SourcePointer;
  readonly passthroughCols: readonly PassthroughColumnSpec[];
}

/** The pipeline's base read (SD-20): its own study (full entries) or none. */
export interface PipelineBase {
  readonly study: NcaStudy | undefined;
}

/** A pipeline save (S5): the planned study (pure planner output) over the base it was planned on.
 *  Its new entries already point at snapshots written with {@link NcaRepository.writeSnapshot}. An
 *  entry whose id is in `base` is updated in place; any other is inserted. */
export interface PipelineWrite {
  readonly base: NcaStudy | undefined;
  readonly planned: NcaStudy;
  readonly owner: string;
  readonly studyKey: string;
}

/** A written snapshot. */
export interface WrittenSnapshot {
  readonly source: SourcePointer;
  readonly signature: SourceSignature;
}

// ─── The contract ────────────────────────────────────────────────────────────────

export interface NcaRepository {
  /** Every visible study and standalone container with its analyses' light columns (SD-17). */
  listStudies(): Promise<StudyHeader[]>;
  /** A study with every analysis column except the ledgers (SD-21). Throws when not visible. */
  loadStudy(studyId: string): Promise<LoadedStudy>;
  /** One analysis, its study and its snapshot frame (C5). Throws when not visible. */
  loadAnalysis(analysisId: string): Promise<LoadedAnalysis>;
  /** The curves of {@link LoadedStudy.pendingCurves} ids, by id (one request; an invisible or
   *  deleted row is absent from the result). */
  loadCurves(ids: readonly string[]): Promise<Map<string, CompoundCurve | null>>;

  createStudy(label: string): Promise<StudyHeader>;
  renameStudy(studyId: string, label: string): Promise<void>;
  /** Final, no restore (user decision 2026-10-21, superseding decisions 11 / 12): the study, its
   *  analyses and caches leave every listing, and their snapshot files are deleted. The rows stay
   *  soft-deleted in the schema only because a package schema has no purge (core gap). */
  deleteStudy(studyId: string): Promise<void>;
  /** View grant on a study (decision 10). Refuses a standalone container. */
  share(studyId: string, groupId: string): Promise<void>;

  /** A new standalone container holding `input` (the first save of an analysis). */
  createStandalone(input: NewAnalysisInput): Promise<{studyId: string; analysisId: string; label: string;
    source: SourcePointer}>;
  /** Append to a study; the label is made unique within it. */
  appendAnalysis(studyId: string, input: NewAnalysisInput): Promise<{analysisId: string; label: string;
    source: SourcePointer}>;
  /** Explicit Save (decision 1, PR-1). `expectedRef` = the source the view opened. Throws
   *  {@link SourceReplacedError} or {@link ContentChangedError}; returns the stored source. */
  saveAnalysis(analysisId: string, expectedRef: string, input: SaveAnalysisInput): Promise<SourcePointer>;
  /** Background cache refresh of a clean view (SD-18, PR-2). Silent no-op on any conflict. */
  writeSummary(analysisId: string, expectedRef: string, write: SummaryWrite): Promise<void>;
  /** Rename one analysis (a standalone container's label follows it, in one transaction). */
  renameAnalysis(analysisId: string, label: string): Promise<void>;
  /** Final, like {@link deleteStudy}; a standalone analysis takes its container with it. */
  deleteAnalysis(analysisId: string): Promise<void>;
  /** Move a standalone analysis into a study (AC-3.5, SD-14). */
  moveToStudy(standaloneId: string, target: MoveTarget, opts?: {
    label?: string; passthroughCols?: readonly PassthroughColumnSpec[]; summaryCache?: SummaryCache;
  }): Promise<MoveResult>;

  /** The pipeline's own study by business key (SD-20), full entries. */
  pipelineBase(owner: string, studyKey: string): Promise<PipelineBase>;
  /** Write a planned pipeline study in one transaction (S5); returns the study id. Throws a
   *  {@link KeepsChangingError} when the base moved (the caller re-plans). */
  savePipelinePlan(write: PipelineWrite): Promise<string>;

  /** Store a frame as a new snapshot (verify-on-write, BD-3). Throws {@link SnapshotFidelityError}. */
  writeSnapshot(df: DG.DataFrame, computedAt: string): Promise<WrittenSnapshot>;
  /** Remove a snapshot no row points at any more (a superseded or unused one). Never throws. */
  deleteSnapshot(source: SourcePointer): Promise<void>;
  /** The snapshot frame a source pointer names (decision 3, SD-4). */
  readSource(source: SourcePointer): Promise<DG.DataFrame>;
}

// ─── The session singleton (spec principle 1) ───────────────────────────────────

let current: NcaRepository | null = null;
let factory: (() => NcaRepository) | null = null;

/** Register how the production repository is built (package init). Kept as a factory so the
 *  adapter does not import the EMS client at module init. */
export function registerNcaRepositoryFactory(f: () => NcaRepository): void {
  factory = f;
}

/** The session's repository. */
export function getNcaRepository(): NcaRepository {
  if (current === null) {
    if (factory === null) throw new Error('NCA storage is not initialised');
    current = factory();
  }
  return current;
}

/** Tests: swap in another repository (an in-memory one); returns the previous one to restore. */
export function setNcaRepository(repo: NcaRepository | null): NcaRepository | null {
  const prev = current;
  current = repo;
  return prev;
}
