/**
 * An in-memory {@link NcaRepository} for grid / workflow / pipeline-planner tests (nca-ems-storage-
 * migration decision 7). It keeps the EMS semantics the callers depend on — row versions and
 * version conflicts, final delete (rows hidden, snapshots removed), label rules,
 * verify-on-write of snapshots — and the shared contract suite (`nca-repository-contract-tests.ts`)
 * runs against it and the EMS repository alike, so it cannot drift (AC-3.8).
 *
 * Test hooks: {@link MemoryNcaRepository.bumpContent} / `replaceSource` stage another session's write.
 */

import * as DG from 'datagrok-api/dg';

import {
  newId, uniqueLabel, computeSnapshotHash, snapshotSignature, snapshotDiff, snapshotFrame, SourceReplacedError,
  ContentChangedError, KeepsChangingError, SnapshotFidelityError,
  STANDALONE_ROUTE_SEGMENT, NCA_STUDIES_CHANGED_EVENT,
} from '../adapter/nca-repository';
import type {
  NcaRepository, StudyHeader, LoadedStudy, LoadedAnalysis, NewAnalysisInput,
  SaveAnalysisInput, SummaryWrite, MoveTarget, MoveResult, PipelineBase, PipelineWrite, WrittenSnapshot,
} from '../adapter/nca-repository';
import type {AnalysisEntry, CompoundCurve, NcaStudy, SourcePointer, SummaryCache} from '../adapter';
import {INSTALLED_ENGINE_VERSION, RESULTS_BASIS} from '../adapter';
import * as grok from 'datagrok-api/grok';

interface MemStudy {
  id: string; label: string; kind: 'study' | 'standalone'; createdAt: string; updatedAt: string; version: number;
  deleted: boolean; movedTo: string | null;
  pipelineKey?: {owner: string; studyKey: string};
}
interface MemAnalysis {
  entry: AnalysisEntry; studyId: string; version: number; deleted: boolean;
  movedFromLabel: string | null; updatedAt: string;
}

const clone = <T>(x: T): T => JSON.parse(JSON.stringify(x)) as T;
const nowIso = (): string => new Date().toISOString();
/** The entry with `resultsComputedAt` moved iff its results differ from `before`'s (decision 2: a
 *  write that recomputes the same numbers does not move it — the EMS repository's rule). */
function stampResults(before: AnalysisEntry | null, after: AnalysisEntry): AnalysisEntry {
  const results = (e: AnalysisEntry | null): string => JSON.stringify(e === null ? null :
    [e.summaryCache, e.engineVersion ?? null, e.resultsBasis ?? null]);
  if (after.summaryCache === null) return {...after, resultsComputedAt: null};
  return results(before) === results(after) ? after : {...after, resultsComputedAt: nowIso()};
}

const content = (e: AnalysisEntry): string => JSON.stringify({
  r: {columns: e.columns, modes: e.modes, ncaRules: e.ncaRules, units: e.units, groupByCols: e.groupByCols,
    passthroughCols: e.passthroughCols, acceptanceRSquared: e.acceptanceRSquared},
  c: e.confirmed ?? null, l: e.lzOverrides ?? [], d: e.reviewDecisions ?? []});

export class MemoryNcaRepository implements NcaRepository {
  studies = new Map<string, MemStudy>();
  analyses = new Map<string, MemAnalysis>();
  snapshots = new Map<string, Uint8Array>();

  /** Another session on the SAME store: shares the rows and snapshots, keeps its own tracker. */
  peer(): MemoryNcaRepository {
    const p = new MemoryNcaRepository(this.fireEvents);
    p.studies = this.studies;
    p.analyses = this.analyses;
    p.snapshots = this.snapshots;
    return p;
  }
  /** Test hook: runs at the start of every `savePipelinePlan` with its 1-based number — where a
   *  test stages "another writer moved the base between the plan's read and its write". */
  beforePipelineWrite: (n: number) => void = () => {};
  private pipelineWriteCount = 0;
  /** How many pipeline plan writes were attempted. */
  get pipelineWrites(): number {return this.pipelineWriteCount;}
  /** Per analysis: the version + content this "session" last saw (the EMS tracker's twin). */
  private readonly seen = new Map<string, {version: number; content: string; ref: string}>();
  /** Fire the platform change event (off in pure unit tests). */
  constructor(private readonly fireEvents = false) {}

  private fire(): void {
    if (this.fireEvents) grok.events.fireCustomEvent(NCA_STUDIES_CHANGED_EVENT, {});
  }

  private see(a: MemAnalysis): void {
    this.seen.set(a.entry.id, {version: a.version, content: content(a.entry), ref: a.entry.source.ref});
  }

  // ── test hooks ─────────────────────────────────────────────────────────────

  /** Another writer adds an analysis to `studyId` — synchronously, so a test hook can land it
   *  between a read and a write (a copy of the study's first analysis under a new id / label). */
  addForeignAnalysis(studyId: string, label: string): void {
    const src = [...this.analyses.values()].find((a) => a.studyId === studyId && !a.deleted)!;
    const id = `foreign-${newId()}`;
    this.analyses.set(id, {entry: {...clone(src.entry), id, label, origin: null},
      studyId, version: 1, deleted: false, movedFromLabel: null, updatedAt: nowIso()});
    this.studies.get(studyId)!.version++;
  }

  /** Put a prebuilt study (a grid / overview fixture) into the store as is. */
  seedStudy(study: NcaStudy): void {
    this.studies.set(study.id, {id: study.id, label: study.label, kind: study.kind, createdAt: study.createdAt,
      updatedAt: study.updatedAt, version: 1, deleted: false, movedTo: null,
      ...(study.pipelineKey !== undefined ? {pipelineKey: study.pipelineKey} : {})});
    for (const e of study.analyses) {
      this.analyses.set(e.id, {entry: clone(e), studyId: study.id, version: 1, deleted: false,
        movedFromLabel: null, updatedAt: study.updatedAt});
    }
  }

  /** Another session saves content of `analysisId` (PR-1). */
  bumpContent(analysisId: string, patch: Partial<AnalysisEntry>): void {
    const a = this.analyses.get(analysisId)!;
    a.entry = {...a.entry, ...patch};
    a.version++;
  }

  /** Another writer replaces the data of `analysisId` (a pipeline re-run, F11). */
  replaceSource(analysisId: string, df: DG.DataFrame): void {
    const a = this.analyses.get(analysisId)!;
    const path = `mem/${newId()}/source.d42`;
    this.snapshots.set(path, df.toByteArray());
    a.entry = {...a.entry, source: {kind: 'stored', ref: path}, sourceSignature: snapshotSignature(df, nowIso())};
    a.version++;
  }

  // ── reads ──────────────────────────────────────────────────────────────────

  private header(s: MemStudy): StudyHeader {
    const list = [...this.analyses.values()].filter((a) => a.studyId === s.id && !a.deleted)
      .sort((x, y) => x.entry.label.localeCompare(y.entry.label));
    return {
      id: s.id, label: s.label, kind: s.kind, updatedAt: s.updatedAt, analysisCount: list.length,
      access: {canEdit: true, canDelete: true, canShare: s.kind === 'study'},
      analyses: list.map((a) => ({
        id: a.entry.id, studyId: s.id, label: a.entry.label, engineVersion: a.entry.engineVersion ?? null,
        resultsBasis: a.entry.resultsBasis ?? null, hasCache: a.entry.summaryCache !== null,
        stageCache: a.entry.stageCache ?? null,
        pipelineFlagCount: (a.entry.pipelineFlags ?? []).filter((f) => f.pipelineOnly.length > 0).length,
        movedFromLabel: a.movedFromLabel,
        updatedAt: a.updatedAt,
      })),
    };
  }

  async listStudies(): Promise<StudyHeader[]> {
    return [...this.studies.values()].filter((s) => !s.deleted).sort((a, b) => a.label.localeCompare(b.label))
      .map((s) => this.header(s));
  }

  private studyOf(s: MemStudy, ledgers: boolean): NcaStudy {
    const entries = [...this.analyses.values()].filter((a) => a.studyId === s.id && !a.deleted)
      .sort((x, y) => x.entry.label.localeCompare(y.entry.label))
      .map((a) => {
        this.see(a);
        const e = clone(a.entry);
        return ledgers ? e : {...e, lzOverrides: [], reviewDecisions: []};
      });
    return {id: s.id, label: s.label, kind: s.kind, createdAt: s.createdAt, updatedAt: s.updatedAt,
      access: {canEdit: true, canDelete: true, canShare: s.kind === 'study'}, analyses: entries,
      ...(s.pipelineKey !== undefined ? {pipelineKey: s.pipelineKey} : {})};
  }

  async loadStudy(studyId: string): Promise<LoadedStudy> {
    const s = this.studies.get(studyId);
    if (s === undefined || s.deleted) throw new Error(`Study '${studyId}' is not available.`);
    return this.studyOf(s, false);
  }

  async loadAnalysis(analysisId: string): Promise<LoadedAnalysis> {
    const a = this.analyses.get(analysisId);
    if (a === undefined || a.deleted) throw new Error(`Analysis '${analysisId}' is not available.`);
    const s = this.studies.get(a.studyId)!;
    this.see(a);
    return {study: {...this.header(s), access: {canEdit: true, canDelete: true, canShare: false}, analyses: []},
      entry: clone(a.entry), df: await this.readSource(a.entry.source)};
  }

  /** The memory store returns every curve with the study (no `pendingCurves`), so nothing is pending. */
  async loadCurves(_ids: readonly string[]): Promise<Map<string, CompoundCurve | null>> {
    return new Map();
  }

  // ── studies ────────────────────────────────────────────────────────────────

  async createStudy(label: string): Promise<StudyHeader> {
    const s: MemStudy = {id: newId(), label, kind: 'study', createdAt: nowIso(), updatedAt: nowIso(), version: 1,
      deleted: false, movedTo: null};
    this.studies.set(s.id, s);
    this.fire();
    return this.header(s);
  }

  async renameStudy(studyId: string, label: string): Promise<void> {
    const s = this.studies.get(studyId)!;
    s.label = label; s.updatedAt = nowIso(); s.version++;
    this.fire();
  }

  /** Final delete: the study and its analyses leave every listing; their snapshots go. */
  private remove(studyId: string): void {
    const s = this.studies.get(studyId);
    if (s === undefined || s.deleted) return;
    for (const a of this.analyses.values()) {
      if (a.studyId === studyId && !a.deleted) {
        a.deleted = true; a.version++;
        this.snapshots.delete(a.entry.source.ref);
      }
    }
    s.deleted = true; s.updatedAt = nowIso(); s.version++;
  }

  async deleteStudy(studyId: string): Promise<void> {
    this.remove(studyId);
    this.fire();
  }

  async share(studyId: string, _groupId: string): Promise<void> {
    const s = this.studies.get(studyId);
    if (s === undefined) throw new Error('The study is not available.');
    if (s.kind !== 'study') throw new Error('Only a study can be shared; move the analysis into a study first.');
  }

  // ── analyses ───────────────────────────────────────────────────────────────

  private entryFrom(id: string, label: string, input: NewAnalysisInput, snap: WrittenSnapshot): AnalysisEntry {
    return {
      ...input.recipe, id, label, source: snap.source, sourceSignature: snap.signature,
      summaryCache: input.summaryCache, engineVersion: input.summaryCache === null ? null : INSTALLED_ENGINE_VERSION,
      resultsBasis: input.summaryCache === null ? null : RESULTS_BASIS,
      confirmed: input.workflow.confirmed, lzOverrides: input.workflow.lzOverrides,
      reviewDecisions: input.workflow.reviewDecisions, stageCache: input.workflow.stageCache,
      origin: input.origin, pipelineFlags: input.pipelineFlags ?? [],
    };
  }

  private add(studyId: string, entry: AnalysisEntry): MemAnalysis {
    const a: MemAnalysis = {entry: stampResults(null, entry), studyId, version: 1, deleted: false, movedFromLabel: null,
      updatedAt: nowIso()};
    this.analyses.set(entry.id, a);
    this.see(a);
    return a;
  }

  async createStandalone(input: NewAnalysisInput): Promise<{studyId: string; analysisId: string; label: string;
    source: SourcePointer}> {
    const taken = new Set([...this.studies.values()].filter((s) => !s.deleted && s.kind === 'standalone')
      .map((s) => s.label));
    const label = uniqueLabel(input.label.trim() || 'Analysis', taken);
    const snap = await this.writeSnapshot(input.df, nowIso());
    const s: MemStudy = {id: newId(), label, kind: 'standalone', createdAt: nowIso(), updatedAt: nowIso(), version: 1,
      deleted: false, movedTo: null};
    this.studies.set(s.id, s);
    const id = input.id ?? newId();
    this.add(s.id, this.entryFrom(id, label, input, snap));
    this.fire();
    return {studyId: s.id, analysisId: id, label, source: snap.source};
  }

  async appendAnalysis(studyId: string, input: NewAnalysisInput): Promise<{analysisId: string; label: string;
    source: SourcePointer}> {
    const taken = new Set([...this.analyses.values()].filter((a) => a.studyId === studyId && !a.deleted)
      .map((a) => a.entry.label));
    const label = uniqueLabel(input.label.trim() || 'Analysis', taken);
    const snap = await this.writeSnapshot(input.df, nowIso());
    const id = input.id ?? newId();
    this.add(studyId, this.entryFrom(id, label, input, snap));
    const s = this.studies.get(studyId)!;
    s.updatedAt = nowIso(); s.version++;
    this.fire();
    return {analysisId: id, label, source: snap.source};
  }

  async saveAnalysis(analysisId: string, expectedRef: string, input: SaveAnalysisInput): Promise<SourcePointer> {
    const a = this.analyses.get(analysisId);
    if (a === undefined || a.deleted) throw new Error('The analysis no longer exists.');
    const seen = this.seen.get(analysisId) ?? {version: a.version, content: content(a.entry), ref: a.entry.source.ref};
    if (a.entry.source.ref !== expectedRef) throw new SourceReplacedError(a.entry.source.ref, expectedRef);
    if (a.version !== seen.version && content(a.entry) !== seen.content)
      throw new ContentChangedError(a.entry.label, null);
    const hash = computeSnapshotHash(input.df);
    const snap = hash === a.entry.sourceSignature.hash ? null : await this.writeSnapshot(input.df, nowIso());
    if (snap !== null) this.snapshots.delete(a.entry.source.ref);
    a.entry = stampResults(a.entry, {
      ...a.entry, ...input.recipe,
      ...(snap !== null ? {source: snap.source, sourceSignature: snap.signature} : {}),
      summaryCache: input.summaryCache, engineVersion: input.summaryCache === null ? null : INSTALLED_ENGINE_VERSION,
      resultsBasis: input.summaryCache === null ? null : RESULTS_BASIS,
      confirmed: input.workflow.confirmed, lzOverrides: input.workflow.lzOverrides,
      reviewDecisions: input.workflow.reviewDecisions, stageCache: input.workflow.stageCache,
    });
    a.version++; a.updatedAt = nowIso();
    this.see(a);
    this.fire();
    return a.entry.source;
  }

  async writeSummary(analysisId: string, expectedRef: string, write: SummaryWrite): Promise<void> {
    const a = this.analyses.get(analysisId);
    const seen = this.seen.get(analysisId);
    if (a === undefined || a.deleted || seen === undefined) return;
    if (a.entry.source.ref !== expectedRef || a.version !== seen.version) return;
    if (a.entry.sourceSignature.hash !== null && write.frameHash !== a.entry.sourceSignature.hash) return;
    const next = stampResults(a.entry, {...a.entry, summaryCache: write.summaryCache,
      engineVersion: INSTALLED_ENGINE_VERSION, resultsBasis: RESULTS_BASIS, confirmed: write.workflow.confirmed,
      lzOverrides: write.workflow.lzOverrides, reviewDecisions: write.workflow.reviewDecisions,
      stageCache: write.workflow.stageCache});
    // Diff-gated like the EMS write: nothing differs → nothing is written, the version stays.
    if (JSON.stringify(next) === JSON.stringify(a.entry)) return;
    a.entry = next;
    a.version++;
    this.see(a);
  }

  async renameAnalysis(analysisId: string, label: string): Promise<void> {
    const a = this.analyses.get(analysisId)!;
    a.entry = {...a.entry, label}; a.version++;
    const s = this.studies.get(a.studyId)!;
    if (s.kind === 'standalone') s.label = label;
    s.updatedAt = nowIso(); s.version++;
    this.fire();
  }

  async deleteAnalysis(analysisId: string): Promise<void> {
    const a = this.analyses.get(analysisId);
    if (a === undefined || a.deleted) return;
    const s = this.studies.get(a.studyId)!;
    if (s.kind === 'standalone') this.remove(s.id);
    else {
      a.deleted = true; a.version++; s.updatedAt = nowIso(); s.version++;
      this.snapshots.delete(a.entry.source.ref);
    }
    this.fire();
  }

  async moveToStudy(standaloneId: string, target: MoveTarget, opts: {
    label?: string; passthroughCols?: MoveResult['passthroughCols']; summaryCache?: SummaryCache;
  } = {}): Promise<MoveResult> {
    const c = this.studies.get(standaloneId);
    if (c === undefined || c.deleted || c.kind !== 'standalone')
      throw new Error('The analysis is no longer in Ungrouped.');
    const a = [...this.analyses.values()].find((x) => x.studyId === standaloneId && !x.deleted);
    if (a === undefined) throw new Error('The standalone analysis is empty.');
    const summaryCache = opts.summaryCache ?? a.entry.summaryCache;
    if (summaryCache === null)
      throw new Error('Confirm the setup first: an analysis without results cannot join a study.');
    let targetStudy: MemStudy;
    if (target.kind === 'existing') targetStudy = this.studies.get(target.studyId)!;
    else {
      const labels = new Set([...this.studies.values()].filter((s) => !s.deleted && s.kind === 'study')
        .map((s) => s.label));
      targetStudy = {id: newId(), label: uniqueLabel(target.label.trim() || 'New study',
        new Set([...labels, STANDALONE_ROUTE_SEGMENT])), kind: 'study', createdAt: nowIso(), updatedAt: nowIso(),
      version: 1, deleted: false, movedTo: null};
      this.studies.set(targetStudy.id, targetStudy);
    }
    const taken = new Set([...this.analyses.values()].filter((x) => x.studyId === targetStudy.id && !x.deleted)
      .map((x) => x.entry.label));
    const label = uniqueLabel((opts.label ?? a.entry.label).trim() || a.entry.label, taken);
    const passthroughCols = opts.passthroughCols ?? a.entry.passthroughCols;
    a.entry = stampResults(a.entry, {...a.entry, label, passthroughCols, summaryCache,
      // A move with a fresh cache stamps it, as the EMS move does.
      ...(opts.summaryCache !== undefined ?
        {engineVersion: INSTALLED_ENGINE_VERSION, resultsBasis: RESULTS_BASIS} : {}),
      origin: {kind: 'standalone', label: c.label, containerId: c.id}});
    a.studyId = targetStudy.id; a.movedFromLabel = c.label; a.version++;
    targetStudy.updatedAt = nowIso(); targetStudy.version++;
    c.deleted = true; c.movedTo = targetStudy.id; c.version++;
    this.fire();
    return {studyId: targetStudy.id, studyLabel: targetStudy.label, analysisId: a.entry.id, label,
      source: a.entry.source, passthroughCols};
  }

  // ── pipeline ───────────────────────────────────────────────────────────────

  /** The study version each pipeline base read saw, keyed by the base it returned (the EMS
   *  repository's twin) — so a later read never refreshes an earlier base. */
  private readonly pipelineSeen = new WeakMap<NcaStudy, number>();

  async pipelineBase(owner: string, studyKey: string): Promise<PipelineBase> {
    const s = [...this.studies.values()].find((x) => !x.deleted && x.pipelineKey?.owner === owner &&
      x.pipelineKey?.studyKey === studyKey);
    if (s === undefined) return {study: undefined};
    const study = this.studyOf(s, true);
    this.pipelineSeen.set(study, s.version);
    return {study};
  }

  async savePipelinePlan(write: PipelineWrite): Promise<string> {
    this.beforePipelineWrite(++this.pipelineWriteCount);
    let s: MemStudy;
    if (write.base === undefined) {
      if ([...this.studies.values()].some((x) => !x.deleted && x.pipelineKey?.owner === write.owner &&
        x.pipelineKey?.studyKey === write.studyKey)) throw new KeepsChangingError(write.planned.label);
      s = {id: newId(), label: write.planned.label, kind: 'study', createdAt: nowIso(), updatedAt: nowIso(), version: 1,
        deleted: false, movedTo: null, pipelineKey: {owner: write.owner, studyKey: write.studyKey}};
      this.studies.set(s.id, s);
    } else {
      s = this.studies.get(write.base.id)!;
      // A plan on a base that moved since its read (another write to the study) re-plans.
      if (this.pipelineSeen.get(write.base) !== s.version) throw new KeepsChangingError(write.planned.label);
      for (const b of write.base.analyses) {
        const live = this.analyses.get(b.id);
        if (live === undefined || live.deleted || JSON.stringify(live.entry) !== JSON.stringify(b))
          throw new KeepsChangingError(write.planned.label);
      }
    }
    for (const e of write.planned.analyses) {
      const live = this.analyses.get(e.id);
      if (live !== undefined) {
        live.entry = stampResults(live.entry, clone(e));
        live.version++;
      } else this.add(s.id, clone(e));
    }
    s.updatedAt = nowIso(); s.version++;
    this.fire();
    return s.id;
  }

  // ── snapshots ──────────────────────────────────────────────────────────────

  async writeSnapshot(frame: DG.DataFrame, computedAt: string): Promise<WrittenSnapshot> {
    const df = snapshotFrame(frame); // PR-13: the source columns only
    const bytes = df.toByteArray();
    const diff = snapshotDiff(df, DG.DataFrame.fromByteArray(bytes));
    if (diff !== null) throw new SnapshotFidelityError(diff);
    const path = `mem/${newId()}/source.d42`;
    this.snapshots.set(path, bytes);
    return {source: {kind: 'stored', ref: path}, signature: snapshotSignature(df, computedAt)};
  }

  async deleteSnapshot(source: SourcePointer): Promise<void> {
    this.snapshots.delete(source.ref);
  }

  async readSource(source: SourcePointer): Promise<DG.DataFrame> {
    const bytes = this.snapshots.get(source.ref);
    if (bytes === undefined) throw new Error(`no snapshot '${source.ref}'`);
    return DG.DataFrame.fromByteArray(bytes);
  }
}

/** Every study's entries (tests). */
export function allEntries(repo: MemoryNcaRepository): AnalysisEntry[] {
  return [...repo.analyses.values()].filter((a) => !a.deleted).map((a) => a.entry);
}
