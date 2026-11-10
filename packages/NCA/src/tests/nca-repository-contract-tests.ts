/**
 * NCA storage contract — ONE suite, run against the in-memory repository and the live EMS one,
 * so the fake the grid / workflow tests use cannot drift from the store the app ships with.
 *
 * Each implementation supplies two hooks that stage ANOTHER session's write: a content change
 * (PR-1) and a data replacement (F11). The EMS run writes real `nca` rows (labels `zz-contract …`)
 * and soft-deletes every study it made in `finally`.
 */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';

import {
  ContentChangedError, SourceReplacedError, KeepsChangingError, SnapshotFidelityError,
  computeSnapshotHash, snapshotDiff, newId, isStale, isEntryStale,
} from '../adapter/nca-repository';
import {INSTALLED_ENGINE_VERSION} from '../adapter/engine-version';
import {RESULTS_BASIS} from '../adapter/unit-conversion';
import type {NcaRepository, NewAnalysisInput, WorkflowState} from '../adapter/nca-repository';
import {EmsNcaRepository, isConflict, isDefinitiveRefusal} from '../adapter/ems-nca-repository';
import {recipeOf} from '../adapter/study-schema';
import {decisionRowsOf, normLedgers, normLz, isoOf} from '../adapter/ems-nca-rows';
import {pendingCurveIds, withCurves} from '../adapter/nca-repository';
import {computeRecipeFingerprint} from '../adapter/fingerprint';
import type {AnalysisRecipe, CompoundTierRow, NcaStudy, PipelineOrigin, SummaryCache} from '../adapter';
import {ncaDb} from '../generated/db';
import {MemoryNcaRepository} from './memory-nca-repository';
import {withRestrictedUser} from './restricted-user';
import GOLDENS from './goldens/ems-pre-migration.json';

const TIMEOUT = 120000;

/** How an implementation stages another session's writes. */
interface Harness {
  readonly name: string;
  make(): NcaRepository;
  /** A second session on the same store as `repo` (its own caches and tracker). */
  peer(repo: NcaRepository): NcaRepository;
  /** Another session saves content (a review decision) of `analysisId`. */
  foreignContent(repo: NcaRepository, analysisId: string): Promise<void>;
  /** Another session renames `analysisId` (a non-content write). */
  foreignRename(repo: NcaRepository, analysisId: string): Promise<void>;
  /** A pipeline re-run replaces the data of `analysisId`. */
  foreignReplace(repo: NcaRepository, analysisId: string): Promise<void>;
}

const RECIPE = {
  columns: {subject: 'SUBJ', time: 'TIME', concentration: 'CONC'}, modes: {}, ncaRules: {}, units: {},
  groupByCols: [], passthroughCols: [], acceptanceRSquared: 0.9,
} as unknown as AnalysisRecipe;
const CACHE = {analysisTier: {nCompounds: 1, nSubjects: 2, route: 'PO', dose: '1', status: 'ok',
  flagSummary: ''}, compoundTier: []} as unknown as SummaryCache;
const WF: WorkflowState = {confirmed: null, lzOverrides: [], reviewDecisions: [],
  stageCache: {stage: 'complete', pending: 0, notApplicable: 0}};

/** A summary with `n` grouped tier rows, each with a curve (the lazily loaded part). */
function tierCache(n: number): SummaryCache {
  const mc = {mean: 1.5, cvPct: 10, n: 2};
  const row = (i: number): CompoundTierRow => ({compoundId: `C${i}`, isWholeAnalysis: false,
    groupValues: [{key: 'Compound', value: `C${i}`}], keyInventory: [], passthrough: [], n: 2, route: 'PO',
    doseNorm: 'absolute', cmax: mc, tmax: mc, aucInf: mc, thalf: mc, cl: mc, vz: mc, mrt: mc,
    vss: {mean: NaN, cvPct: NaN, n: 0}, fPct: null, pctExtrap: mc,
    curve: {series: [{group: `C${i}`, time: [0, 1, 2], geoMean: [0, 5, 2.5 + i], geoLower: [0, 4, 2],
      geoUpper: [0, 6, 3]}]},
    warnings: []});
  return {analysisTier: {nCompounds: n, nSubjects: 2 * n, route: 'PO', dose: '1', status: 'ok', flagSummary: ''},
    compoundTier: Array.from({length: n}, (_v, i) => row(i))};
}

const frame = (csv = 'SUBJ,TIME,CONC\nA,0,0\nA,1,5\nA,2,3\nB,0,0\nB,1,6\nB,2,2\n'): DG.DataFrame =>
  DG.DataFrame.fromCsv(csv);

function input(label: string, df = frame(), summaryCache: SummaryCache | null = CACHE): NewAnalysisInput {
  return {label, df, recipe: RECIPE, summaryCache, workflow: WF, origin: null};
}

const DECISION = [{unitKey: 'A', current: {action: 'accept' as const, reasons: [{code: 'X', text: 'x'}], user: 'u',
  at: '2026-01-01T00:00:00Z', fingerprint: 'f'}, history: []}];

const HARNESSES: Harness[] = [
  {
    name: 'memory',
    make: () => new MemoryNcaRepository(),
    peer: (repo) => (repo as MemoryNcaRepository).peer(),
    foreignContent: async (repo, id) => {(repo as MemoryNcaRepository).bumpContent(id, {reviewDecisions: DECISION});},
    foreignRename: async (repo, id) => {
      const m = repo as MemoryNcaRepository;
      const a = m.analyses.get(id)!;
      a.version++;
    },
    foreignReplace: async (repo, id) => {(repo as MemoryNcaRepository).replaceSource(id, frame());},
  },
  {
    name: 'EMS',
    make: () => new EmsNcaRepository(),
    peer: () => new EmsNcaRepository(),
    // A raw write is exactly what another session's repository would send: the decision rows and
    // the guarded analysis update that carries them (decision 14).
    foreignContent: async (_repo, id) => {
      const r = await ncaDb.analyses.get(id);
      await grok.dapi.domains.transaction('nca', [
        {op: 'update', table: 'analysis', id, values: {updated_at: new Date().toISOString()}},
        ...decisionRowsOf(id, r.study_id, DECISION).map((values) => ({op: 'insert', table: 'review_decision', values})),
      ] as never);
    },
    foreignRename: async (_repo, id) => {
      const r = await ncaDb.analyses.get(id);
      await ncaDb.analyses.update(id, {label: `${r.label}`, updated_at: new Date().toISOString()} as never);
    },
    // A pipeline run re-points the row and removes the snapshot it replaced (else the file would
    // outlive every row: a test-made orphan).
    foreignReplace: async (repo, id) => {
      const before = await ncaDb.analyses.get(id);
      await ncaDb.analyses.update(id,
        {source_file: `file://System:DomainFiles/nca/analysis/${newId()}/source.d42`} as never);
      await repo.deleteSnapshot({kind: 'stored', ref: String(before.source_file).replace(/^file:\/\//, '')});
    },
  },
];

async function thrown(p: Promise<unknown>): Promise<unknown> {
  try {await p;} catch (e) {return e;}
  return null;
}

for (const h of HARNESSES) {
  category(`NCA storage contract: ${h.name}`, () => {
    /** Run `body` on a fresh repository; soft-delete every study it made. */
    const withRepo = async (body: (repo: NcaRepository, made: string[]) => Promise<string | void>):
      Promise<string | void> => {
      const repo = h.make();
      const made: string[] = [];
      try {
        return await body(repo, made);
      } finally {
        for (const id of made) try {await repo.deleteStudy(id);} catch {/* gone */}
      }
    };
    const tag = (s: string): string => `zz-contract ${s} ${newId().slice(0, 6)}`;

    test('create, append, list, load: the frame and the entry come back', async () => withRepo(async (repo, made) => {
      const s = await repo.createStudy(tag('s'));
      made.push(s.id);
      const df = frame();
      const a = await repo.appendAnalysis(s.id, input('A'));
      const b = await repo.appendAnalysis(s.id, input('A'));
      expect(b.label, 'A-2', 'labels are unique within a study');
      const header = (await repo.listStudies()).find((x) => x.id === s.id)!;
      expect(header.analysisCount, 2);
      expect(header.analyses.every((x) => x.hasCache && x.stageCache?.stage === 'complete'), true);
      const study = await repo.loadStudy(s.id);
      expect(study.analyses.length, 2);
      expect(study.analyses[0].summaryCache !== null, true, 'the overview carries the cache');
      const opened = await repo.loadAnalysis(a.analysisId);
      expect(snapshotDiff(df, opened.df), null, 'the snapshot is the saved frame');
      expect(opened.entry.columns.subject, 'SUBJ', 'the recipe round-trips');
      expect(opened.study.id, s.id);
    }), {timeout: TIMEOUT});

    test('standalone: unique among standalone labels; one container per analysis', async () =>
      withRepo(async (repo, made) => {
        const label = tag('u');
        const x = await repo.createStandalone(input(label));
        const y = await repo.createStandalone(input(label));
        made.push(x.studyId, y.studyId);
        expect(y.label, `${label}-2`);
        const headers = (await repo.listStudies()).filter((s) => s.id === x.studyId || s.id === y.studyId);
        expect(headers.every((s) => s.kind === 'standalone' && s.analysisCount === 1), true);
      }), {timeout: TIMEOUT});

    test('Save: same data keeps the snapshot; changed data writes a new one', async () =>
      withRepo(async (repo, made) => {
        const s = await repo.createStudy(tag('save'));
        made.push(s.id);
        const a = await repo.appendAnalysis(s.id, input('A'));
        await repo.loadAnalysis(a.analysisId);
        const same = await repo.saveAnalysis(a.analysisId, a.source.ref,
          {df: frame(), recipe: RECIPE, summaryCache: CACHE, workflow: WF});
        expect(same.ref, a.source.ref, 'unchanged data → the same snapshot');
        const changed = frame('SUBJ,TIME,CONC\nA,0,0\nA,1,9\nA,2,3\n');
        const next = await repo.saveAnalysis(a.analysisId, a.source.ref,
          {df: changed, recipe: RECIPE, summaryCache: CACHE, workflow: WF});
        expect(next.ref !== a.source.ref, true, 'changed data → a new snapshot');
        expect(snapshotDiff(changed, await repo.readSource(next)), null);
        expect(await thrown(repo.readSource(a.source)) !== null, true, 'the replaced snapshot is gone');
      }), {timeout: TIMEOUT});

    test('PR-1: content saved by another session → ContentChangedError, nothing written', async () =>
      withRepo(async (repo, made) => {
        const s = await repo.createStudy(tag('pr1'));
        made.push(s.id);
        const a = await repo.appendAnalysis(s.id, input('A'));
        await repo.loadAnalysis(a.analysisId);
        await h.foreignContent(repo, a.analysisId);
        const e = await thrown(repo.saveAnalysis(a.analysisId, a.source.ref,
          {df: frame(), recipe: RECIPE, summaryCache: CACHE, workflow: WF}));
        expect(e instanceof ContentChangedError, true, `expected ContentChangedError, got ${e}`);
        const kept = await h.peer(repo).loadAnalysis(a.analysisId);
        expect((kept.entry.reviewDecisions ?? []).length, 1, 'the other session\'s decision is kept');
      }), {timeout: TIMEOUT});

    test('PR-1 after an own rename: the guard still holds, and a plain Save still lands', async () =>
      withRepo(async (repo, made) => {
        const s = await repo.createStudy(tag('rename'));
        made.push(s.id);
        const a = await repo.appendAnalysis(s.id, input('A'));
        const b = await repo.appendAnalysis(s.id, input('B'));
        await repo.loadAnalysis(a.analysisId);
        await repo.loadAnalysis(b.analysisId);
        await repo.renameAnalysis(a.analysisId, 'A renamed');
        await repo.renameAnalysis(b.analysisId, 'B renamed');
        // A: another session changes the content after the rename — the open's base must still refuse.
        await h.foreignContent(repo, a.analysisId);
        const e = await thrown(repo.saveAnalysis(a.analysisId, a.source.ref,
          {df: frame(), recipe: RECIPE, summaryCache: CACHE, workflow: WF}));
        expect(e instanceof ContentChangedError, true, `expected ContentChangedError, got ${e}`);
        // B: only this session's own rename moved the version — the Save lands.
        await repo.saveAnalysis(b.analysisId, b.source.ref,
          {df: frame(), recipe: RECIPE, summaryCache: CACHE, workflow: {...WF, reviewDecisions: DECISION}});
        const kept = await h.peer(repo).loadAnalysis(b.analysisId);
        expect((kept.entry.reviewDecisions ?? []).length, 1, 'the Save after an own rename landed');
        expect(kept.entry.label, 'B renamed');
      }), {timeout: TIMEOUT});

    test('SD-19: an edit in the middle of a 40-row frame survives Save → reopen', async () =>
      withRepo(async (repo, made) => {
        const s = await repo.createStudy(tag('mid'));
        made.push(s.id);
        const rows = ['SUBJ,TIME,CONC'];
        for (let i = 0; i < 40; i++) rows.push(`S${Math.floor(i / 8)},${i % 8},${(10 - (i % 8)) * 1.5}`);
        const df = frame(rows.join('\n'));
        const a = await repo.appendAnalysis(s.id, input('A', df));
        await repo.loadAnalysis(a.analysisId);
        const edited = df.clone();
        edited.col('CONC')!.set(15, 123.456);
        const next = await repo.saveAnalysis(a.analysisId, a.source.ref,
          {df: edited, recipe: RECIPE, summaryCache: CACHE, workflow: WF});
        expect(next.ref !== a.source.ref, true, 'a mid-table edit is a data change (full-content hash)');
        const reopened = await h.peer(repo).loadAnalysis(a.analysisId);
        expect(reopened.df.col('CONC')!.get(15), edited.col('CONC')!.get(15), 'the edited cell came back');
        expect(Math.abs(reopened.df.col('CONC')!.get(15) - 123.456) < 1e-3, true, 'it is the edit, not the old value');
        expect(snapshotDiff(edited, reopened.df), null, 'the whole frame came back');
      }), {timeout: TIMEOUT});

    test('PR-13: a derived `~` column is not stored and does not make a Save write a new snapshot', async () =>
      withRepo(async (repo, made) => {
        const s = await repo.createStudy(tag('derived'));
        made.push(s.id);
        const a = await repo.appendAnalysis(s.id, input('A'));
        await repo.loadAnalysis(a.analysisId);
        const computed = frame();
        computed.columns.add(DG.Column.fromInt32Array('~RouteCode',
          Int32Array.from({length: computed.rowCount}, () => 1)));
        const same = await repo.saveAnalysis(a.analysisId, a.source.ref,
          {df: computed, recipe: RECIPE, summaryCache: CACHE, workflow: WF});
        expect(same.ref, a.source.ref, 'the source data did not change → the same snapshot');
        const b = await repo.appendAnalysis(s.id, input('B', computed));
        expect((await repo.readSource(b.source)).col('~RouteCode'), null,
          'a new snapshot holds the source columns only');
      }), {timeout: TIMEOUT});

    test('decision 1: a non-content write elsewhere → Save retries and lands', async () =>
      withRepo(async (repo, made) => {
        const s = await repo.createStudy(tag('retry'));
        made.push(s.id);
        const a = await repo.appendAnalysis(s.id, input('A'));
        await repo.loadAnalysis(a.analysisId);
        await h.foreignRename(repo, a.analysisId);
        const wf: WorkflowState = {...WF, reviewDecisions: DECISION};
        await repo.saveAnalysis(a.analysisId, a.source.ref,
          {df: frame(), recipe: RECIPE, summaryCache: CACHE, workflow: wf});
        const kept = await h.peer(repo).loadAnalysis(a.analysisId);
        expect((kept.entry.reviewDecisions ?? []).length, 1, 'the retried Save landed');
      }), {timeout: TIMEOUT});

    test('F11: data replaced elsewhere → SourceReplacedError', async () =>
      withRepo(async (repo, made) => {
        const s = await repo.createStudy(tag('f11'));
        made.push(s.id);
        const a = await repo.appendAnalysis(s.id, input('A'));
        await repo.loadAnalysis(a.analysisId);
        await h.foreignReplace(repo, a.analysisId);
        const e = await thrown(repo.saveAnalysis(a.analysisId, a.source.ref,
          {df: frame(), recipe: RECIPE, summaryCache: CACHE, workflow: WF}));
        // The repository may see the replacement on its first write (409 → re-read) or not at all
        // (memory) — either way the outcome is the same refusal.
        expect(e instanceof SourceReplacedError, true, `expected SourceReplacedError, got ${e}`);
      }), {timeout: TIMEOUT});

    test('PR-12 / PR-2: an own background write does not fork the next Save; unsaved data is not cached', async () =>
      withRepo(async (repo, made) => {
        const s = await repo.createStudy(tag('pr12'));
        made.push(s.id);
        const df = frame();
        const a = await repo.appendAnalysis(s.id, input('A', df));
        await repo.loadAnalysis(a.analysisId);
        const wf: WorkflowState = {...WF, reviewDecisions: DECISION};
        await repo.writeSummary(a.analysisId, a.source.ref, {summaryCache: CACHE, workflow: wf,
          frameHash: computeSnapshotHash(df)});
        await repo.saveAnalysis(a.analysisId, a.source.ref, {df, recipe: RECIPE, summaryCache: CACHE, workflow: wf});
        // PR-2: a compute on an edited (unsaved) frame must not refresh the cache.
        const edited = frame('SUBJ,TIME,CONC\nA,0,0\nA,1,7\n');
        await repo.writeSummary(a.analysisId, a.source.ref, {summaryCache: {...CACHE, compoundTier: []},
          workflow: {...WF, stageCache: {stage: 'review', pending: 3, notApplicable: 0}},
          frameHash: computeSnapshotHash(edited)});
        const header = (await h.peer(repo).listStudies()).find((x) => x.id === s.id)!;
        expect(header.analyses[0].stageCache?.stage, 'complete', 'a cache from unsaved data was not written');
      }), {timeout: TIMEOUT});

    test('AC-3.5 / SD-14: Move into a study; the container is gone, the snapshot moves with the analysis', async () =>
      withRepo(async (repo, made) => {
        const u = await repo.createStandalone(input(tag('mv')));
        const s = await repo.createStudy(tag('target'));
        made.push(s.id, u.studyId);
        const moved = await repo.moveToStudy(u.studyId, {kind: 'existing', studyId: s.id});
        expect(moved.studyId, s.id);
        const study = await repo.loadStudy(s.id);
        const e = study.analyses.find((x) => x.id === u.analysisId)!;
        expect(e.origin?.kind === 'standalone' && e.origin.containerId === u.studyId, true, 'origin recorded');
        const list = await repo.listStudies();
        expect(list.some((x) => x.id === u.studyId), false, 'the container left Ungrouped');
        // The container's delete must not take the moved analysis's data with it.
        expect(snapshotDiff(frame(), (await repo.loadAnalysis(u.analysisId)).df), null, 'the snapshot still opens');
        const header = list.find((x) => x.id === s.id)!;
        expect(header.analyses[0].movedFromLabel, u.label, 'the analyses/<label> redirect key');
      }), {timeout: TIMEOUT});

    test('delete is final: a study, its analyses and their snapshots are gone; the label is free', async () =>
      withRepo(async (repo, made) => {
        const label = tag('del');
        const s = await repo.createStudy(label);
        made.push(s.id);
        const a = await repo.appendAnalysis(s.id, input('A'));
        const b = await repo.appendAnalysis(s.id, input('B'));
        await repo.deleteStudy(s.id);
        expect((await repo.listStudies()).some((x) => x.id === s.id), false, 'not listed');
        expect(await thrown(repo.loadStudy(s.id)) !== null, true, 'the study does not open');
        expect(await thrown(repo.loadAnalysis(a.analysisId)) !== null, true, 'its analysis does not open');
        for (const x of [a, b])
          expect(await thrown(repo.readSource(x.source)) !== null, true, 'its snapshot is deleted');
        const again = await repo.createStudy(label);
        made.push(again.id);
        expect(again.label, label, 'the label is free again');
      }), {timeout: TIMEOUT});

    test('delete one analysis of a live study: gone with its snapshot; the others stay', async () =>
      withRepo(async (repo, made) => {
        const s = await repo.createStudy(tag('one'));
        made.push(s.id);
        const a = await repo.appendAnalysis(s.id, input('A'));
        const b = await repo.appendAnalysis(s.id, input('B'));
        await repo.deleteAnalysis(a.analysisId);
        const study = await repo.loadStudy(s.id);
        expect(JSON.stringify(study.analyses.map((x) => x.id)), JSON.stringify([b.analysisId]), 'only B remains');
        expect(await thrown(repo.readSource(a.source)) !== null, true, 'A\'s snapshot is deleted');
        expect(snapshotDiff(frame(), await repo.readSource(b.source)), null, 'B\'s snapshot is intact');
        const again = await repo.appendAnalysis(s.id, input('A'));
        expect(again.label, 'A', 'the label is free again');
      }), {timeout: TIMEOUT});

    test('delete a standalone analysis: its container goes with it', async () =>
      withRepo(async (repo, made) => {
        const u = await repo.createStandalone(input(tag('solo')));
        made.push(u.studyId);
        await repo.deleteAnalysis(u.analysisId);
        expect((await repo.listStudies()).some((x) => x.id === u.studyId), false, 'Ungrouped no longer lists it');
        expect(await thrown(repo.readSource(u.source)) !== null, true, 'its snapshot is deleted');
      }), {timeout: TIMEOUT});

    test('S5: a pipeline plan creates, then updates in place; a stale base re-plans', async () =>
      withRepo(async (repo, made) => {
        const owner = `zz-owner-${newId().slice(0, 8)}`;
        // Boolean form: `expect(x, undefined)` defaults the expected value to true.
        expect((await repo.pipelineBase(owner, 'K')).study === undefined, true, 'no study under a fresh key');
        const snap = await repo.writeSnapshot(frame(), new Date().toISOString());
        const id = newId();
        const entry = {...RECIPE, id, label: 'P', source: snap.source, sourceSignature: snap.signature,
          summaryCache: CACHE,
          engineVersion: null, confirmed: null, lzOverrides: [], reviewDecisions: [], stageCache: null, origin: null,
          pipelineFlags: []};
        const planned: NcaStudy = {id: 'new', label: tag('pipe'), kind: 'study', createdAt: '', updatedAt: '',
          access: {canEdit: true, canDelete: true, canShare: true}, analyses: [entry]};
        const sid = await repo.savePipelinePlan({base: undefined, planned, owner, studyKey: 'K'});
        made.push(sid);
        const base = (await repo.pipelineBase(owner, 'K')).study!;
        expect(base.analyses.length, 1);
        const again: NcaStudy = {...base, analyses: [{...base.analyses[0], label: 'P', stageCache: {stage: 'review',
          pending: 1, notApplicable: 0}}]};
        await repo.savePipelinePlan({base, planned: again, owner, studyKey: 'K'});
        const after = (await repo.pipelineBase(owner, 'K')).study!;
        expect(after.analyses.length, 1, 'updated in place, not duplicated');
        expect(after.analyses[0].stageCache?.stage, 'review');
        const e = await thrown(repo.savePipelinePlan({base, planned: again, owner, studyKey: 'K'}));
        expect(e instanceof KeepsChangingError, true, `a plan on a stale base must re-plan, got ${e}`);
      }), {timeout: TIMEOUT});

    // nca-cl-bw-normalised-units AD-8 / F-6 (iii), architect B4: a write that does not recompute the
    // cache (a supersede relabel through the row values, a rename) round-trips BOTH stamps, so
    // neither half of the stale cue is cleared by it.
    test('AD-8: a relabel or rename keeps the engine and the basis stamp; a stale entry stays stale', async () =>
      withRepo(async (repo, made) => {
        const owner = `zz-owner-${newId().slice(0, 8)}`;
        const OLD_ENGINE = '0.0.1-old';
        const mk = async (label: string, engineVersion: string, resultsBasis: string | null) => {
          const snap = await repo.writeSnapshot(frame(), new Date().toISOString());
          return {...RECIPE, id: newId(), label, source: snap.source, sourceSignature: snap.signature,
            summaryCache: CACHE, engineVersion, resultsBasis, confirmed: null, lzOverrides: [], reviewDecisions: [],
            stageCache: null, origin: null, pipelineFlags: []};
        };
        // E: fresh basis, old engine. B: installed engine, pre-basis (null) cache. N: an unknown (null)
        // engine — the LATENT `engine_version ?? INSTALLED` fallback (architect B4) would have filled it.
        const planned: NcaStudy = {id: 'new', label: tag('stamps'), kind: 'study', createdAt: '', updatedAt: '',
          access: {canEdit: true, canDelete: true, canShare: true},
          analyses: [await mk('E', OLD_ENGINE, RESULTS_BASIS), await mk('B', INSTALLED_ENGINE_VERSION, null),
            await mk('N', null as unknown as string, RESULTS_BASIS)]};
        const sid = await repo.savePipelinePlan({base: undefined, planned, owner, studyKey: 'K'});
        made.push(sid);
        const base = (await repo.pipelineBase(owner, 'K')).study!;
        const relabelled: NcaStudy = {...base, analyses: base.analyses.map((a) => ({...a, label: `${a.label} old`}))};
        await repo.savePipelinePlan({base, planned: relabelled, owner, studyKey: 'K'});
        const after = (await repo.pipelineBase(owner, 'K')).study!;
        for (const a of after.analyses) await repo.renameAnalysis(a.id, `${a.label} renamed`);
        const headers = (await repo.listStudies()).find((s) => s.id === sid)!.analyses;
        const byLabel = (l: string) => headers.find((h) => h.label === `${l} old renamed`)!;
        expect(byLabel('E').engineVersion, OLD_ENGINE, 'E: the old engine stamp survives');
        expect(byLabel('E').resultsBasis, RESULTS_BASIS, 'E: the basis stamp survives');
        expect(byLabel('B').engineVersion, INSTALLED_ENGINE_VERSION, 'B: the engine stamp survives');
        expect(byLabel('B').resultsBasis === null, true, 'B: the missing basis is not back-filled');
        expect(byLabel('N').engineVersion === null, true, 'N: the unknown engine is not back-filled (B4)');
        for (const l of ['E', 'B', 'N']) expect(isStale(byLabel(l)), true, `${l}: still stale after relabel + rename`);
        const full = (await repo.loadStudy(sid)).analyses;
        for (const a of full) expect(isEntryStale(a), true, `${a.label}: the full entry is still stale`);
      }), {timeout: TIMEOUT});

    test('BD-3: a frame d42 cannot store exactly is refused before upload', async () => withRepo(async (repo) => {
      // The PR-18 probe's frame — the shape the d42 defect is measured on.
      const df = DG.DataFrame.fromColumns([DG.Column.fromBigInt64Array('Big',
        BigInt64Array.from([BigInt(0), BigInt('9007199254740993'), BigInt(1), BigInt(-2), BigInt(0)]))]);
      df.col('Big')!.set(0, null);
      const e = await thrown(repo.writeSnapshot(df, new Date().toISOString()));
      expect(e instanceof SnapshotFidelityError, true, `expected SnapshotFidelityError, got ${e}`);
    }), {timeout: TIMEOUT});

    // ems-modification decision 2: the time the results changed, not the time they were recomputed.
    test('resultsComputedAt: a recompute of the same numbers keeps it; changed numbers move it', async () =>
      withRepo(async (repo, made) => {
        const s = await repo.createStudy(tag('rca'));
        made.push(s.id);
        const df = frame();
        const a = await repo.appendAnalysis(s.id, input('A', df, tierCache(2)));
        const t0 = (await repo.loadAnalysis(a.analysisId)).entry.resultsComputedAt;
        expect(typeof t0 === 'string' && t0 !== '', true, 'a computed analysis has the time');
        await DG.delay(20);
        await repo.writeSummary(a.analysisId, a.source.ref, {summaryCache: tierCache(2), workflow: WF,
          frameHash: computeSnapshotHash(df)});
        const t1 = (await h.peer(repo).loadAnalysis(a.analysisId)).entry.resultsComputedAt;
        expect(t1, t0, 'the same numbers → the time stays');
        await repo.writeSummary(a.analysisId, a.source.ref, {summaryCache: tierCache(1), workflow: WF,
          frameHash: computeSnapshotHash(df)});
        const t2 = (await h.peer(repo).loadAnalysis(a.analysisId)).entry.resultsComputedAt;
        expect(t2 !== null && t2 !== undefined && t2 > t0!, true, `changed numbers → the time moves (${t0} → ${t2})`);
      }), {timeout: TIMEOUT});

    test('loadStudy + loadCurves: the overview gets every row, and every curve after the load', async () =>
      withRepo(async (repo, made) => {
        const s = await repo.createStudy(tag('curves'));
        made.push(s.id);
        const cache = tierCache(3);
        await repo.appendAnalysis(s.id, input('A', frame(), cache));
        const study = await repo.loadStudy(s.id);
        const ids = pendingCurveIds(study);
        const full = withCurves(study, await repo.loadCurves(ids));
        expect(pendingCurveIds(full).length, 0, 'nothing pending after the load');
        expect(canonical(full.analyses[0].summaryCache), canonical(cache), 'the summary, curves included');
        return `${ids.length} curves loaded lazily`;
      }), {timeout: TIMEOUT});

    test('decision 4: Confirm + Save, then an unconfirmed Setup edit + Save — reopen restores the confirmed recipe',
      async () => withRepo(async (repo, made) => {
        const s = await repo.createStudy(tag('confirm'));
        made.push(s.id);
        const a = await repo.appendAnalysis(s.id, input('A'));
        await repo.loadAnalysis(a.analysisId);
        const confirmed = {recipe: RECIPE, fingerprint: computeRecipeFingerprint(RECIPE),
          at: '2026-01-01T00:00:00.000Z'};
        await repo.saveAnalysis(a.analysisId, a.source.ref, {df: frame(), recipe: RECIPE, summaryCache: CACHE,
          workflow: {...WF, confirmed}});
        const edited = {...RECIPE, units: {time: 'min'}} as unknown as AnalysisRecipe;
        await repo.saveAnalysis(a.analysisId, a.source.ref, {df: frame(), recipe: edited, summaryCache: CACHE,
          workflow: {...WF, confirmed}});
        const back = (await h.peer(repo).loadAnalysis(a.analysisId)).entry;
        expect(canonical(recipeOf(back)), canonical(edited), 'the saved (edited) recipe');
        expect(canonical(back.confirmed?.recipe), canonical(RECIPE),
          'Discard after a reopen restores the confirmed one');
        expect(back.confirmed?.fingerprint, confirmed.fingerprint);
      }), {timeout: TIMEOUT});

    test('AC-3.9: every golden recipe and its workflow state round-trip through the store', async () =>
      withRepo(async (repo, made) => {
        const s = await repo.createStudy(tag('recipes'));
        made.push(s.id);
        const cases = Object.entries((GOLDENS as {cases: Record<string, string>}).cases);
        let n = 0;
        for (const [name, json] of cases) {
          const g = unmasked(JSON.parse(json)) as {recipe: AnalysisRecipe; workflow?: WorkflowState};
          const wf: WorkflowState = g.workflow ?? WF;
          const a = await repo.appendAnalysis(s.id, {label: name, df: frame(), recipe: g.recipe, summaryCache: CACHE,
            workflow: wf, origin: null});
          const back = (await repo.loadAnalysis(a.analysisId)).entry;
          expect(canonical(recipeOf(back)), canonical(g.recipe), `${name}: the recipe`);
          // The stores keep sets, not order or a time format: compare in the normalised form.
          expect(canonical(normLz(back.lzOverrides ?? [])), canonical(normLz(wf.lzOverrides)),
            `${name}: the λz overrides`);
          expect(canonical(normLedgers(back.reviewDecisions ?? [])), canonical(normLedgers(wf.reviewDecisions)),
            `${name}: the decisions`);
          // decision 4: a snapshot whose fingerprint is the recipe's is the recipe (Setup inputs equal).
          const c = wf.confirmed;
          const expected = c === null ? null : {recipe: c.fingerprint === computeRecipeFingerprint(g.recipe) ?
            g.recipe : c.recipe, fingerprint: c.fingerprint, at: isoOf(c.at) ?? ''};
          expect(canonical(back.confirmed ?? null), canonical(expected), `${name}: the confirmed snapshot`);
          expect(canonical(back.stageCache ?? null), canonical(wf.stageCache), `${name}: the stage cache`);
          n++;
        }
        return `${n} recipes round-trip`;
      }), {timeout: 300000});
  });
}

/** A golden with its masked times (`at` / `supersededAt` are `<masked>` in the capture) replaced by
 *  distinct, increasing ISO times — a datetime column stores an instant, not a placeholder, and the
 *  increase keeps each ledger's history in its captured order. */
function unmasked(x: unknown): unknown {
  let n = 0;
  const walk = (v: unknown, key: string): unknown => {
    if (Array.isArray(v)) return v.map((y) => walk(y, ''));
    if (v !== null && typeof v === 'object')
      return Object.fromEntries(Object.entries(v as Record<string, unknown>).map(([k, y]) => [k, walk(y, k)]));
    if (v === '<masked>' && (key === 'at' || key === 'supersededAt'))
      return new Date(Date.UTC(2000, 0, 1) + 1000 * (++n)).toISOString();
    return v;
  };
  return walk(x, '');
}

/** JSON with sorted keys and `undefined` dropped — the store keeps values, not key order. */
function canonical(x: unknown): string {
  const walk = (v: unknown): unknown => {
    if (Array.isArray(v)) return v.map(walk);
    if (v !== null && typeof v === 'object') {
      return Object.fromEntries(Object.keys(v as Record<string, unknown>).sort()
        .filter((k) => (v as Record<string, unknown>)[k] !== undefined)
        .map((k) => [k, walk((v as Record<string, unknown>)[k])]));
    }
    return v;
  };
  return JSON.stringify(walk(x));
}

/** Count the store reads `fn` makes: domain-table reads and snapshot file reads (spec §5). */
async function countReads<T>(fn: () => Promise<T>): Promise<{result: T; reads: number}> {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const tableProto = Object.getPrototypeOf(ncaDb.studies) as any;
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const filesProto = Object.getPrototypeOf(grok.dapi.files) as any;
  const names = ['query', 'get', 'getByKey', 'aggregate', 'count', 'first'];
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const saved: Array<[any, string, any]> = [];
  let reads = 0;
  const wrap = (proto: object, name: string): void => {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    const orig = (proto as any)[name];
    if (typeof orig !== 'function') return;
    saved.push([proto, name, orig]);
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    (proto as any)[name] = function(...args: unknown[]): unknown {reads++; return orig.apply(this, args);};
  };
  for (const n of names) wrap(tableProto, n);
  wrap(filesProto, 'readAsBytes');
  try {
    return {result: await fn(), reads};
  } finally {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    for (const [p, n, f] of saved) (p as any)[n] = f;
  }
}

// eslint-disable-next-line @typescript-eslint/no-explicit-any
type Raw = Record<string, any>;
/** One raw analysis row (the EMS-only tests read what the repository wrote). */
const rawAnalysis = async (id: string): Promise<Raw> => (await ncaDb.analyses.get(id)) as Raw;
const rawRows = async (table: 'analysisResults' | 'lzOverrides' | 'reviewDecisions' | 'pipelineFlags',
  property: string, value: string): Promise<Raw[]> =>
  (await (ncaDb[table] as unknown as {query: (s: unknown) => Promise<unknown[]>}).query(
    {filter: {property, operator: '=', value}, limit: 1000})) as Raw[];

const LZ_ONE = [{profileKey: 'A', timeKey: '2', ordinal: 0, concKey: '3', action: 'force-out' as const,
  user: 'Analyst One', at: '2026-01-02T03:04:05.000Z'}];

// ems-modification step 8: what only the stored rows can show (versions, raw columns, child rows).
category('NCA storage: EMS layout (live)', () => {
  const made: string[] = [];
  const repo = new EmsNcaRepository();
  const cleanup = async (): Promise<void> => {
    for (const id of made.splice(0)) await repo.deleteStudy(id).catch(() => {});
  };
  const tag = (s: string): string => `zz-contract ${s} ${newId().slice(0, 6)}`;

  test('decision 2: a plain reopen writes nothing; a results change moves the version and the time', async () => {
    try {
      const s = await repo.createStudy(tag('reopen'));
      made.push(s.id);
      const df = frame();
      const a = await repo.appendAnalysis(s.id, input('A', df, tierCache(3)));
      await repo.loadAnalysis(a.analysisId);
      const before = await rawAnalysis(a.analysisId);
      await repo.writeSummary(a.analysisId, a.source.ref, {summaryCache: tierCache(3), workflow: WF,
        frameHash: computeSnapshotHash(df)});
      const same = await rawAnalysis(a.analysisId);
      expect(same.version, before.version, 'a reopen that recomputes the same numbers does not bump the version');
      // A recompute with fewer groups: the extra result rows go, the version and the time move.
      await repo.writeSummary(a.analysisId, a.source.ref, {summaryCache: tierCache(1), workflow: WF,
        frameHash: computeSnapshotHash(df)});
      const after = await rawAnalysis(a.analysisId);
      expect(after.version, before.version + 1, 'a results-only change bumps the version once');
      expect(isoOf(after.results_computed_at)! > isoOf(before.results_computed_at)!, true, 'results_computed_at moved');
      const rows = await rawRows('analysisResults', 'analysis_id', a.analysisId);
      expect(rows.length, 1, 'the extra result rows were deleted');
      return `version ${before.version} → ${same.version} (reopen) → ${after.version} (fewer groups)`;
    } finally {await cleanup();}
  }, {timeout: TIMEOUT});

  test('decision 15: result rows carry the stamps; has_cache with zero rows reads back as an empty summary',
    async () => {
      try {
        const s = await repo.createStudy(tag('stamps'));
        made.push(s.id);
        const a = await repo.appendAnalysis(s.id, input('A', frame(), tierCache(2)));
        const b = await repo.appendAnalysis(s.id, input('B'));
        const rows = await rawRows('analysisResults', 'analysis_id', a.analysisId);
        expect(rows.every((r) => r.engine_version === INSTALLED_ENGINE_VERSION && r.results_basis === RESULTS_BASIS),
          true, 'every result row carries the analysis stamps');
        const study = await repo.loadStudy(s.id);
        const eb = study.analyses.find((x) => x.id === b.analysisId)!;
        expect(eb.summaryCache !== null && eb.summaryCache.compoundTier.length === 0, true,
          'computed with zero rows → an empty summary, not "never computed"');
        expect(isEntryStale(eb), false);
        // An engine change makes the entry stale exactly as before (the predicate reads the analysis).
        expect(isEntryStale(study.analyses.find((x) => x.id === a.analysisId)!, '999.0.0'), true);
      } finally {await cleanup();}
    }, {timeout: TIMEOUT});

  test('a move keeps the analysis id, re-parents its result and decision rows, keeps set_* / decided_*', async () => {
    try {
      const wf: WorkflowState = {...WF, lzOverrides: LZ_ONE, reviewDecisions: DECISION};
      const u = await repo.createStandalone({...input(tag('mv'), frame(), tierCache(2)), workflow: wf});
      const s = await repo.createStudy(tag('target'));
      made.push(s.id, u.studyId);
      const lz0 = await rawRows('lzOverrides', 'analysis_id', u.analysisId);
      expect(lz0[0].set_by_label === 'Analyst One' && isoOf(lz0[0].set_at) === LZ_ONE[0].at, true,
        'the override keeps its analyst and time (decision 10, the save-beside guarantee)');
      const moved = await repo.moveToStudy(u.studyId, {kind: 'existing', studyId: s.id});
      expect(moved.analysisId, u.analysisId, 'the same analysis');
      const results = await rawRows('analysisResults', 'analysis_id', u.analysisId);
      const decisions = await rawRows('reviewDecisions', 'analysis_id', u.analysisId);
      expect(results.length === 2 && results.every((r) => r.study_id === s.id), true, 'result rows follow the move');
      expect(decisions.length === 1 && decisions[0].study_id === s.id, true, 'decision rows follow the move');
      expect(decisions[0].decided_by_label === 'u' &&
        isoOf(decisions[0].decided_at) === '2026-01-01T00:00:00.000Z', true,
      'decided_by / decided_at unchanged');
      const lz1 = await rawRows('lzOverrides', 'analysis_id', u.analysisId);
      expect(lz1[0].id === lz0[0].id && lz1[0].set_by_label === 'Analyst One', true,
        'the override row is the same, unchanged');
      const study = await repo.loadStudy(s.id);
      expect(study.analyses[0].summaryCache?.compoundTier.length, 2, 'the study reads the moved results');
    } finally {await cleanup();}
  }, {timeout: TIMEOUT});

  test('decision 2: a move with a fresh summary moves results_computed_at only when the results differ', async () => {
    try {
      const s = await repo.createStudy(tag('mv-target'));
      made.push(s.id);
      const same = await repo.createStandalone(input(tag('mv-same'), frame(), tierCache(2)));
      const other = await repo.createStandalone(input(tag('mv-other'), frame(), tierCache(2)));
      made.push(same.studyId, other.studyId);
      const t0 = isoOf((await rawAnalysis(same.analysisId)).results_computed_at);
      const t1 = isoOf((await rawAnalysis(other.analysisId)).results_computed_at);
      await repo.moveToStudy(same.studyId, {kind: 'existing', studyId: s.id}, {summaryCache: tierCache(2)});
      await repo.moveToStudy(other.studyId, {kind: 'existing', studyId: s.id}, {summaryCache: tierCache(1)});
      expect(isoOf((await rawAnalysis(same.analysisId)).results_computed_at), t0,
        'the same numbers: the compute time stays');
      expect(isoOf((await rawAnalysis(other.analysisId)).results_computed_at)! > t1!, true,
        'different numbers: the compute time moves');
      expect((await rawRows('analysisResults', 'analysis_id', other.analysisId)).length, 1,
        'the moved results are the fresh ones');
    } finally {await cleanup();}
  }, {timeout: TIMEOUT});

  test('save beside (F11): the new analysis keeps each override\'s and decision\'s own analyst and time', async () => {
    try {
      const s = await repo.createStudy(tag('beside'));
      made.push(s.id);
      // updateAnalysisInPlace saves beside through appendAnalysis with the view's workflow state.
      // IR-12: the records carry the analyst's user id beside the name; a save beside keeps both.
      const me = grok.shell.user.id;
      const lzMine = LZ_ONE.map((o) => ({...o, userId: me}));
      const decMine = DECISION.map((l) => ({...l, current: {...l.current, userId: me}}));
      const wf: WorkflowState = {...WF, lzOverrides: lzMine, reviewDecisions: decMine};
      const a = await repo.appendAnalysis(s.id, {...input('A', frame(), tierCache(2)), workflow: wf, origin: null});
      const [lz] = await rawRows('lzOverrides', 'analysis_id', a.analysisId);
      const [d] = await rawRows('reviewDecisions', 'analysis_id', a.analysisId);
      expect(lz.set_by === me && lz.set_by_label === LZ_ONE[0].user && isoOf(lz.set_at) === LZ_ONE[0].at, true,
        'set_by (the user) / set_by_label / set_at are the override\'s own, not the saver\'s or the save time');
      expect(d.decided_by === me && d.decided_by_label === DECISION[0].current!.user &&
        isoOf(d.decided_at) === isoOf(DECISION[0].current!.at), true,
      'decided_by (the user) / decided_by_label / decided_at are the decision\'s own');
      const back = (await repo.loadAnalysis(a.analysisId)).entry;
      expect(back.lzOverrides?.[0].userId === me && back.reviewDecisions?.[0].current?.userId === me, true,
        'the user id reads back on both records');
    } finally {await cleanup();}
  }, {timeout: TIMEOUT});

  test('decision 4 + 7: the stored confirmed_recipe and the recipe search columns', async () => {
    try {
      const s = await repo.createStudy(tag('cols'));
      made.push(s.id);
      const recipe = {...RECIPE, modes: {sparse: {on: true, samplingType: 'batch'}},
        ncaRules: {aucMethod: 'log-linear', blq: {preFirstMeasurable: 'exclude', embedded: 'set-zero',
          afterLast: 'missing', consecutiveAfterLast: 'set-half-lloq'}}} as unknown as AnalysisRecipe;
      const a = await repo.appendAnalysis(s.id, {...input('A'), recipe});
      await repo.loadAnalysis(a.analysisId);
      const confirmed = {recipe, fingerprint: computeRecipeFingerprint(recipe), at: '2026-01-01T00:00:00.000Z'};
      await repo.saveAnalysis(a.analysisId, a.source.ref, {df: frame(), recipe, summaryCache: CACHE,
        workflow: {...WF, confirmed}});
      const r1 = await rawAnalysis(a.analysisId);
      expect(r1.confirmed_recipe === null || r1.confirmed_recipe === undefined, true, 'Confirm + Save → no copy');
      expect(r1.confirmed_fingerprint, confirmed.fingerprint);
      expect([r1.is_sparse, r1.sparse_sampling_type, r1.auc_method, r1.blq_pre_first, r1.blq_embedded,
        r1.blq_after_last, r1.blq_consecutive].join(','),
      'true,batch,log-linear,exclude,set-zero,missing,set-half-lloq',
      'the search columns follow the recipe');
      const edited = {...recipe, ncaRules: {...recipe.ncaRules, aucMethod: 'linear'}} as unknown as AnalysisRecipe;
      await repo.saveAnalysis(a.analysisId, a.source.ref, {df: frame(), recipe: edited, summaryCache: CACHE,
        workflow: {...WF, confirmed}});
      const r2 = await rawAnalysis(a.analysisId);
      expect(canonical(r2.confirmed_recipe), canonical(recipe), 'an unconfirmed edit + Save keeps the snapshot');
      expect(r2.auc_method, 'linear', 'the search column follows the saved recipe');
    } finally {await cleanup();}
  }, {timeout: TIMEOUT});

  test('decision 12: one run over two studies → two pipeline_run rows; a re-run replaces the flag set', async () => {
    try {
      const owner = `zz-owner-${newId().slice(0, 8)}`;
      const runId = newId();
      const rules = '{"lambdaZ":  {"minPoints": 3}, "a": 1}';
      const mapping = {subject: 'ID', time: 'T', concentration: 'C', nested: {z: 1, a: [1, 'x']}};
      const entry = async (flags: string[], run: string) => {
        const snap = await repo.writeSnapshot(frame(), new Date().toISOString());
        const origin: PipelineOrigin = {kind: 'pipeline', runId: run, runAt: '2026-01-01T00:00:00.000Z', user: 'svc',
          engineVersion: INSTALLED_ENGINE_VERSION, mapping, rules, partition: 'P', recipeHash: 'h',
          recipeHashVersion: 1};
        return {...RECIPE, id: newId(), label: 'P', source: snap.source, sourceSignature: snap.signature,
          summaryCache: CACHE, engineVersion: INSTALLED_ENGINE_VERSION, resultsBasis: RESULTS_BASIS, confirmed: null,
          lzOverrides: [], reviewDecisions: [], stageCache: null, origin,
          pipelineFlags: flags.map((p) => ({profileKey: p, verdict: 'Review' as const, reasons: ['r'], info: [],
            pipelineOnly: ['r']}))};
      };
      const plan = async (key: string, flags: string[]): Promise<string> => {
        const planned: NcaStudy = {id: 'new', label: tag(key), kind: 'study', createdAt: '', updatedAt: '',
          access: {canEdit: true, canDelete: true, canShare: true}, analyses: [await entry(flags, runId)]};
        const sid = await repo.savePipelinePlan({base: undefined, planned, owner, studyKey: key});
        made.push(sid);
        return sid;
      };
      await plan('K1', ['S1', 'S2']);
      await plan('K2', ['S3']);
      const runs = await ncaDb.pipelineRuns.query(
        {filter: {property: 'run_id', operator: '=', value: runId}, limit: 10});
      expect(runs.length, 2, 'one run row per (run, study)');
      const base = (await repo.pipelineBase(owner, 'K1')).study!;
      const a = base.analyses[0];
      const opened = (await repo.loadAnalysis(a.id)).entry;
      const o = opened.origin as PipelineOrigin;
      expect(o.rules, rules, 'rules come back byte-identical');
      expect(canonical(o.mapping), canonical(mapping), 'the mapping comes back');
      expect(opened.pipelineFlagCount, 2, 'FIELD-11 count');
      const before = await rawRows('pipelineFlags', 'analysis_id', a.id);
      const rerun = await entry(['S9'], newId());
      await repo.savePipelinePlan({base, planned: {...base, analyses: [{...rerun, id: a.id}]}, owner, studyKey: 'K1'});
      const after = await rawRows('pipelineFlags', 'analysis_id', a.id);
      expect(after.length === 1 && after[0].profile_key === 'S9', true, 'the re-run replaced the flag set');
      expect(before.every((b) => !after.some((x) => x.id === b.id)), true, 'the old flag rows are gone');
      const back = (await repo.loadAnalysis(a.id)).entry;
      expect(back.pipelineFlagCount === 1 && (back.pipelineFlags ?? []).length === 1, true);
      // Step 7: an update of an untouched entry does not touch lz_override / review_decision.
      expect((await rawRows('lzOverrides', 'analysis_id', a.id)).length +
        (await rawRows('reviewDecisions', 'analysis_id', a.id)).length, 0,
      'the untouched entry\'s update wrote no override or decision rows');
    } finally {await cleanup();}
  }, {timeout: TIMEOUT});

  test('AC-8.1: a written decision refuses a change of its core; its supersede update lands', async () => {
    try {
      const s = await repo.createStudy(tag('immut'));
      made.push(s.id);
      const a = await repo.appendAnalysis(s.id, {...input('A'), workflow: {...WF, reviewDecisions: DECISION}});
      const [d] = await rawRows('reviewDecisions', 'analysis_id', a.analysisId);
      let refused = false;
      try {await ncaDb.reviewDecisions.update(d.id, {action: 'exclude'} as never);} catch {refused = true;}
      expect(refused, true, 'the decision core is write-once');
      await repo.loadAnalysis(a.analysisId);
      const reset = [{...DECISION[0], current: null, history: [{...DECISION[0].current,
        supersededBy: 'reset' as const, supersededAt: '2026-01-02T00:00:00.000Z'}]}];
      await repo.saveAnalysis(a.analysisId, a.source.ref, {df: frame(), recipe: RECIPE, summaryCache: CACHE,
        workflow: {...WF, reviewDecisions: reset}});
      const rows = await rawRows('reviewDecisions', 'analysis_id', a.analysisId);
      expect(rows.length === 1 && rows[0].id === d.id && rows[0].superseded_by === 'reset', true,
        'the reset superseded the same row (append-only)');
    } finally {await cleanup();}
  }, {timeout: TIMEOUT});
});

// SD-5: a snapshot is deleted after a failed write only when the server definitively refused it.
// Classified on the typed error, never on message text (a label like "Run 404" is not a refusal).
category('NCA storage: write-failure classification', () => {
  test('a 4xx DomainError refuses; a timeout, a throttle, a transport or a plain error does not', async () => {
    // The server envelope carries the discriminant in `error` (DomainError.code reads body.error).
    const domain = (status: number, code = ''): DG.DomainError => new DG.DomainError('x', status, {error: code});
    expect(isDefinitiveRefusal(domain(403, 'forbidden')), true, '403');
    expect(isDefinitiveRefusal(domain(422, 'validation')), true, '422');
    expect(isDefinitiveRefusal(new DG.DomainVersionConflictError('x', 409, {error: 'version-conflict'})), true, '409');
    expect(isDefinitiveRefusal(domain(408)), false, '408 timeout: the write may have landed');
    expect(isDefinitiveRefusal(domain(429)), false, '429');
    expect(isDefinitiveRefusal(domain(500)), false, '500');
    expect(isDefinitiveRefusal(domain(0)), false, 'transport');
    expect(isDefinitiveRefusal(new Error('Run 404 failed: validation')), false, 'message text alone');
    expect(isConflict(new Error('Run 409')), false, 'a 409 in a label is not a conflict');
    expect(isConflict(new DG.DomainVersionConflictError('x', 409, {error: 'version-conflict'})), true);
    expect(isConflict(domain(409, 'version-conflict')), true, 'a transaction op conflict');
  });
});

category('NCA storage: request counts (live)', () => {
  test('spec §5 + ems-modification: lists 2, a study 3, an open 6, a Save 0 reads, a Move O(1)', async () => {
    const repo = new EmsNcaRepository();
    const made: string[] = [];
    try {
      const s = await repo.createStudy(`zz-contract counts ${newId().slice(0, 6)}`);
      made.push(s.id);
      const a = await repo.appendAnalysis(s.id, input('A'));
      const list = await countReads(() => repo.listStudies());
      expect(list.reads, 2, 'Hub + tree: studies + analyses');
      const study = await countReads(() => repo.loadStudy(s.id));
      expect(study.reads, 3, 'a study: the study, its analyses and its result rows (in parallel)');
      const open = await countReads(() => repo.loadAnalysis(a.analysisId));
      expect(open.reads, 6, 'an open: the row + its 4 child tables (in parallel) + its snapshot');
      const save = await countReads(() => repo.saveAnalysis(a.analysisId, a.source.ref,
        {df: frame(), recipe: RECIPE, summaryCache: CACHE, workflow: WF}));
      expect(save.reads, 0, 'a Save after the open reads nothing');
      const u = await repo.createStandalone(input(`zz-contract mv ${newId().slice(0, 6)}`));
      made.push(u.studyId);
      const move = await countReads(() => repo.moveToStudy(u.studyId, {kind: 'existing', studyId: s.id}));
      expect(move.reads <= 5, true, `a Move is O(1): ${move.reads} reads`);
      return `list ${list.reads}, study ${study.reads}, open ${open.reads}, save ${save.reads}, move ${move.reads}`;
    } finally {
      for (const id of made) await repo.deleteStudy(id).catch(() => {});
    }
  }, {timeout: TIMEOUT});
});

category('NCA storage: permissions (live, two users)', () => {
  test('a stranger sees nothing; a View grant shows the study, its analysis and snapshot; edits stay refused',
    async () => {
      const owner = new EmsNcaRepository();
      const s = await owner.createStudy(`zz-contract perm ${newId().slice(0, 6)}`);
      try {
        const a = await owner.appendAnalysis(s.id, input('A'));
        const outcome = await withRestrictedUser('ncaperm', async (u) => {
          const notes: string[] = [];
          await u.asUser(async () => {
            const repo = new EmsNcaRepository(); // a fresh session cache, as another user has
            expect((await repo.listStudies()).some((x) => x.id === s.id), false, 'a stranger does not list it');
            expect(await thrown(repo.loadAnalysis(a.analysisId)) !== null, true, 'a stranger cannot open it');
            // Decision 9 known gap: the snapshot file is readable by path until the platform binds
            // file access to row permissions. This assertion fails the day it lands.
            const fileReadable = await thrown(grok.dapi.files.readAsBytes(a.source.ref)) === null;
            expect(fileReadable, true, 'file access now follows row permissions — close TODO ' +
              'EXT-DEP-EMS-FILE-ROW-ACCESS and drop the Share-dialog caveat');
            notes.push('stranger reads the snapshot by path (decision 9 known gap)');
          });
          await owner.share(s.id, u.group);
          await u.asUser(async () => {
            const repo = new EmsNcaRepository();
            const seen = (await repo.listStudies()).find((x) => x.id === s.id);
            expect(seen !== undefined, true, 'after the View grant the study is listed');
            expect(seen!.access.canEdit, false, 'View only (decision 10)');
            const opened = await repo.loadAnalysis(a.analysisId);
            expect(snapshotDiff(frame(), opened.df), null, 'the analysis and its snapshot open');
            const e = await thrown(repo.saveAnalysis(a.analysisId, a.source.ref,
              {df: frame(), recipe: RECIPE, summaryCache: CACHE, workflow: WF}));
            expect(e !== null, true, 'a View-only user cannot save');
          });
          // Decision 10 known gap: a row Edit grant made from package code stores the grant but the
          // row stays read-only. This assertion fails the day the platform honours it.
          const [entity] = await grok.dapi.getEntities([s.id]);
          await grok.dapi.permissions.grant(entity, await grok.dapi.groups.find(u.group), true);
          await u.asUser(async () => {
            const seen = (await new EmsNcaRepository().listStudies()).find((x) => x.id === s.id);
            expect(seen?.access.canEdit, false, 'row Edit grants now work — close TODO ' +
              'EXT-DEP-EMS-ROW-EDIT-GRANT and revisit decision 10 (the Share dialog grants View only)');
          });
          notes.push('an Edit grant leaves the row read-only (decision 10 known gap)');
          return notes.join('; ');
        });
        return outcome ?? 'skipped: self-signup is not available on this server';
      } finally {
        await owner.deleteStudy(s.id).catch(() => {});
      }
    }, {timeout: TIMEOUT});
});
