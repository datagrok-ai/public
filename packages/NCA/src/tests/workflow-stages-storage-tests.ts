/**
 * nca-workflow-stages S2.3–S2.5 — the standalone (Ungrouped) analyses on the REAL store (the EMS
 * schema `nca`), live `grok test --host local`.
 *
 * - AC-S2.4.1 (store half): a Setup-stage save (no cache, never confirmed, with overrides and a
 *   superseded decision) round-trips.
 * - AC-S2.5.1: Move to study — the study lists it, Ungrouped no longer does, and the entry's
 *   decisions / overrides / history move identical; `origin` is recorded; the snapshot is not copied.
 * - AC-S1.7.2: a never-confirmed analysis is not moved.
 * (Creator-only visibility is the two-user permissions test of the storage contract suite.)
 * Every study made here is soft-deleted in `finally`.
 */

import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';

import {getNcaRepository, newId, INSTALLED_ENGINE_VERSION} from '../adapter';
import type {AnalysisRecipe, NewAnalysisInput, SummaryCache} from '../adapter';
import {DEFAULT_RULES} from '../adapter/defaults';
import {moveStandaloneIntoStudy} from '../views/nca-add-to-study';

const TIMEOUT = 60000;

/** JSON with sorted keys: EMS json columns keep the values, not the key order. */
function sorted(x: unknown): string {
  return JSON.stringify(x, (_k, v: unknown) => v !== null && typeof v === 'object' && !Array.isArray(v) ?
    Object.fromEntries(Object.entries(v as Record<string, unknown>).sort(([a], [b]) => a < b ? -1 : 1)) : v);
}

const RECIPE: AnalysisRecipe = {
  columns: {
    subject: 'Subject', time: 'Time', concentration: 'conc',
    doseSource: 'value', doseCol: null, doseScalar: 1,
    routeSource: 'value', routeScalar: 'PO', routeCol: null,
    infusionSource: 'none', infusionDurationCol: null, infusionDurationScalar: 0,
    bodyWeightSource: 'none', bodyWeightCol: null, bodyWeightScalar: 0,
    blqFlagCol: null, lloqSource: 'none', lloqCol: null, lloqScalar: 0,
    sexCol: null, sexValueMap: {},
  },
  modes: {nonmem: false, confirmedRoute: null,
    sparse: {on: false, nominalTimeCol: null, animalIdCol: null, samplingType: '(auto)'}},
  ncaRules: DEFAULT_RULES,
  units: {time: 'h', concentration: 'mg/L', dose: 'mg', bodyWeight: 'kg'},
  groupByCols: [], passthroughCols: [], acceptanceRSquared: 0.9,
};

const CACHE: SummaryCache = {
  analysisTier: {nCompounds: 1, nSubjects: 1, route: 'PO', dose: '1 mg',
    status: 'ok', flagSummary: 'No flags'},
  compoundTier: [],
};

function input(label: string, cache: SummaryCache | null): NewAnalysisInput {
  return {
    label, df: DG.DataFrame.fromCsv('Subject,Time,conc\nA,0,1\nA,1,5\n'), recipe: RECIPE, summaryCache: cache,
    origin: null,
    workflow: {
      confirmed: null,
      lzOverrides: [{profileKey: 'A', timeKey: '1', ordinal: 0, concKey: '5', action: 'force-out',
        user: 'probe', at: '2026-09-26T00:00:00.000Z'}],
      reviewDecisions: [{unitKey: 'A', current: null, history: [{action: 'accept',
        reasons: [{code: 'ACC_VISUAL_OK', text: 'ok'}], user: 'probe', at: '2026-09-26T00:00:00.000Z',
        fingerprint: 'f', supersededBy: 'reset', supersededAt: '2026-09-26T00:01:00.000Z'}]}],
      stageCache: {stage: cache === null ? 'setup' : 'complete', pending: 0, notApplicable: 0},
    },
  };
}

category('Workflow stages: standalone store (S2.3–S2.5, live)', () => {
  test('AC-S2.4.1: a Setup-stage standalone save round-trips (no cache, overrides, history)', async () => {
    const repo = getNcaRepository();
    const inp = input(`probe-${newId().slice(0, 8)}`, null);
    const saved = await repo.createStandalone(inp);
    try {
      const {entry, study} = await repo.loadAnalysis(saved.analysisId);
      expect(study.kind, 'standalone');
      expect(entry.summaryCache, null, 'no cache');
      expect(entry.engineVersion ?? null, null, 'no engine stamp without numbers');
      expect(sorted(entry.lzOverrides), sorted(inp.workflow.lzOverrides), 'overrides identical');
      expect(sorted(entry.reviewDecisions), sorted(inp.workflow.reviewDecisions),
        'decisions + superseded history identical');
      const header = (await repo.listStudies()).find((s) => s.id === saved.studyId);
      expect(header !== undefined && header.analyses[0].stageCache?.stage === 'setup', true, 'listed in Ungrouped');
    } finally {
      await repo.deleteStudy(saved.studyId).catch(() => {});
    }
  }, {timeout: TIMEOUT});

  test('AC-S2.5.1: Move to study keeps decisions/overrides/history identical and the snapshot', async () => {
    const repo = getNcaRepository();
    const inp = input(`probe-move-${newId().slice(0, 8)}`, CACHE);
    const saved = await repo.createStandalone(inp);
    let studyId = '';
    try {
      const moved = await moveStandaloneIntoStudy(saved.studyId,
        {kind: 'new', label: `probe-study-${newId().slice(0, 8)}`});
      expect(moved !== null, true, 'move succeeded');
      studyId = moved!.studyId;
      const {entry, study} = await repo.loadAnalysis(saved.analysisId);
      expect(study.id, studyId, 'the analysis now belongs to the study');
      expect(sorted(entry.reviewDecisions), sorted(inp.workflow.reviewDecisions), 'decisions + history');
      expect(sorted(entry.lzOverrides), sorted(inp.workflow.lzOverrides), 'overrides');
      expect(entry.origin?.kind === 'standalone' ? entry.origin.containerId : undefined, saved.studyId,
        'origin recorded for the link redirect');
      expect(entry.source.ref, saved.source.ref, 'the snapshot moves with the row — nothing copied');
      expect(entry.engineVersion, INSTALLED_ENGINE_VERSION);
      const list = await repo.listStudies();
      expect(list.some((s) => s.id === saved.studyId), false, 'Ungrouped no longer lists it');
    } finally {
      if (studyId !== '') await repo.deleteStudy(studyId).catch(() => {});
      await repo.deleteStudy(saved.studyId).catch(() => {});
    }
  }, {timeout: TIMEOUT});

  test('AC-S1.7.2 guard: a never-confirmed analysis is not moved into a study', async () => {
    const repo = getNcaRepository();
    const saved = await repo.createStandalone(input(`probe-guard-${newId().slice(0, 8)}`, null));
    try {
      const moved = await moveStandaloneIntoStudy(saved.studyId, {kind: 'new', label: 'never'});
      expect(moved, null, 'refused');
      expect((await repo.listStudies()).some((s) => s.id === saved.studyId), true, 'the standalone one is untouched');
    } finally {
      await repo.deleteStudy(saved.studyId).catch(() => {});
    }
  }, {timeout: TIMEOUT});
});
