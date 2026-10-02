import * as grok from 'datagrok-api/grok';
import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {
  AutoAdvanceConfig, AutoAdvanceInput, describeAutoAdvance, planAutoAdvance, summarizeAutoAdvance, validateAutoAdvance,
} from '../model/auto-advance';
import {computeReadiness, ReadinessCriterion} from '../model/readiness';
import {CriterionState} from '../model/gates';
import {validatePortfolio} from '../model/config';
import {describeChange, withAdvances, worthNotifying} from '../model/sync-runs';
import {portalDb, ProgramRow} from '../generated/db';
import {loadAllPortfolios, portfolioByCode, readPortfolioFiles, syncPortfolio} from '../portal';
import {autoAdvance} from '../auto-advance';
import {recordedSync} from '../sync-runs';
import {syncHistoryPanel} from '../views/sync-history';
import {_package} from '../package-test';

const stages = [
  {id: 's1', code: 'A', name: 'Alpha', ordinal: 1},
  {id: 's2', code: 'B', name: 'Beta', ordinal: 2},
  {id: 's3', code: 'C', name: 'Gamma', ordinal: 3},
];
const criteria: ReadinessCriterion[] = [
  {id: 'k1', code: 'potent', name: 'Potent', kind: 'data_condition', required: true,
    condition: {where: [{column: 'ic50', op: '<', value: 1}]}},
  {id: 'k2', code: 'panel', name: 'Panel', kind: 'data_condition', required: true, condition: {measured: ['sol']}},
  {id: 'k3', code: 'extra', name: 'Extra', kind: 'data_condition', required: false, condition: {measured: ['herg']}},
];
const good = [{compound: 'X1', ic50: 0.3, sol: 40}];
const bad = [{compound: 'X1', ic50: 3, sol: null}];

const input = (over: Partial<AutoAdvanceInput> & {rows?: any[] | null, recorded?: Record<string, CriterionState>,
  crit?: ReadinessCriterion[], dataError?: string} = {}): AutoAdvanceInput => ({
  program: {code: 'P1', status: 'active', stageId: 's1'},
  stages, config: {through: 'B'},
  readiness: computeReadiness('s2', over.crit ?? criteria, over.recorded ?? {}, over.rows === undefined ? good : over.rows,
    undefined, over.dataError),
  ...over,
});
const why = (i: AutoAdvanceInput) => {
  const o = planAutoAdvance(i);
  return o.to ? `→ ${o.to.code}` : o.reason;
};

category('Portal: auto-advance rules', () => {
  test('a program moves one stage when the data alone meets every required criterion', async () => {
    const o = planAutoAdvance(input());
    expect(`${o.from?.code} → ${o.to?.code}`, 'A → B');
    expect(o.evidence.split('\n')[0], '2 of 2 required criteria of Beta met from data:');
    expect(o.evidence.includes('potent: '), true, o.evidence);
    expect(describeAutoAdvance(o), 'P1: A → B by data');
  });

  test('never past through: the gates above it need an approver', async () => {
    expect(why(input({program: {code: 'P1', status: 'active', stageId: 's2'},
      readiness: computeReadiness('s3', criteria, {}, good)})),
    'Gamma is past Beta, the last stage data may enter: an approver decides');
    expect(why(input({program: {code: 'P1', status: 'active', stageId: 's3'}, readiness: null})),
      'at the last stage of its model');
  });

  test('a waiver counts only when the portfolio says so; a recorded "met" never does', async () => {
    expect(why(input({rows: bad, recorded: {k1: 'waived', k2: 'waived'}})),
      'Beta is not met by data: potent is met only by waiver; panel is met only by waiver');
    expect(why(input({rows: bad, recorded: {k1: 'waived', k2: 'waived'}, config: {through: 'B', countWaivers: true}})), '→ B');
    const manual: ReadinessCriterion[] = [{id: 'm', code: 'assay', name: 'Assay', kind: 'data_condition', required: true, rule: 'Z\' ≥ 0.5'}];
    expect(why(input({crit: manual, recorded: {m: 'met'}})), 'Beta is not met by data: assay is recorded met, not evaluated from data');
    expect(why(input({rows: bad})), 'Beta is not met by data: potent is not met; panel is not met');
  });

  test('held, ended and unstaged programs stay; so does one whose data cannot be read', async () => {
    expect(why(input({program: {code: 'P1', status: 'on_hold', stageId: 's1'}})), 'on hold: only active programs advance by data');
    expect(why(input({program: {code: 'P1', status: 'terminated', stageId: 's1'}})), 'terminated: only active programs advance by data');
    expect(why(input({program: {code: 'P1', status: 'active', stageId: null}})), 'unstaged: a program enters its model by a decision');
    expect(why(input({rows: null, dataError: 'connection refused'})), 'program data could not be read (connection refused)');
    expect(why(input({crit: []})), 'Beta has no required criteria for data to decide');
  });

  test('the setting is checked: portal authority, a stage of the model', async () => {
    const codes = stages.map((s) => s.code);
    expect(validateAutoAdvance(undefined, 'feed', codes).length, 0, 'off by default');
    expect(validateAutoAdvance({through: 'B'}, 'portal', codes).length, 0);
    expect(validateAutoAdvance({through: 'Z'}, 'portal', codes)[0], 'autoAdvance: through "Z" is not a stage of the portfolio\'s model');
    expect(validateAutoAdvance({through: 'B'}, 'feed', codes)[0].startsWith('autoAdvance needs stageAuthority "portal"'), true);
    expect(validateAutoAdvance({through: 'B', countWaivers: 'yes'} as any, 'portal', codes)[0],
      'autoAdvance: countWaivers is true or false');
  });

  test('a run that advanced a program says so, and is worth telling the owners', async () => {
    const changes = withAdvances([{program: 'P2', fields: ['name']}], [{program: 'P1', from: 'A', to: 'B'}]);
    expect(changes.map(describeChange).join(' | '), 'P1: advanced by data A → B | P2: changed name');
    const quiet = {created: 0, updated: 0, flagged: 0, unresolved: 0, unresolvedSpellings: [], overlaps: []};
    expect(worthNotifying({counts: quiet, changes}, []), true);
    expect(worthNotifying({counts: quiet, changes: []}, []), false);
    const config: AutoAdvanceConfig = {through: 'B'};
    expect(summarizeAutoAdvance('X', config, [planAutoAdvance(input()), planAutoAdvance(input({rows: bad}))], true)
      .split('\n').slice(0, 2).join(' / '), 'X: 1 of 2 program(s) would advance by data (through B) /   P1: A → B by data');
    expect(summarizeAutoAdvance('X', null, [], true), 'X: stages advance only by decision (no autoAdvance)');
  });
});

/** Puts EXS-GPC3 (unstaged, its data meets Hit-to-Lead) at Hit Finding for [fn], then puts it back:
 * its stage, and the transitions, criterion statuses and evidence written meanwhile. */
export async function gpc3AtHitFinding(fn: (p: ProgramRow, stage: (code: string) => string) => Promise<void>): Promise<void> {
  const smd = await portfolioByCode('EX-SMD');
  const p = (await portalDb.programs.query({filter: `portfolio_id = "${smd.id}" and code = "EXS-GPC3"`}))[0];
  const stages = await portalDb.stages.query({filter: `stage_model_id = "${smd.stage_model_id}"`});
  const stage = (code: string) => stages.find((s) => s.code === code)!.id;
  const transitionsBefore = new Set((await portalDb.stageTransitions.query({filter: `program_id = "${p.id}"`})).map((t) => t.id));
  const statusesBefore = new Set((await portalDb.criterionStatuses.query({filter: `program_id = "${p.id}"`})).map((t) => t.id));
  const evidenceBefore = new Set((await portalDb.evidences.query({filter: `program_id = "${p.id}"`})).map((t) => t.id));
  await portalDb.programs.update(p.id, {stage_id: stage('HF')});
  try {
    await fn(p, stage);
  } finally {
    await portalDb.programs.update(p.id, {stage_id: p.stage_id ?? null} as any);
    for (const e of await portalDb.evidences.query({filter: `program_id = "${p.id}"`}))
      if (!evidenceBefore.has(e.id)) await portalDb.evidences.delete(e.id).catch(() => undefined);
    for (const t of await portalDb.stageTransitions.query({filter: `program_id = "${p.id}"`}))
      if (!transitionsBefore.has(t.id)) await portalDb.stageTransitions.delete(t.id).catch(() => undefined);
    for (const s of await portalDb.criterionStatuses.query({filter: `program_id = "${p.id}"`}))
      if (!statusesBefore.has(s.id)) await portalDb.criterionStatuses.delete(s.id).catch(() => undefined);
  }
}

category('Portal: auto-advance on the examples', () => {
  before(async () => {
    await loadAllPortfolios(_package);
    await syncPortfolio('EX-SMD');
  });

  test('EX-SMD lets data advance through Hit-to-Lead; the other examples do not opt in', async () => {
    const files = await readPortfolioFiles(_package);
    const smd = files.find((c) => c.code === 'EX-SMD')!;
    expect(JSON.stringify(smd.autoAdvance), '{"through":"H2L","countWaivers":false}');
    expect(validatePortfolio(smd).join('; '), '');
    expect(files.filter((c) => c.code !== 'EX-SMD' && c.code.startsWith('EX-')).every((c) => !c.autoAdvance), true);
    expect(((await portfolioByCode('EX-SMD')).auto_advance as any)?.through, 'H2L', 'loaded');
    expect(await grok.functions.call('ProjectPortalCore:autoAdvancePortfolio', {code: 'EX-BIOC', dryRun: true}),
      'EX-BIOC: stages advance only by decision (no autoAdvance)');
  });

  test('a dry run says what would move and why the rest stay, and writes nothing', async () => {
    await gpc3AtHitFinding(async (p, stage) => {
      const r = await autoAdvance('EX-SMD', true);
      const line = (code: string) => r.outcomes.find((o) => o.program === code);
      expect(describeAutoAdvance(line('EXS-GPC3')!), 'EXS-GPC3: HF → H2L by data');
      expect(line('EXS-KIN1')!.reason, 'unstaged: a program enters its model by a decision');
      expect(line('EXS-ION4')!.reason, 'terminated: only active programs advance by data');
      expect((await portalDb.programs.get(p.id)).stage_id, stage('HF'), 'nothing written');
    });
  });

  test('a real pass advances by data: gate log with source data and the evidence, one stage only', async () => {
    await gpc3AtHitFinding(async (p, stage) => {
      const report: string = await grok.functions.call('ProjectPortalCore:autoAdvancePortfolio', {code: 'EX-SMD', dryRun: false});
      expect(report.split('\n')[0].includes('advanced by data (through H2L)'), true, report);
      expect((await portalDb.programs.get(p.id)).stage_id, stage('H2L'));
      const t = (await portalDb.stageTransitions.query({filter: `program_id = "${p.id}"`}))
        .find((x) => x.to_stage_id === stage('H2L') && x.source === 'data');
      expect(!!t, true, 'the transition is logged with source data');
      expect(`${t!.decision} ${t!.from_stage_id === stage('HF')}`, 'advance true');
      expect(t!.evidence!.includes('validated_series') && t!.evidence!.includes('potency'), true, t!.evidence);
      const again = await autoAdvance('EX-SMD', false);
      expect(again.outcomes.find((o) => o.program === 'EXS-GPC3')!.reason,
        'Lead Optimization is past Hit-to-Lead, the last stage data may enter: an approver decides');
      expect((await portalDb.programs.get(p.id)).stage_id, stage('H2L'), 'one stage, never past through');
    });
  });

  test('a sync advances by data and its history says so', async () => {
    await gpc3AtHitFinding(async (p, stage) => {
      const r = await recordedSync('EX-SMD', 'manual');
      expect(r.error, null);
      expect(r.report.includes('EX-SMD auto-advance: EXS-GPC3 HF → H2L'), true, r.report);
      expect(describeChange(r.changes.find((c) => c.program === 'EXS-GPC3')!).includes('advanced by data HF → H2L'), true);
      expect((await portalDb.programs.get(p.id)).stage_id, stage('H2L'));
      const panel = await syncHistoryPanel(await portfolioByCode('EX-SMD'));
      expect(panel.querySelector('[data-u2-name="sync-run-0"]')!.textContent!.includes('EXS-GPC3: advanced by data HF → H2L'), true);
      expect(panel.querySelector('[data-u2-name="sync-auto-advance"]')!.textContent,
        'Stages advance by data through Hit-to-Lead; waivers do not count; later gates need an approver');
    });
  });

  test('the example is put back: EXS-GPC3 unstaged, no transition by data left', async () => {
    const smd = await portfolioByCode('EX-SMD');
    const p = (await portalDb.programs.query({filter: `portfolio_id = "${smd.id}" and code = "EXS-GPC3"`}))[0];
    expect(p.stage_id ?? null, null);
    const left = (await portalDb.stageTransitions.query({filter: `program_id = "${p.id}"`})).filter((t) => t.source === 'data');
    expect(left.length, 0);
  });
});
