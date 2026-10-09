import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {describeSeedPlan, ExampleHistory, planExampleSeed, SeedState, validateExampleHistory} from '../model/examples';
import {PortfolioConfig} from '../model/config';
import {readPortfolioFiles} from '../portal';
import {readExampleHistory, seedExamples} from '../examples';
import {portalDb} from '../generated/db';
import {_package} from '../package-test';

const model = (authority: 'portal' | 'feed', code = 'EX-A'): PortfolioConfig => ({code, name: code, setting: 'discovery',
  stageAuthority: authority, feeds: [], stageModel: {code: `${code}-M`, name: 'M', setting: 'discovery', stages: [
    {code: 'S1', name: 'One', ordinal: 1}, {code: 'S2', name: 'Two', ordinal: 2, criteria: [{code: 'c', name: 'C', kind: 'decision'}]},
    {code: 'S3', name: 'Three', ordinal: 3}]}} as any);

const history = (over: Partial<ExampleHistory> = {}): ExampleHistory => ({
  criteria: [{portfolio: 'EX-A', program: 'P1', stage: 'S2', criterion: 'c', status: 'met', on: '2026-02-01', evidence: 'ok'}],
  decisions: [
    {portfolio: 'EX-A', program: 'P1', decision: 'advance', to: 'S1', on: '2026-01-05', rationale: 'kick-off'},
    {portfolio: 'EX-A', program: 'P1', decision: 'advance', to: 'S2', on: '2026-02-09', rationale: 'go'},
    {portfolio: 'EX-A', program: 'P1', decision: 'hold', on: '2026-03-02', rationale: 'budget'},
  ],
  observations: [{portfolio: 'EX-A', program: 'P1', subject: 'program', on: '2026-02-10', title: 'T', finding: 'F'}],
  follows: [{portfolio: 'EX-A', program: 'P1'}],
  ...over,
});

const state = (over: Partial<SeedState> = {}): SeedState => ({
  programs: new Map([['EX-A P1', {id: 'p1', stage: null, status: 'active'}]]),
  transitions: [], recorded: new Set(), observations: [], follows: new Set(), ...over,
});

category('Portal: example history', () => {
  test('the shipped history fits the portfolio files: staged programs, go, hold and no-go', async () => {
    const h = await readExampleHistory(_package);
    expect(validateExampleHistory(h, await readPortfolioFiles(_package)).join('; '), '');
    expect(h.decisions.every((d) => d.portfolio.startsWith('EX-')), true);
    const kinds = new Set<string>(h.decisions.map((d) => d.outcome ?? d.decision));
    expect(['advance', 'hold', 'terminate', 'go'].every((k) => kinds.has(k)), true, [...kinds].join());
    expect(h.decisions.some((d) => d.program === 'EXS-KIN1' && d.to === 'H2L'), true, 'EXS-KIN1 reaches Hit-to-Lead');
  });

  test('the checks: one stage at a time, records only where the registry owns the stage, EX-* only', async () => {
    const files = [model('portal'), model('feed', 'EX-F')];
    expect(validateExampleHistory(history(), files).length, 0);
    const skip = history({decisions: [{portfolio: 'EX-A', program: 'P1', decision: 'advance', to: 'S2', on: '2026-01-05', rationale: 'x'}]});
    expect(validateExampleHistory(skip, files)[0], 'EX-A P1: advances one stage at a time; next is S1, not S2');
    const record = history({decisions: [{portfolio: 'EX-A', program: 'P1', decision: 'record', outcome: 'go', from: 'S1', to: 'S2',
      on: '2026-01-05', rationale: 'x'}]});
    expect(validateExampleHistory(record, files)[0], 'EX-A P1: a record is for a portfolio whose registry owns the stage');
    const other = history({follows: [{portfolio: 'ACME', program: 'P1'}]});
    expect(validateExampleHistory(other, files).includes('ACME P1: examples go to EX-* portfolios only'), true);
    const criterion = history({criteria: [{portfolio: 'EX-A', program: 'P1', stage: 'S3', criterion: 'c', status: 'met', on: '2026-01-01', evidence: ''}]});
    expect(validateExampleHistory(criterion, files)[0], 'EX-A P1: stage S3 has no criterion c');
  });

  test('the plan: all of it at first, nothing once written, a moved program and a recorded criterion left alone', async () => {
    const first = planExampleSeed(history(), state());
    expect(`${first.decisions.length} ${first.criteria.length} ${first.observations.length} ${first.follows.length}`, '3 1 1 1');
    expect(JSON.stringify(first.moves), '[{"key":"EX-A P1","stage":"S2","status":"on_hold"}]');
    const written = state({
      programs: new Map([['EX-A P1', {id: 'p1', stage: 'S2', status: 'on_hold'}]]),
      transitions: history().decisions.map((d) => ({programId: 'p1', decision: d.decision, rationale: d.rationale})),
      recorded: new Set(['p1 S2 c']), observations: [{programId: 'p1', title: 'T'}], follows: new Set(['p1'])});
    const again = planExampleSeed(history(), written);
    expect(describeSeedPlan(again, true), 'Example history: nothing to write; it is all there');
    // half written, then someone moved the program on: the rest is not forced
    const moved = state({programs: new Map([['EX-A P1', {id: 'p1', stage: 'S3', status: 'active'}]]),
      transitions: [{programId: 'p1', decision: 'advance', rationale: 'kick-off'}]});
    const plan = planExampleSeed(history(), moved);
    expect(plan.decisions.length, 0);
    expect(plan.skipped[0], 'EX-A P1: at S3, not S1 as its example history leaves it; 2 decision(s) left alone');
    expect(planExampleSeed(history(), state({programs: new Map()})).skipped[0], 'EX-A P1: no such program (sync the portfolio first)');
  });

  test('a dry run on the stand writes nothing', async () => {
    const before = await portalDb.stageTransitions.count({});
    const r = await seedExamples(_package, true);
    expect(/^Example history: (would write|nothing to write)/.test(r.report), true, r.report);
    expect(r.plan.decisions.every((d) => d.portfolio.startsWith('EX-')), true);
    expect(await portalDb.stageTransitions.count({}), before);
  });
});
