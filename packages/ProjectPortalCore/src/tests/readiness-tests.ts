import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {DataCondition, evaluateCondition, holds, validateCondition} from '../model/conditions';
import {computeReadiness, ReadinessCriterion} from '../model/readiness';
import {PortfolioConfig, validatePortfolio} from '../model/config';
import {portalDb} from '../generated/db';
import {loadAllPortfolios, portfolioByCode, programContext, programReadiness, recordDecision, setCriterion,
  syncPortfolio, verdictFor} from '../portal';
import {_package} from '../package-test';

const rows = [
  {compound: 'C1', series: 'A', role: 'lead', ic50_um: 0.2, sol: 40, papp: 10, clint: 8},
  {compound: 'C2', series: 'A', role: 'lead', ic50_um: 0.6, sol: 30, papp: null, clint: 9},
  {compound: 'C3', series: 'A', role: 'backup', ic50_um: 2.0, sol: 50, papp: 12, clint: 7},
  {compound: 'C4', series: 'B', role: 'explored', ic50_um: 5.0, sol: null, papp: null, clint: null},
  {compound: 'C5', series: null, role: 'explored', ic50_um: '0.9', sol: null, papp: null, clint: null},
];

category('Portal: data conditions', () => {
  test('operators: numbers, text, lists, presence; blanks never compare', async () => {
    expect(holds({x: 0.5}, {column: 'x', op: '<', value: 1}), true);
    expect(holds({x: '0.5'}, {column: 'x', op: '<', value: 1}), true, 'a numeric string compares as a number');
    expect(holds({x: null}, {column: 'x', op: '<', value: 1}), false);
    expect(holds({x: null}, {column: 'x', op: '!=', value: 'a'}), false);
    expect(holds({X: 'Validated'}, {column: 'x', op: 'in', value: ['validated', 'lead']}), true,
      'columns and text match case-insensitively');
    expect(holds({x: 'hit'}, {column: 'x', op: 'not_in', value: ['validated']}), true);
    expect(holds({x: true}, {column: 'x', op: '=', value: true}), true);
    expect(holds({x: ''}, {column: 'x', op: 'absent'}), true);
    expect(holds({x: 0}, {column: 'x', op: 'present'}), true, 'zero is a value');
  });

  test('program scope: count rows that pass where', async () => {
    const r = evaluateCondition({where: [{column: 'ic50_um', op: '<', value: 1}], minCount: 3}, rows);
    expect(r.qualifying, 3);
    expect(r.met, true, r.summary);
    expect(evaluateCondition({where: [{column: 'ic50_um', op: '<', value: 0.1}]}, rows).met, false);
  });

  test('lead-set scope and a measured panel: coverage and the rows that miss it', async () => {
    const c: DataCondition = {scope: 'lead_set', measured: ['sol', 'papp', 'clint'], minCount: 3};
    const r = evaluateCondition(c, rows);
    expect(r.inScope, 3);
    expect(r.qualifying, 2);
    expect(r.met, false);
    expect(r.coverage.papp, 2);
    expect(r.gaps.map((g) => `${g.id}:${g.missing.join('+')}`).join(), 'C2:papp');
    expect(evaluateCondition(c, rows, ['lead']).inScope, 2, 'the portfolio decides what a lead role is');
  });

  test('series scope: enough series must each reach the count', async () => {
    const c: DataCondition = {scope: 'series', where: [{column: 'ic50_um', op: '<', value: 3}], minCount: 2};
    const r = evaluateCondition(c, rows);
    expect(r.seriesMet!.join(), 'A');
    expect(r.met, true);
    expect(evaluateCondition({...c, minSeries: 2}, rows).met, false);
    expect(r.inScope, 4, 'a row without a series is outside a series scope');
  });

  test('a condition that cannot be evaluated is refused at load', async () => {
    expect(validateCondition({where: [{column: 'x', op: '<', value: 'high'}]}).length, 1);
    expect(validateCondition({where: [{column: 'x', op: 'in', value: 'a'}]}).length, 1);
    expect(validateCondition({where: [{column: 'x', op: 'like' as any, value: 'a'}]}).length, 1);
    expect(validateCondition({minSeries: 2, measured: ['a']}).length, 1, 'minSeries without a series scope');
    expect(validateCondition({}).length, 1, 'a condition any row satisfies');
    expect(validateCondition({scope: 'lead_set', measured: ['sol'], minCount: 10}).length, 0);
  });

  test('a portfolio with conditions needs a program dataset that selects by program', async () => {
    const p: PortfolioConfig = {code: 'P', name: 'P', setting: 'discovery', stageAuthority: 'portal', feeds: [],
      stageModel: {code: 'M', name: 'M', setting: 'discovery', stages: [{code: 'A', name: 'A', ordinal: 1,
        criteria: [{code: 'c', name: 'c', kind: 'data_condition', condition: {measured: ['x']}}]}]}};
    expect(validatePortfolio(p).some((e) => e.includes('no programData')), true);
    p.programData = {connection: 'System:Datagrok', sql: 'select 1'};
    expect(validatePortfolio(p).some((e) => e.includes('@program')), true);
    p.programData.sql = 'select * from t where program = @program';
    expect(validatePortfolio(p).join('; '), '');
    p.stageModel.stages[0].criteria![0].kind = 'document';
    expect(validatePortfolio(p).some((e) => e.includes('only data conditions take a condition')), true);
  });
});

const criteria: ReadinessCriterion[] = [
  {id: 'k1', code: 'potent', name: 'Potent', kind: 'data_condition', required: true,
    condition: {where: [{column: 'ic50_um', op: '<', value: 1}]}},
  {id: 'k2', code: 'panel', name: 'Panel', kind: 'data_condition', required: true,
    condition: {scope: 'lead_set', measured: ['sol', 'papp', 'clint'], minCount: 3}},
  {id: 'k3', code: 'doc', name: 'Profile', kind: 'document', required: true},
  {id: 'k4', code: 'manual', name: 'Assay', kind: 'data_condition', required: false, rule: 'Z\' ≥ 0.5'},
];

category('Portal: readiness', () => {
  test('computed from data, recorded for documents, and a waiver wins', async () => {
    const r = computeReadiness('s', criteria, {k3: 'in_progress'}, rows);
    expect(r.criteria.map((c) => `${c.code}=${c.state}/${c.basis}`).join(),
      'potent=met/evaluated,panel=not_met/evaluated,doc=in_progress/recorded,manual=missing/recorded');
    expect(`${r.requiredMet}/${r.requiredTotal}`, '1/3');
    expect(r.gaps.length, 2, 'the optional criterion is not a gap');
    expect(r.gaps[0].includes('C2 lacks papp'), true, r.gaps[0]);
    const waived = computeReadiness('s', criteria, {k2: 'waived', k3: 'met'}, rows);
    expect(waived.score, 1);
    expect(waived.states.k1, 'met', 'evaluated states feed the gate rules');
  });

  test('a recorded "met" does not outvote the data', async () => {
    const r = computeReadiness('s', criteria, {k2: 'met'}, rows);
    expect(r.criteria.find((c) => c.code === 'panel')!.state, 'not_met');
  });

  test('without the program data, evaluated criteria fall back to what is recorded', async () => {
    const r = computeReadiness('s', criteria, {k1: 'met'}, null, undefined, 'connection refused');
    expect(r.criteria.find((c) => c.code === 'potent')!.basis, 'recorded');
    expect(r.criteria.find((c) => c.code === 'panel')!.gap.includes('connection refused'), true);
  });

  test('a stage without required criteria is ready', async () => {
    expect(computeReadiness('s', [], {}, []).score, 1);
  });
});

let ready: Promise<void> | undefined;
const setUp = () => ready ??= (async () => {
  await loadAllPortfolios(_package);
  for (const code of ['EX-SMD', 'EX-BIOD', 'EX-BIOC'])
    await syncPortfolio(code);
})();

const line = async (portfolio: string, program: string, stage: string) => {
  const r = await programReadiness(portfolio, program, stage);
  return r.readiness.criteria.map((c) => `${c.code}=${c.state}`).join();
};

category('Portal: readiness on the example portfolios', () => {
  before(async () => { await setUp(); });

  test('a validated, potent series passes the hit-to-lead data conditions', async () => {
    expect(await line('EX-SMD', 'EXS-KIN1', 'H2L'), 'validated_series=met,potency=met');
    expect(await line('EX-SMD', 'EXS-GPC3', 'H2L'), 'validated_series=met,potency=met');
  });

  test('hits only: neither condition holds, and the gap list says why', async () => {
    const r = (await programReadiness('EX-SMD', 'EXS-PRO2', 'H2L')).readiness;
    expect(r.criteria.map((c) => `${c.code}=${c.state}`).join(), 'validated_series=not_met,potency=not_met');
    expect(r.gaps.length, 2);
    expect(r.score, 0);
  });

  test('lead optimization: the ADME panel counts complete lead-set compounds and names the rest', async () => {
    const r = (await programReadiness('EX-SMD', 'EXS-KIN1', 'LO')).readiness;
    const panel = r.criteria.find((c) => c.code === 'adme_panel')!;
    expect(`${panel.result!.qualifying}/${panel.result!.inScope}`, '5/8');
    expect(panel.result!.gaps.map((g) => g.id).join(), 'EXC-1006,EXC-1007,EXC-1008');
    expect(r.criteria.find((c) => c.code === 'herg_margin')!.state, 'met', 'optional, and met');
    expect(`${r.requiredMet}/${r.requiredTotal}`, '0/1');
  });

  test('a program with no data is told so', async () => {
    const r = (await programReadiness('EX-SMD', 'EXS-ION4', 'H2L')).readiness;
    expect(r.gaps.every((g) => g.includes('no data')), true, r.gaps.join('; '));
  });

  test('the same language reads development lots', async () => {
    expect(await line('EX-BIOC', 'EXB-201', 'D3'), 'tox_lot=met,ind_quality=missing');
    expect(await line('EX-BIOC', 'EXB-104', 'D3'), 'tox_lot=not_met,ind_quality=missing',
      'a tox lot that is not released does not count');
  });

  test('the gate uses the evaluated conditions; a waiver opens it', async () => {
    const sm = await portfolioByCode('EX-SMD');
    const kin = (await portalDb.programs.query({filter: `portfolio_id = "${sm.id}" and code = "EXS-KIN1"`}))[0];
    // a manual twin of EXS-KIN1 shares its data through the program code: take it TV → HF → H2L
    const code = `TEST-RDY-${Date.now() % 1e9}`;
    const [ins] = await portalDb.programs.insert({portfolio_id: sm.id, code, name: code,
      stage_model_id: sm.stage_model_id, source: 'manual'});
    try {
      await recordDecision(ins.id, 'advance', null, 'kick-off');
      let ctx = await programContext(ins.id);
      const hf = ctx.stages.find((s) => s.code === 'HF')!;
      await setCriterion(ins.id, ctx.criteria.find((c) => c.stageId === hf.id)!.id, 'met');
      await recordDecision(ins.id, 'advance', null, 'assay validated');
      ctx = await programContext(ins.id);
      expect(ctx.readiness!.criteria.every((c) => c.state === 'not_met'), true, 'a program with no data');
      expect(verdictFor(ctx, 'advance').allowed, false);
      // the same program under KIN1's code would pass: readiness reads the data by code
      expect((await programReadiness('EX-SMD', kin.code, 'H2L')).readiness.score, 1);
      const h2l = ctx.stages.find((s) => s.code === 'H2L')!;
      for (const c of ctx.criteria.filter((c) => c.stageId === h2l.id))
        await setCriterion(ins.id, c.id, 'waived', 'test');
      ctx = await programContext(ins.id);
      expect(verdictFor(ctx, 'advance').allowed, true, verdictFor(ctx, 'advance').reasons.join('; '));
    } finally {
      await portalDb.programs.delete(ins.id);
    }
  });
});
