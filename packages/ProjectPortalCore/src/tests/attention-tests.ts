import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {loadAllPortfolios, portfolioByCode, readPortfolioFiles} from '../portal';
import {validatePortfolio} from '../model/config';
import {AttentionRule, attentionMessage, evaluateAttention, validateAttentionRules} from '../model/attention';
import {attentionRulesOf, openAttentionItem, programAttentionOf} from '../attention';
import {COL, portfolioTables} from '../portfolio-table';
import {_package} from '../package-test';

const rows = [
  {compound: 'C1', role: 'lead', ic50_um: 0.05, sol: 40, herg: 22},
  {compound: 'C2', role: 'lead', ic50_um: 0.2, sol: null, herg: null},
  {compound: 'C3', role: 'backup', ic50_um: 0.4, sol: 50, herg: 4},
  {compound: 'C4', role: 'explored', ic50_um: 5.0, sol: null, herg: null},
];

const potent: AttentionRule = {code: 'potent_no_herg', name: 'Potent, no hERG',
  condition: {where: [{column: 'ic50_um', op: '<', value: 0.5}], measured: ['herg']}};
const noSol: AttentionRule = {code: 'lead_no_sol', name: 'Lead without solubility', severity: 'info',
  condition: {scope: 'lead_set', measured: ['sol']}};
const liability: AttentionRule = {code: 'herg_low', name: 'hERG under 10', severity: 'critical',
  condition: {where: [{column: 'herg', op: '<', value: 10}]}};

category('Portal: attention rules', () => {
  test('a row trips a rule when in scope, where holds, and a measured column is missing', async () => {
    const a = evaluateAttention([potent, noSol, liability], rows);
    expect(a.rules.map((r) => `${r.code}:${r.items.join('+')}`).join(), 'potent_no_herg:C2,lead_no_sol:C2,herg_low:C3');
    expect(a.count, 2, 'C2 trips two rules and counts once');
    expect(a.worst, 'critical');
    expect(a.items.map((i) => `${i.item}:${i.severity}:${i.rules.join('+')}`).join(),
      'C3:critical:herg_low,C2:warning:potent_no_herg+lead_no_sol', 'most severe first; a rule\'s default severity is warning');
    expect(a.summary, '2 items need attention: Potent, no hERG 1, Lead without solubility 1, hERG under 10 1');
  });

  test('the lead set follows the portfolio\'s lead roles; no data, nothing flagged', async () => {
    expect(evaluateAttention([noSol], rows, ['lead', 'explored']).rules[0].items.join(), 'C2,C4');
    const none = evaluateAttention([potent], []);
    expect(`${none.count} ${none.worst} ${none.summary}`, '0 null Nothing needs attention');
  });

  test('a message fills in the row\'s values and what is missing; the default names the rule', async () => {
    const r = {...potent, message: '{compound}: IC50 {ic50_um} uM, missing {missing}, sol {sol}'};
    expect(attentionMessage(r, rows[1], ['herg']), 'C2: IC50 0.2 uM, missing herg, sol —');
    expect(attentionMessage(potent, rows[1], ['herg']), 'Potent, no hERG; ic50_um 0.2; no herg');
    expect(evaluateAttention([potent], rows).flags[0].values.ic50_um, 0.2);
  });

  test('rules are checked: code, name, severity, a where or measured, no counts', async () => {
    expect(validateAttentionRules([potent, noSol, liability]).length, 0);
    const errors = validateAttentionRules([potent, {...potent, name: ''}, {code: 'x', name: 'X', severity: 'high' as any,
      condition: {minCount: 2}}, {code: 'y', name: 'Y', condition: {where: [{column: 'a', op: '<', value: 'b'}]}},
    {code: 'z', name: 'Z'} as any]);
    for (const e of ['potent_no_herg: code is used twice', 'name is required', 'severity "high"', 'x: minCount and minSeries',
      'x: give it where or measured', 'y: predicate on a: "<" needs a number', 'z: condition is required'])
      expect(errors.some((x) => x.includes(e)), true, `${e} in: ${errors.join('; ')}`);
    const smd = (await readPortfolioFiles(_package)).find((c) => c.code === 'EX-SMD')!;
    expect(validatePortfolio({...smd, programData: undefined, stageModel: {...smd.stageModel,
      stages: smd.stageModel.stages.map((s) => ({...s, criteria: []}))}}).some((e) => e.includes('attention rules are declared but')), true);
  });

  test('the examples: three rules on EX-SMD, one on EX-BIOC lots, all valid', async () => {
    const all = await readPortfolioFiles(_package);
    const smd = all.find((c) => c.code === 'EX-SMD')!;
    const bioc = all.find((c) => c.code === 'EX-BIOC')!;
    expect((smd.attention ?? []).map((r) => r.code).join(), 'herg_unmeasured,lead_no_solubility,herg_liability');
    expect((bioc.attention ?? []).map((r) => r.code).join(), 'tox_lot_unreleased');
    expect(validatePortfolio(smd).join('; '), '');
    expect(validatePortfolio(bioc).join('; '), '');
  });
});

const opened: DG.ViewBase[] = [];

category('Portal: attention', () => {
  before(async () => {
    await loadAllPortfolios(_package);
  });

  after(async () => {
    for (const v of opened.splice(0))
      v.close();
  });

  test('EXS-KIN1: four potent compounds lack hERG, one lead lacks solubility', async () => {
    const r = await programAttentionOf('EX-SMD', 'EXS-KIN1');
    expect(attentionRulesOf(r.portfolio).length, 3, 'the rules are loaded with the portfolio');
    expect(r.attention.rules.map((x) => `${x.code}:${x.items.join('+')}`).join(),
      'herg_unmeasured:EXC-1002+EXC-1004+EXC-1006+EXC-1007,lead_no_solubility:EXC-1006,herg_liability:');
    expect(r.attention.count, 4);
    expect(r.attention.items.find((i) => i.item === 'EXC-1006')!.messages.join(' | '),
      'EXC-1006: IC50 0.3 uM, no hERG IC50 | EXC-1006 (lead): no solubility');
  });

  test('the programAttention function: one row per item and rule', async () => {
    const df: DG.DataFrame = await grok.functions.call('ProjectPortalCore:programAttention', {portfolio: 'EX-SMD', program: 'EXS-GPC3'});
    expect(df.rowCount, 1);
    expect(`${df.get('item', 0)} ${df.get('rule', 0)} ${df.get('severity', 0)}`, 'EXC-3002 herg_liability critical');
    const bioc = await programAttentionOf('EX-BIOC', 'EXB-104');
    expect(bioc.attention.items.map((i) => i.item).join(), 'EXL-104-02', 'the unreleased tox lot');
  });

  test('a portfolio without rules is refused', async () => {
    let message = '';
    try {
      await programAttentionOf('EX-BIOD', 'EXB-101');
    } catch (e: any) {
      message = e?.message ?? String(e);
    }
    expect(message, 'EX-BIOD declares no attention rules');
  });

  test('the portfolio table counts attention per program; no rules, no count', async () => {
    const t = await portfolioTables('EX-SMD');
    const df = t.programs;
    const count = (code: string) => {
      const i = Array.from({length: df.rowCount}, (_, k) => k).find((k) => df.get(COL.program, k) === code)!;
      return `${df.get(COL.attention, i)}|${df.get(COL.attentionRules, i) ?? ''}`;
    };
    expect(count('EXS-KIN1'), '4|Potent, hERG not measured 4; Lead without solubility 1');
    expect(count('EXS-GPC3'), '1|hERG IC50 under 10 uM 1');
    expect(count('EXS-PRO2'), '0|');
    const biod = (await portfolioTables('EX-BIOD')).programs;
    expect(biod.col(COL.attention)!.stats.missingValueCount, biod.rowCount);
  });

  test('a flagged compound opens the program data at its row, the flagged rows selected', async () => {
    const tv = await openAttentionItem('EX-SMD', 'EXS-KIN1', 'EXC-1006', ['EXC-1002', 'EXC-1004', 'EXC-1006', 'EXC-1007']);
    opened.push(tv);
    const df = tv.dataFrame;
    expect(df.get('compound', df.currentRowIdx), 'EXC-1006');
    expect(df.selection.trueCount, 4);
  });
});
