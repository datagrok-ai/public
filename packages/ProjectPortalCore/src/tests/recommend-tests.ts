import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {portalDb} from '../generated/db';
import {loadAllPortfolios, readPortfolioFiles} from '../portal';
import {validateCatalogue} from '../model/catalogue';
import {CandidateModel, missingAssays, prefillInputs, recommend, validateApplicability} from '../model/recommend';
import {CriterionReadiness} from '../model/readiness';
import {DataCondition} from '../model/conditions';
import {recommendByCode} from '../recommend';
import {_package} from '../package-test';

const SOL = 'ProjectPortalCore:SolubilityEstimate';
const HERG = 'ProjectPortalCore:HergRisk';

const open = (code: string, condition: DataCondition | null, over: Partial<CriterionReadiness & {rule: string}> = {}) =>
  ({id: code, code, name: code, kind: 'data_condition', required: true, state: 'not_met', basis: 'evaluated', gap: '',
    condition, ...over}) as CriterionReadiness & {condition: DataCondition | null};

const rows = [
  {compound: 'C1', series: 'A', role: 'lead', ic50_um: 0.1, sol: 10, herg: null},
  {compound: 'C2', series: 'A', role: 'lead', ic50_um: null, sol: null, herg: 12},
  {compound: 'C3', series: 'B', role: 'lead', ic50_um: 0.4, sol: null, herg: null},
  {compound: 'C4', series: 'B', role: 'explored', ic50_um: 3, sol: null, herg: null},
];

const model = (code: string, over: Partial<CandidateModel> = {}): CandidateModel =>
  ({code, name: code, portfolios: ['P'], approval: 'draft', ...over});

category('Portal: model recommendation rules', () => {
  test('a card\'s applicability keys are lists; work units are kinds; anything else is refused', async () => {
    expect(validateApplicability({tasks: ['physchem'], stages: ['LO'], workUnits: ['series'], rank: 2}).join('; '), '');
    const errors = validateApplicability({tasks: 'physchem' as any, workUnits: ['plate' as any], colour: ['red'], rank: 'x'} as any);
    expect(errors.some((e) => e.includes('appliesTo.tasks must be a list')), true, errors.join('; '));
    expect(errors.some((e) => e.includes('"plate" is not one of')), true);
    expect(errors.some((e) => e.includes('appliesTo.colour is not a key')), true);
    expect(errors.some((e) => e.includes('rank must be a number')), true);
    expect(validateCatalogue({models: [{code: 'Pkg:M', name: 'M', card: {version: '1', appliesTo: {stages: [''] as any}}}]})
      .some((e) => e.includes('appliesTo.stages')), true, 'the catalogue check reads it');
  });

  test('missing assays: measured columns and numbered where-columns, per criterion, open ones only', async () => {
    const m = missingAssays([
      open('panel', {scope: 'lead_set', measured: ['sol'], minCount: 3}),
      open('potent', {where: [{column: 'ic50_um', op: '<', value: 1}, {column: 'series', op: 'in', value: ['A', 'B']}]}),
      open('done', {measured: ['herg']}, {state: 'met'}),
      open('waived', {measured: ['herg']}, {state: 'waived'}),
      open('words', null, {rule: 'assay validated', required: false}),
    ], rows);
    expect(m.map((x) => `${x.criterion}:${x.column}:${x.items.join('+')}`).join(), 'panel:sol:C2+C3,potent:ic50_um:C2,words:null:',
      'series is a classification, not an assay; C4 is not in the lead set');
    const b = missingAssays([open('panel', {scope: 'lead_set', measured: ['sol']})], rows, {series: 'B'});
    expect(b.map((x) => x.items.join('+')).join(), 'C3', 'only the series asked for');
    expect(missingAssays([open('panel', {measured: ['sol']})], null).length, 0, 'no data, no counts');
  });

  test('recommendation: declared keys must match, absent keys match any, a predicted gap ranks first', async () => {
    const ctx = {portfolio: 'P', modality: 'Small molecule', stages: ['H2L', 'LO'], workUnits: ['campaign', 'series'],
      missingColumns: ['sol']};
    const recs = recommend([
      model('any'),
      model('elsewhere', {portfolios: ['Q']}),
      model('lo', {portfolios: ['Q'], appliesTo: {stages: ['LO'], modalities: ['Small Molecule']}}),
      model('cs', {appliesTo: {stages: ['CS']}}),
      model('bio', {appliesTo: {modalities: ['Antibody']}}),
      model('study', {appliesTo: {workUnits: ['study']}}),
      model('sol', {portfolios: ['Q'], appliesTo: {predicts: ['SOL'], tasks: ['physchem']}}),
      model('old', {status: 'deprecated'}),
      model('retired', {approval: 'retired'}),
      model('approved', {approval: 'approved', appliesTo: {stages: ['lo']}}),
    ], ctx);
    expect(recs.map((r) => r.code).join(), 'sol,lo,approved,any');
    expect(recs[0].reasons.join('; '), 'physchem; predicts SOL, missing for the gate');
    expect(recs[1].reasons.join('; '), 'stage LO; Small molecule');
    expect(recs[3].reasons.join('; '), 'in the P catalogue');
    expect(recommend([model('sol', {appliesTo: {tasks: ['physchem']}}), model('any')], {...ctx, task: 'PhysChem'})
      .map((r) => r.code).join(), 'sol', 'a task asked for needs a model that declares it');
  });

  test('a model\'s dialog opens with the program\'s data and codes', async () => {
    const p = prefillInputs([{name: 'data', type: 'dataframe'}, {name: 'program', type: 'string'},
      {name: 'Series', type: 'string'}, {name: 'threshold', type: 'double'}], {portfolio: 'P', program: 'X1', series: 'A'});
    expect(JSON.stringify(p), '{"values":{"program":"X1","Series":"A"},"needsData":["data"]}');
  });
});

category('Portal: model recommendation', () => {
  before(async () => {
    await loadAllPortfolios(_package);
  });

  test('EX-SMD declares two synthetic predictive models; their cards keep where they apply', async () => {
    const smd = (await readPortfolioFiles(_package)).find((c) => c.code === 'EX-SMD')!;
    expect(validateCatalogue(smd.catalogue!).join('; '), '');
    expect(smd.catalogue!.models!.map((m) => m.code).join(), `${SOL},${HERG}`);
    const m = await portalDb.models.getByKey({code: SOL});
    const card = (await portalDb.modelCards.query({filter: `model_id = "${m!.id}"`}))[0];
    expect(JSON.stringify((card.applies_to as any)?.predicts), '["solubility_um"]');
  });

  test('EXS-KIN1 towards Lead Optimization: missing assays and the models that predict them', async () => {
    const r = await recommendByCode('EX-SMD', 'EXS-KIN1', {stage: 'LO'});
    expect(r.missing.map((m) => `${m.column}:${m.items.join('+')}`).join(),
      'solubility_um:EXC-1006,papp:EXC-1006+EXC-1007,clint:EXC-1006+EXC-1007+EXC-1008',
      'the ADME panel on the lead set; hERG margin is already met');
    expect(r.models.map((x) => x.code).join(), `${SOL},${HERG}`, 'a model predicting a gap first; no bioprocess model');
    expect(r.predictors.get('solubility_um')!.map((x) => x.code).join(), SOL);
    expect(r.models.every((x) => x.installed), true, 'the example scripts ship with the package');
  });

  test('the recommendModels and missingAssays functions; a stage no model is for suggests none', async () => {
    const df: DG.DataFrame = await grok.functions.call('ProjectPortalCore:recommendModels',
      {portfolio: 'EX-SMD', program: 'EXS-KIN1', stage: 'LO'});
    expect(`${df.get('model', 0)} ${df.get('rank', 0)} ${df.get('predicts', 0)}`, `${SOL} 1 solubility_um`);
    const tv: DG.DataFrame = await grok.functions.call('ProjectPortalCore:recommendModels',
      {portfolio: 'EX-SMD', program: 'EXS-KIN1', stage: 'TV'});
    expect(tv.rowCount, 0);
    const adme: DG.DataFrame = await grok.functions.call('ProjectPortalCore:recommendModels',
      {portfolio: 'EX-SMD', program: 'EXS-KIN1', stage: 'LO', task: 'ADMET'});
    expect(adme.rowCount === 1 && adme.get('model', 0) === HERG, true, 'only the model declaring the task');
    const gaps: DG.DataFrame = await grok.functions.call('ProjectPortalCore:missingAssays',
      {portfolio: 'EX-SMD', program: 'EXS-KIN1', stage: 'LO'});
    expect(`${gaps.rowCount} ${gaps.get('assay', 0)} ${gaps.get('missing', 0)} ${gaps.get('predicted by', 0)}`,
      `3 solubility_um 1 ${SOL}`);
    const series: DG.DataFrame = await grok.functions.call('ProjectPortalCore:missingAssays',
      {portfolio: 'EX-SMD', program: 'EXS-KIN1', stage: 'LO', series: 'KIN-B'});
    expect(series.rowCount, 0, 'KIN-B has no lead-set compounds');
  });

  test('a biologics program is offered its own catalogue, not the small-molecule models', async () => {
    const r = await recommendByCode('EX-BIOC', 'EXB-203');
    expect(r.models.map((x) => x.code).join(), 'ProjectPortalCore:FedBatchTiter');
  });

  test('the example models run on the program\'s data', async () => {
    const t = await grok.functions.call('ProjectPortalCore:programDataTable', {portfolio: 'EX-SMD', program: 'EXS-KIN1'});
    const sol: DG.DataFrame = await grok.functions.call(SOL, {data: t});
    expect(sol.rowCount, t.rowCount);
    expect(sol.col('estimated_solubility_um')!.stats.missingValueCount, 0);
    const herg: DG.DataFrame = await grok.functions.call(HERG, {data: t, threshold: 10});
    expect(herg.col('estimated_herg_ic50_um') !== null, true);
  });
});
