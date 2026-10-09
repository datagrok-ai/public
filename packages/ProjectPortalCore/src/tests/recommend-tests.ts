import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {portalDb} from '../generated/db';
import {loadAllPortfolios, portfolioByCode, readPortfolioFiles} from '../portal';
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

  test('recommendation: only the portfolio\'s catalogue; declared keys must match; a predicted gap ranks first', async () => {
    const ctx = {portfolio: 'P', modality: 'Small molecule', stages: ['H2L', 'LO'], workUnits: ['campaign', 'series'],
      missingColumns: ['sol']};
    const recs = recommend([
      model('any'),
      model('elsewhere', {portfolios: ['Q']}),
      model('lo', {portfolios: ['Q'], appliesTo: {stages: ['LO'], modalities: ['Small Molecule']}}),
      model('cs', {appliesTo: {stages: ['CS']}}),
      model('bio', {appliesTo: {modalities: ['Antibody']}}),
      model('study', {appliesTo: {workUnits: ['study']}}),
      model('sol', {appliesTo: {predicts: ['SOL'], tasks: ['physchem']}, question: 'How soluble is it?'}),
      model('far', {installed: false, appliesTo: {predicts: ['SOL']}}),
      model('old', {status: 'deprecated'}),
      model('retired', {approval: 'retired'}),
      model('approved', {approval: 'approved', appliesTo: {stages: ['lo']}}),
    ], ctx);
    expect(recs.map((r) => r.code).join(), 'sol,approved,any,far', 'another catalogue\'s model never; not installed last');
    expect(recs[0].reasons.join('; '), 'physchem; predicts SOL, missing for the gate');
    expect(recs[0].question, 'How soluble is it?');
    expect(recs[2].reasons.join('; '), 'in the P catalogue');
    expect(recs[3].reasons.join('; '), 'predicts SOL, missing for the gate; not installed on this server');
    expect(recommend([model('sol', {appliesTo: {tasks: ['physchem']}}), model('any')], {...ctx, task: 'PhysChem'})
      .map((r) => r.code).join(), 'sol', 'a task asked for needs a model that declares it');
  });

  test('recommendation: a unit operation counts when the model and the program both name one', async () => {
    const ctx = {portfolio: 'P', stages: [], workUnits: [], missingColumns: []};
    const ms = [model('harvest', {unitOperations: ['HARVEST']}), model('ferm', {unitOperations: ['FERM', 'SEED']}), model('any')];
    expect(recommend(ms, {...ctx, unitOperations: ['seed']}).map((r) => r.code).join(), 'ferm,any');
    expect(recommend(ms, {...ctx, unitOperations: ['seed']})[0].reasons.join(), 'unit operation seed');
    expect(recommend(ms, ctx).length, 3, 'unknown for the program: the key matches any');
  });

  test('a model\'s dialog opens with the program\'s data and codes', async () => {
    const p = prefillInputs([{name: 'data', type: 'dataframe'}, {name: 'program', type: 'string'},
      {name: 'Series', type: 'string'}, {name: 'threshold', type: 'double'}], {portfolio: 'P', program: 'X1', series: 'A'});
    expect(JSON.stringify(p), '{"values":{"program":"X1","Series":"A"},"needsData":["data"]}');
  });
});

category('Portal: model recommendation', () => {
  // another package's catalogue on the same stand may declare where its models apply (a modality,
  // say) and be offered here too: the tests look at this package's example models only
  const own = (codes: string[]) => codes.filter((c) => c.startsWith('ProjectPortalCore:'));
  const ownRows = (df: DG.DataFrame) => Array.from({length: df.rowCount}, (_, i) => i)
    .filter((i) => String(df.get('model', i)).startsWith('ProjectPortalCore:'));

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
    expect(own(r.models.map((x) => x.code)).join(), `${SOL},${HERG}`, 'a model predicting a gap first; no bioprocess model');
    expect(r.models[0].code, SOL, 'the model predicting a gap comes first');
    expect(r.predictors.get('solubility_um')!.map((x) => x.code).join(), SOL);
    expect(r.models.filter((x) => own([x.code]).length).every((x) => x.installed), true, 'the example scripts ship with the package');
  });

  test('the recommendModels and missingAssays functions; a stage no model is for suggests none', async () => {
    const df: DG.DataFrame = await grok.functions.call('ProjectPortalCore:recommendModels',
      {portfolio: 'EX-SMD', program: 'EXS-KIN1', stage: 'LO'});
    expect(`${df.get('model', 0)} ${df.get('rank', 0)} ${df.get('predicts', 0)}`, `${SOL} 1 solubility_um`);
    // EXS-PRO2 is at most at Target Validation, where no example model applies
    const tv: DG.DataFrame = await grok.functions.call('ProjectPortalCore:recommendModels',
      {portfolio: 'EX-SMD', program: 'EXS-PRO2', stage: 'TV'});
    expect(ownRows(tv).length, 0);
    const adme: DG.DataFrame = await grok.functions.call('ProjectPortalCore:recommendModels',
      {portfolio: 'EX-SMD', program: 'EXS-KIN1', stage: 'LO', task: 'ADMET'});
    expect(ownRows(adme).map((i) => adme.get('model', i)).join(), HERG, 'only the model declaring the task');
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
    expect(r.models.map((x) => x.code).join(), 'ProjectPortalCore:FedBatchTiter', 'nothing from another catalogue');
  });

  test('every model offered is in the portfolio\'s catalogue; a portfolio without one is offered none', async () => {
    const biod = await recommendByCode('EX-BIOD', (await portalDb.programs.query({filter:
      `portfolio_id = "${(await portfolioByCode('EX-BIOD')).id}"`, limit: 1}))[0].code);
    expect(biod.models.map((m) => m.code).join(), '', 'EX-BIOD has no catalogue');
    const kin1 = await recommendByCode('EX-SMD', 'EXS-KIN1', {stage: 'LO'});
    const catalogue = new Set((await portalDb.models.query({limit: 10000}))
      .filter((m) => (m.portfolios ?? '').split('|').includes('EX-SMD')).map((m) => m.code));
    expect(kin1.models.filter((m) => !catalogue.has(m.code)).map((m) => m.code).join(), '');
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
