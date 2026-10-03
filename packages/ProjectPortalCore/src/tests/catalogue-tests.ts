import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {portalDb, ProgramRow} from '../generated/db';
import {loadAllPortfolios, portfolioByCode, readPortfolioFiles} from '../portal';
import {
  cardMismatches, cardParams, mergeCatalogues, ModelCatalogueConfig, outsideValidity, unionCodes, validateCatalogue,
} from '../model/catalogue';
import {modelCardOf} from '../catalogue';
import {modelCardElement, workspaceModelCardWidget} from '../views/model-card';
import {saveModelRun, savedRunsOf} from '../model-runs';
import {_package} from '../package-test';

const MODEL = 'ProjectPortalCore:FedBatchTiter';

const catalogue = (over: Partial<ModelCatalogueConfig> = {}): ModelCatalogueConfig => ({
  unitOperations: [{code: 'OP', name: 'Operation'}],
  runTypes: [{code: 'SIM', name: 'Simulated'}],
  models: [{code: 'Pkg:Model', name: 'Model', unitOperations: ['OP'], defaultRunType: 'SIM', card: {version: '1'}}],
  ...over,
});

category('Portal: model catalogue rules', () => {
  test('the shipped catalogue declares the example model with a card', async () => {
    const bioc = (await readPortfolioFiles(_package)).find((c) => c.code === 'EX-BIOC')!;
    expect(validateCatalogue(bioc.catalogue!).join('; '), '');
    const m = bioc.catalogue!.models!.find((x) => x.code === MODEL)!;
    expect(!!m.card?.question && !!m.card?.intendedUse, true, 'the card states its question and intended use');
    expect(m.unitOperations![0], 'EX-PROD');
    expect(m.defaultRunType, 'EX-SIM');
  });

  test('a model must name its function and the codes the catalogue declares', async () => {
    const errors = validateCatalogue(catalogue({models: [
      {code: 'NoPackage', name: 'A'},
      {code: 'Pkg:B', name: 'B', unitOperations: ['NOPE'], defaultRunType: 'NADA'}]}));
    expect(errors.some((e) => e.includes('Package:Name')), true, errors.join('; '));
    expect(errors.some((e) => e.includes('unit operation NOPE is not declared')), true);
    expect(errors.some((e) => e.includes('run type NADA is not declared')), true);
  });

  test('a card\'s validity range and approval state are checked', async () => {
    const errors = validateCatalogue(catalogue({models: [{code: 'Pkg:M', name: 'M', card: {version: '2',
      approval: 'blessed' as any, supersedes: '2', inputs: [{name: 'x', min: 5, max: 1}, {name: 'x'}]}}]}));
    expect(errors.some((e) => e.includes('approval "blessed"')), true, errors.join('; '));
    expect(errors.some((e) => e.includes('min 5 above max 1')), true);
    expect(errors.some((e) => e.includes('listed twice')), true);
    expect(errors.some((e) => e.includes('supersede itself')), true);
  });

  test('catalogues merge across portfolios: shared codes must agree, a shared model lists both', async () => {
    const a = catalogue();
    const b: ModelCatalogueConfig = {models: [{...a.models![0]}]};
    const merged = mergeCatalogues([{portfolio: 'A', catalogue: a}, {portfolio: 'B', catalogue: b}]);
    expect(merged.errors.join('; '), '', 'B may use the unit operation A declares');
    expect(merged.catalogue.models[0].portfolios.join(), 'A,B');
    const clash = mergeCatalogues([{portfolio: 'A', catalogue: a},
      {portfolio: 'B', catalogue: {unitOperations: [{code: 'OP', name: 'Another name'}]}}]);
    expect(clash.errors.some((e) => e.includes('declared differently by A and B')), true, clash.errors.join('; '));
    expect(unionCodes('A|B', ['B', 'C']), 'A|B|C');
    expect(unionCodes(null, []), '');
  });

  test('card inputs come from the function unless the card lists them, and flag extrapolation', async () => {
    const func = [{name: 'x', caption: 'X', units: 'mg', min: 0, max: 100}, {name: 'y', caption: null}];
    expect(cardParams(null, func).map((p) => `${p.name}:${p.min}-${p.max}`).join(), 'x:0-100,y:undefined-undefined');
    const listed = cardParams([{name: 'x', min: 10, max: 20}], func);
    expect(`${listed.length} ${listed[0].caption} ${listed[0].units} ${listed[0].min}-${listed[0].max}`, '1 X mg 10-20');
    expect(cardMismatches({inputs: [{name: 'z'}]}, {inputs: ['x', 'y'], outputs: []}).join(), 'input z is not an input of the function');
    expect(outsideValidity(listed, {x: 25}).join(), 'X 25 is above the validity range (10–20 mg)');
    expect(outsideValidity(listed, {x: 15}).length, 0);
  });
});

/** Everything a test of this category made, removed in reverse in `after`. */
const undo: (() => Promise<unknown>)[] = [];
let program: ProgramRow;

category('Portal: model catalogue', () => {
  before(async () => {
    await loadAllPortfolios(_package);
    const bioc = await portfolioByCode('EX-BIOC');
    const d5 = (await portalDb.stages.getByKey({stage_model_id: bioc.stage_model_id, code: 'D5'}))!;
    const tag = `MC-${Date.now() % 1e9}`;
    const [ins] = await portalDb.programs.insert({portfolio_id: bioc.id, code: tag, name: 'Catalogue test',
      stage_model_id: bioc.stage_model_id, stage_id: d5.id, source: 'manual', modality: 'mAb'});
    undo.push(() => portalDb.programs.delete(ins.id));
    program = await portalDb.programs.get(ins.id);
  });

  after(async () => {
    for (const u of undo.splice(0).reverse())
      await u().catch(() => {});
  });

  test('loading the portfolios loads the catalogue, and loading again changes nothing', async () => {
    const model = (await portalDb.models.getByKey({code: MODEL}))!;
    expect(!!model, true, `${MODEL} is catalogued`);
    expect(model.portfolios?.split('|').includes('EX-BIOC'), true, model.portfolios);
    expect(model.current_version, '1.0');
    for (const code of ['EX-SEED', 'EX-PROD', 'EX-HARV', 'EX-PROA', 'EX-UFDF'])
      expect(!!await portalDb.unitOperations.getByKey({code}), true, `unit operation ${code}`);
    const cards = await portalDb.modelCards.query({filter: `model_id = "${model.id}"`});
    expect(cards.length, 1);
    await loadAllPortfolios(_package);
    const again = (await portalDb.models.getByKey({code: MODEL}))!;
    const cardsAgain = await portalDb.modelCards.query({filter: `model_id = "${model.id}"`});
    expect(`${again.version} ${cardsAgain.length} ${cardsAgain[0].version}`, `${model.version} 1 ${cards[0].version}`,
      'the second load wrote nothing');
  });

  test('the card joins the catalogue with the function: validity ranges from the card, outputs from the function', async () => {
    const v = (await modelCardOf(MODEL))!;
    expect(v.installed, true);
    expect(v.card?.approval, 'draft');
    const feed = v.inputs.find((p) => p.name === 'feedRate')!;
    expect(`${feed.caption} ${feed.units} ${feed.min}-${feed.max}`, 'Feed rate % v/v per day 0-6');
    const titer = v.outputs.find((p) => p.name === 'titer')!;
    expect(`${titer.caption} ${titer.units}`, 'Final titer g/L');
    expect(v.unitOperations.map((u) => u.code).join(), 'EX-PROD,EX-SEED');
    expect(v.defaultRunType?.code, 'EX-SIM');
    expect(await modelCardOf('NoSuch:Model'), null);
  });

  test('the card renders its question, inputs and approval', async () => {
    const el = modelCardElement((await modelCardOf(MODEL))!);
    const q = el.querySelector('[data-u2-name="model-card-question"]');
    expect(q?.textContent?.startsWith('What final titer'), true, q?.textContent ?? 'no question');
    expect(el.querySelectorAll('[data-u2-name="model-card-inputs"] tbody tr').length, 4);
    expect(el.querySelector('[data-u2-name="model-card-input-feedrate"]')?.textContent?.includes('0 – 6'), true);
    expect(el.querySelector('[data-u2-name="model-card-approval"]')?.textContent, 'draft');
  });

  test('the Workspace on Home shows a catalogued model\'s card beside its preview, nothing for another function',
    async () => {
      const f = DG.Func.find({package: 'ProjectPortalCore', name: 'FedBatchTiter'})[0];
      const panel = DG.Func.find({package: 'ProjectPortalCore', meta: {role: 'workspaceModelPanel'}});
      expect(panel.map((x) => x.name).join(), 'workspaceModelCard', 'found by its role, as PowerPack finds it');
      const w: DG.Widget | null = await panel[0].apply({model: f});
      const root = w?.root;
      expect(root?.querySelector('[data-u2-name="model-card-approval"]')?.textContent, 'draft');
      expect(root?.querySelector('[data-u2-name="model-card-facts"]')?.textContent?.includes('EX-BIOC'), true,
        'the portfolios it belongs to');
      const other = DG.Func.find({package: 'ProjectPortalCore', name: 'pinnedPrograms'})[0];
      expect(await workspaceModelCardWidget(other), null, 'a function the catalogue does not list');
    });

  test('a run saved with the example catalogue\'s codes resolves its model, unit operation and run type', async () => {
    const f = DG.Func.find({package: 'ProjectPortalCore', name: 'FedBatchTiter'})[0];
    const call = await f.prepare({feedRate: 2, temperature: 36, days: 12, seedingDensity: 0.5}).call();
    const r = await saveModelRun(call, 'EX-BIOC', program.code, {unitOperation: 'EX-PROD', runType: 'EX-SIM',
      purpose: 'screening', link: false});
    for (const id of r.tableIds)
      undo.push(async () => grok.dapi.tables.delete(await grok.dapi.tables.find(id)));
    undo.push(async () => grok.dapi.functions.calls.delete(await grok.dapi.functions.calls.find(r.runId)));
    undo.push(() => portalDb.modelRuns.delete(r.modelRun.id));
    expect(r.notes.some((n) => n.includes('catalogue')), false, r.notes.join('; '));
    const [item] = await savedRunsOf(program.id);
    expect(`${item.model?.code} ${item.stage} ${item.unitOperation} ${item.runType} ${item.run.purpose}`,
      `${MODEL} D5 Production bioreactor Simulated screening`);
  });
});
