import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {portalDb, ProgramRow} from '../generated/db';
import {loadAllPortfolios, portfolioByCode, readPortfolioFiles} from '../portal';
import {
  cardMismatches, cardParams, cardPortfolios, catalogueCodesOf, mergeCatalogues, ModelCatalogueConfig, outsideValidity,
  unionCodes, validateCatalogue,
} from '../model/catalogue';
import {modelCardOf} from '../catalogue';
import {catalogueCodes} from '../save-run';
import {cardLabel, modelCardElement, modelCardPanelWidget, workspaceModelCardWidget} from '../views/model-card';
import {deleteRunTable, saveModelRun, savedRunsOf} from '../model-runs';
import {_package} from '../package-test';
import {withStandUser} from './stand-users';

const MODEL = 'ProjectPortal:FedBatchTiter';

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

  test('the info panel shows a card under the user\'s portfolios only; admins see all', async () => {
    const mine = new Set(['EX-SMD', 'EX-BIOC']);
    expect(cardPortfolios(['EX-BIOD', 'EX-BIOC', 'EX-SMD'], mine)?.join(), 'EX-BIOC,EX-SMD', 'only the user\'s, in catalogue order');
    expect(cardPortfolios(['EX-BIOD'], mine), null, 'no portfolio of the user\'s lists it');
    expect(cardPortfolios(['EX-BIOD', 'EX-SMD'], null)?.join(), 'EX-BIOD,EX-SMD', 'an admin sees every portfolio');
    expect(cardPortfolios([], mine)?.length, 0, 'a model no portfolio declares belongs to all, unlabelled');
    expect(cardLabel(['EX-SMD']), 'In the catalogue of EX-SMD');
    expect(cardLabel(['EX-SMD', 'EX-BIOC']), 'In the catalogue of EX-SMD, EX-BIOC');
    expect(cardLabel([]), '');
  });

  test('a function is looked up under its nqName and Package:Name, only when catalogued', async () => {
    const known = new Set(['Pkg:Model', 'Other:Model']);
    expect(catalogueCodesOf('Pkg:Model', 'Pkg', 'Model', known).join(), 'Pkg:Model', 'the same code once');
    expect(catalogueCodesOf('Model', 'Pkg', 'Model', known).join(), 'Pkg:Model');
    expect(catalogueCodesOf('Pkg:Other', 'Pkg', 'Other', known).length, 0, 'not catalogued: nothing to read');
    expect(catalogueCodesOf(null, null, 'Model', known).length, 0);
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
      const f = DG.Func.find({package: 'ProjectPortal', name: 'FedBatchTiter'})[0];
      const panel = DG.Func.find({package: 'ProjectPortal', meta: {role: 'workspaceModelPanel'}});
      expect(panel.map((x) => x.name).join(), 'workspaceModelCard', 'found by its role, as PowerPack finds it');
      const w: DG.Widget | null = await panel[0].apply({model: f});
      const root = w?.root;
      expect(root?.querySelector('[data-u2-name="model-card-approval"]')?.textContent, 'draft');
      expect(root?.querySelector('[data-u2-name="model-card-facts"]')?.textContent?.includes('EX-BIOC'), true,
        'the portfolios it belongs to');
      const other = DG.Func.find({package: 'ProjectPortal', name: 'pinnedPrograms'})[0];
      expect(await workspaceModelCardWidget(other), null, 'a function the catalogue does not list');
    });

  test('the Model card info panel shows each catalogue entry under its portfolios, nothing for another function',
    async () => {
      const panel = DG.Func.find({package: 'ProjectPortal', meta: {role: 'panel'}}).find((f) => f.name === 'modelCardPanel');
      expect(panel?.friendlyName, 'Model card', 'an info panel the platform finds by its role');
      expect(panel?.inputs[0].propertyType, DG.TYPE.FUNC);
      const f = DG.Func.find({package: 'ProjectPortal', name: 'FedBatchTiter'})[0];
      const root = (await panel!.apply({model: f}) as DG.Widget | null)?.root;
      expect(root?.querySelectorAll('[data-u2-name="model-card-panel-entry"]').length, 1);
      expect(root?.querySelector('[data-u2-name="model-card-panel-portfolios"]')?.textContent, 'In the catalogue of EX-BIOC');
      expect(root?.querySelector('[data-u2-name="model-card-approval"]')?.textContent, 'draft');

      // a function two portfolios catalogue: one card (cards are per model), labelled with the user's portfolios
      const shared = DG.Func.find({package: 'ProjectPortal', name: 'pinnedPrograms'})[0];
      expect(await modelCardPanelWidget(shared), null, 'a function the catalogue does not list');
      const [m] = await portalDb.models.insert({code: shared.nqName, name: 'Shared test model', portfolios: 'EX-SMD|EX-BIOD'});
      const label = (w: DG.Widget | null): string | null =>
        w?.root.querySelector('[data-u2-name="model-card-panel-portfolios"]')?.textContent ?? null;
      try {
        await catalogueCodes(true);
        const both = await modelCardPanelWidget(shared);
        expect(label(both), 'In the catalogue of EX-SMD, EX-BIOD', 'an admin sees every portfolio');
        expect(!!both?.root.querySelector('[data-u2-name="model-card-empty"]'), true, 'declared without a card');
        await withStandUser('viewer', (u) => u.asUser(async () => {
          expect(label(await modelCardPanelWidget(shared)), 'In the catalogue of EX-SMD', 'the viewer is in EX-SMD, not EX-BIOD');
          expect(label(await modelCardPanelWidget(f)), 'In the catalogue of EX-BIOC');
        }));
        await withStandUser('outsider', (u) => u.asUser(async () => {
          expect(label(await modelCardPanelWidget(shared)), 'In the catalogue of EX-BIOD');
          expect(await modelCardPanelWidget(f), null, 'EX-BIOC\'s model is not the outsider\'s');
        }));
      } finally {
        await portalDb.models.delete(m.id);
        await catalogueCodes(true);
      }
      expect(await modelCardPanelWidget(shared), null, 'gone with its catalogue entry');
    });

  test('selecting a catalogued model adds the Model card pane to its context panel',
    async () => {
      const headers = (): string[] =>
        [...document.querySelectorAll('.grok-prop-panel .d4-accordion-pane-header')].map((e) => e.textContent ?? '');
      grok.shell.o = DG.Func.find({package: 'ProjectPortal', name: 'FedBatchTiter'})[0];
      for (let i = 0; i < 40 && !headers().includes('Model card'); i++)
        await DG.delay(250);
      expect(headers().includes('Model card'), true, headers().join());
      // another function gets none: checked through modelCardPanelWidget above, since the test runner's
      // context panel does not follow a second change of the current object
      expect(document.querySelector('.grok-prop-panel [data-u2-name="model-card-panel-portfolios"]')?.textContent,
        'In the catalogue of EX-BIOC');
    });

  test('a run saved with the example catalogue\'s codes resolves its model, unit operation and run type', async () => {
    const f = DG.Func.find({package: 'ProjectPortal', name: 'FedBatchTiter'})[0];
    const call = await f.prepare({feedRate: 2, temperature: 36, days: 12, seedingDensity: 0.5}).call();
    const r = await saveModelRun(call, 'EX-BIOC', program.code, {unitOperation: 'EX-PROD', runType: 'EX-SIM',
      purpose: 'screening', link: false});
    for (const id of r.tableIds)
      undo.push(() => deleteRunTable(id));
    undo.push(async () => grok.dapi.functions.calls.delete(await grok.dapi.functions.calls.find(r.runId)));
    undo.push(() => portalDb.modelRuns.delete(r.modelRun.id));
    expect(r.notes.some((n) => n.includes('catalogue')), false, r.notes.join('; '));
    const [item] = await savedRunsOf(program.id);
    expect(`${item.model?.code} ${item.stage} ${item.unitOperation} ${item.runType} ${item.run.purpose}`,
      `${MODEL} D5 Production bioreactor Simulated screening`);
  });
});
