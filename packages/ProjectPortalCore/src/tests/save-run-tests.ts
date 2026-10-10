import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {portalDb, ProgramRow} from '../generated/db';
import {loadAllPortfolios, portfolioByCode} from '../portal';
import {describeSavedRun, programChoices, saveRunDefaults, SaveRunChoices, validateSaveRun} from '../model/save-run';
import {programsOf, saveRunChoices} from '../save-run';
import {openSaveRunDialog} from '../views/save-run-dialog';
import {deleteRunTable, modelRunsOf} from '../model-runs';
import {_package} from '../package-test';

const choices = (over: Partial<SaveRunChoices> = {}): SaveRunChoices => ({
  portfolios: [{code: 'A', name: 'A'}, {code: 'B', name: 'B'}],
  modelPortfolios: ['B'],
  unitOperations: [{code: 'OP1', name: 'One'}, {code: 'OP2', name: 'Two'}],
  modelUnitOperations: ['OP2', 'OP1'],
  runTypes: [{code: 'SIM', name: 'Simulated'}],
  defaultRunType: 'SIM',
  ...over,
});

category('Portal: save to program rules', () => {
  test('the portfolio comes from the page, else the call, else the model\'s catalogue', async () => {
    expect(saveRunDefaults(choices()).portfolio, 'B', 'from the catalogue');
    const fromPage = saveRunDefaults(choices({context: {portfolio: 'A', program: 'P1'}}));
    expect(`${fromPage.portfolio} ${fromPage.program}`, 'A P1');
    const saved = saveRunDefaults(choices({annotations: {portfolio: 'A', program: 'P2', title: 'T',
      unitOperation: 'OP1', runType: null, purpose: 'verification'}}));
    expect(`${saved.portfolio} ${saved.program} ${saved.unitOperation} ${saved.runType} ${saved.purpose} ${saved.title}`,
      'A P2 OP1 SIM verification T');
  });

  test('a code the lists do not hold is never filled in; one portfolio fills itself', async () => {
    const f = saveRunDefaults(choices({modelPortfolios: ['GONE'], modelUnitOperations: ['GONE'], defaultRunType: 'GONE'}));
    expect(`[${f.portfolio}] [${f.unitOperation}] [${f.runType}]`, '[] [] []');
    expect(saveRunDefaults(choices({portfolios: [{code: 'ONLY', name: 'Only'}], modelPortfolios: []})).portfolio, 'ONLY');
    expect(saveRunDefaults(choices()).unitOperation, 'OP2', 'the model\'s first unit operation');
  });

  test('programs: running ones first, then by code; the form needs a portfolio and a visible program', async () => {
    const ps = programChoices([{code: 'B2', name: 'B2', status: 'terminated'}, {code: 'B9', name: 'Nine'}, {code: 'B1', name: 'One'}]);
    expect(ps.map((p) => p.name).join(), 'B1 · One,B9 · Nine,B2');
    const form = {...saveRunDefaults(choices()), program: 'X'};
    expect(validateSaveRun(form, ps)[0], 'X is not a program of B you can see');
    expect(validateSaveRun({...form, program: 'B9'}, ps).length, 0);
    expect(validateSaveRun({...form, portfolio: '', program: ''}, ps).join('; '), 'Choose a portfolio; Choose a program');
  });

  test('the result says where the run went', async () => {
    expect(describeSavedRun('EXB-1', 'D5', {path: 'Root/EXB-1 · X/Drug Substance/Process Lock'}, []),
      'Saved to EXB-1 at D5, linked into Root/EXB-1 · X/Drug Substance/Process Lock');
    expect(describeSavedRun('EXB-1', null, null, []), 'Saved to EXB-1, not linked into a workspace');
    expect(describeSavedRun('EXB-1', 'D5', null, ['already recorded; nothing changed']), 'Already saved to EXB-1; nothing changed');
  });
});

const undo: (() => Promise<unknown>)[] = [];
let program: ProgramRow;

const fedBatch = () => DG.Func.find({package: 'ProjectPortalCore', name: 'FedBatchTiter'})[0];

category('Portal: save to program dialog', () => {
  before(async () => {
    await loadAllPortfolios(_package);
    const bioc = await portfolioByCode('EX-BIOC');
    const d5 = (await portalDb.stages.getByKey({stage_model_id: bioc.stage_model_id, code: 'D5'}))!;
    const [ins] = await portalDb.programs.insert({portfolio_id: bioc.id, code: `SR-${Date.now() % 1e9}`,
      name: 'Save dialog test', stage_model_id: bioc.stage_model_id, stage_id: d5.id, source: 'manual', modality: 'mAb'});
    undo.push(() => portalDb.programs.delete(ins.id));
    program = await portalDb.programs.get(ins.id);
  });

  after(async () => {
    for (const u of undo.splice(0).reverse())
      await u().catch(() => {});
    document.querySelectorAll('[data-u2-name="save-run-dialog"]').forEach((d) => d.remove());
  });

  test('a Fed-batch titer run is offered EX-BIOC, its production bioreactor and simulated run type', async () => {
    const call = await fedBatch().prepare({feedRate: 3, temperature: 36.5, days: 14, seedingDensity: 0.5}).call();
    const f = saveRunDefaults(await saveRunChoices(call));
    expect(`${f.portfolio} ${f.unitOperation} ${f.runType}`, 'EX-BIOC EX-PROD EX-SIM');
    expect((await programsOf('EX-BIOC')).some((p) => p.code === program.code), true, 'the test program is offered');
  });

  test('the dialog saves the run to the chosen program and says so; outside the validity range it warns', async () => {
    const call = await fedBatch().prepare({feedRate: 8, temperature: 36.5, days: 14, seedingDensity: 0.5}).call();
    let line = '';
    const dlg = await openSaveRunDialog({call, context: {portfolio: 'EX-BIOC', program: program.code},
      onSaved: (r, l) => {
        line = l;
        for (const id of r.tableIds)
          undo.push(() => deleteRunTable(id));
        undo.push(async () => grok.dapi.functions.calls.delete(await grok.dapi.functions.calls.find(r.runId)));
        undo.push(() => portalDb.modelRuns.delete(r.modelRun.id));
      }});
    const root = dlg.root;
    const sel = (name: string) => root.querySelector(`[data-u2-name="${name}-input"] select`) as HTMLSelectElement;
    expect(`${sel('portfolio').value} ${sel('program').value} ${sel('unit-operation').value} ${sel('run-type').value}`,
      `EX-BIOC ${program.code} EX-PROD EX-SIM`);
    await DG.delay(500);
    expect(root.querySelector('[data-u2-name="validity-warnings"]')?.textContent?.includes('Feed rate 8 is above'), true,
      root.querySelector('[data-u2-name="validity-warnings"]')?.textContent ?? 'no warning');
    sel('purpose').value = 'screening';
    sel('purpose').dispatchEvent(new Event('change'));
    (root.querySelector('[data-u2-name="save-run-button"]') as HTMLElement).click();
    for (let i = 0; i < 60 && !line; i++)
      await DG.delay(250);
    expect(line, `Saved to ${program.code} at D5, not linked into a workspace`);
    const [run] = await modelRunsOf(program.id);
    expect(`${run?.purpose} ${run?.model_name}`, 'screening ProjectPortalCore:FedBatchTiter');
    expect(document.body.contains(root), false, 'the dialog closed');
  });
});
