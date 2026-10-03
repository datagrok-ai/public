import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {portalDb, ProgramRow} from '../generated/db';
import {loadAllPortfolios, portfolioByCode} from '../portal';
import {COMPARE_RUNS_FUNCTION, CompareCandidate, comparePlan, describePlan} from '../model/compare-runs';
import {criterionRuns, saveModelRun} from '../model-runs';
import {compareRunsTool, criterionComparePlan} from '../compare-runs';
import {openCompareRunsDialog} from '../views/compare-runs-dialog';
import {_package} from '../package-test';

const run = (id: string, model: string, over: Partial<CompareCandidate> = {}): CompareCandidate => ({id, runId: `call-${id}`,
  title: `Run ${id}`, model: `Pkg:${model}`, modelName: model, on: '2026-10-01T00:00:00Z', selected: false, ...over});

category('Portal: run comparison rules', () => {
  test('runs grouped by model; across models the tool matches by name and units', async () => {
    const p = comparePlan([run('1', 'Titer'), run('2', 'Kinetics'), run('3', 'Titer', {runId: null})], true);
    expect(p.groups.map((g) => `${g.modelName}:${g.runs.map((r) => r.id).join('+')}`).join(), 'Titer:1+3,Kinetics:2');
    expect(`${p.crossModel} ${p.notAddable} ${p.unavailable}`, 'true 1 null');
    expect(p.steps[1], 'Pick Titer, tick "Run 1", then Add selected runs', 'a run without a saved call is not to be ticked');
    expect(p.steps.length, 4);
    expect(describePlan(p), '3 runs of 2 models');
  });

  test('not offered without the tool, or with fewer than two runs it can add', async () => {
    expect(comparePlan([run('1', 'A'), run('2', 'A')], false).unavailable, 'Compute2\'s Compare Runs is not installed on this server');
    expect(comparePlan([run('1', 'A'), run('2', 'A', {runId: null})], true).unavailable,
      'Comparing needs two runs with a saved call; one here');
    expect(comparePlan([], true).steps.length, 0);
    expect(describePlan(comparePlan([], true)), 'No candidate runs');
    expect(COMPARE_RUNS_FUNCTION, 'Compute2:CompareRuns');
  });
});

category('Portal: run comparison', () => {
  let program: ProgramRow;
  let criterionId = '';
  const undo: (() => Promise<any>)[] = [];

  before(async () => {
    await loadAllPortfolios(_package);
    const bioc = await portfolioByCode('EX-BIOC');
    const d4 = (await portalDb.stages.getByKey({stage_model_id: bioc.stage_model_id, code: 'D4'}))!;
    const d5 = (await portalDb.stages.getByKey({stage_model_id: bioc.stage_model_id, code: 'D5'}))!;
    criterionId = (await portalDb.gateCriterions.getByKey({stage_id: d5.id, code: 'process_description'}))!.id;
    const [ins] = await portalDb.programs.insert({portfolio_id: bioc.id, code: `CR-${Date.now() % 1e9}`,
      name: 'Run comparison test', stage_model_id: bioc.stage_model_id, stage_id: d4.id, source: 'manual', modality: 'mAb'});
    undo.push(() => portalDb.programs.delete(ins.id));
    program = await portalDb.programs.get(ins.id);
    const f = DG.Func.find({package: 'ProjectPortalCore', name: 'FedBatchTiter'})[0];
    for (const [feedRate, title] of [[2, 'Feed 2'], [4, 'Feed 4']] as [number, string][]) {
      const call = await f.prepare({feedRate, temperature: 36.5, days: 14, seedingDensity: 0.5}).call();
      const r = await saveModelRun(call, 'EX-BIOC', program.code, {unitOperation: 'EX-PROD', runType: 'EX-SIM', title, link: false});
      for (const id of r.tableIds)
        undo.push(async () => grok.dapi.tables.delete(await grok.dapi.tables.find(id)));
      undo.push(async () => grok.dapi.functions.calls.delete(await grok.dapi.functions.calls.find(r.runId)));
      undo.push(() => portalDb.modelRuns.delete(r.modelRun.id));
    }
  });

  after(async () => {
    for (const u of undo.splice(0).reverse())
      await u().catch(() => {});
    document.querySelectorAll('[data-u2-name="compare-runs-dialog"]').forEach((d) => d.remove());
  });

  test('the candidates of process description, by model, as a table', async () => {
    const plan = criterionComparePlan(await criterionRuns(program.id, criterionId));
    expect(`${plan.runs} ${plan.models} ${plan.crossModel}`, '2 1 false');
    expect(plan.groups[0].modelName, 'Fed-batch titer (example)');
    expect(plan.unavailable === null, compareRunsTool() !== null, plan.unavailable ?? '');
    const df: DG.DataFrame = await grok.functions.call('ProjectPortalCore:criterionRunsToCompare',
      {portfolio: 'EX-BIOC', program: program.code, criterion: 'process_description'});
    expect(df.columns.names().join(), 'Model,Run,Run id,Saved on,Selected');
    expect(df.col('Run')!.toList().sort().join(), 'Feed 2,Feed 4');
    expect(df.col('Run id')!.toList().every((x) => !!x), true, 'each run has its saved call');
  });

  test('the dialog lists the runs and opens Compare Runs only where Compute2 is installed', async () => {
    const dlg = openCompareRunsDialog({program: program.code, criterion: {code: 'process_description', name: 'Process description'},
      runs: await criterionRuns(program.id, criterionId)});
    try {
      const root = dlg.root;
      expect(root.querySelector('[data-u2-name="compare-runs-summary"]')!.textContent!.includes('2 runs of 1 model'), true);
      expect(root.querySelector('[data-u2-name="compare-runs-groups"]')!.textContent!.includes('Feed 4'), true);
      const open = root.querySelector('[data-u2-name="compare-runs-open"]') as HTMLButtonElement;
      expect(open.disabled, compareRunsTool() === null);
      expect(!!root.querySelector(compareRunsTool() ? '[data-u2-name="compare-runs-steps"]' :
        '[data-u2-name="compare-runs-unavailable"]'), true);
    } finally {
      dlg.close();
    }
  });
});
