import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {portalDb, ProgramRow, StageRow} from '../generated/db';
import {loadAllPortfolios, portfolioByCode, programReadiness, readPortfolioFiles} from '../portal';
import {validatePortfolio} from '../model/config';
import {
  checkRunSelection, compareRuns, describeRunFilter, formatRunValue, matchesRunFilter, selectionHistory,
  validateRunCriterion,
} from '../model/run-selection';
import {criterionRuns, runComparison, saveModelRun, SavedModelRun, selectRunForCriterion} from '../model-runs';
import {openChooseRunDialog} from '../views/choose-run-dialog';
import {_package} from '../package-test';

category('Portal: run selection rules', () => {
  test('EX-BIOC\'s process description rests on a production bioreactor run', async () => {
    const bioc = (await readPortfolioFiles(_package)).find((c) => c.code === 'EX-BIOC')!;
    expect(validatePortfolio(bioc).join('; '), '');
    const c = bioc.stageModel.stages.find((s) => s.code === 'D5')!.criteria!.find((x) => x.code === 'process_description')!;
    expect(`${c.kind} ${c.artifact} ${c.runFilter?.unitOperations?.join()}`, 'document model_run EX-PROD');
  });

  test('a criterion takes model runs only, and a filter names lists of codes', async () => {
    expect(validateRunCriterion(undefined, undefined).length, 0);
    expect(validateRunCriterion('model_run', {runTypes: ['SIM']}).length, 0);
    expect(validateRunCriterion('document', undefined)[0], 'artifact "document" is not one of model_run');
    expect(validateRunCriterion(undefined, {models: ['A:B']})[0], 'runFilter needs artifact "model_run"');
    const bad = validateRunCriterion('model_run', {stages: ['D5'], runTypes: 'SIM'} as any);
    expect(bad.some((e) => e.includes('"stages" is not one of')) && bad.some((e) => e.includes('runTypes must be a list')),
      true, bad.join('; '));
  });

  test('the filter narrows by every list it has; a run without the code does not pass', async () => {
    const run = {model: 'P:M', unitOperation: 'PROD', runType: 'SIM'};
    expect(matchesRunFilter(run, null), true);
    expect(matchesRunFilter(run, {unitOperations: ['SEED', 'PROD']}), true);
    expect(matchesRunFilter(run, {unitOperations: ['PROD'], runTypes: ['BENCH']}), false);
    expect(matchesRunFilter({...run, unitOperation: null}, {unitOperations: ['PROD']}), false);
    expect(matchesRunFilter(run, {models: [], runTypes: []}), true, 'an empty list does not narrow');
    expect(describeRunFilter({unitOperations: ['PROD'], runTypes: ['SIM']}, {PROD: 'Production bioreactor'}),
      'Runs on Production bioreactor, of type SIM');
    expect(describeRunFilter(null), 'Every run saved to the program');
  });

  test('contributors and approvers select runs, also where the registry owns the stage; a rationale is needed', async () => {
    const base = {isApprover: false, isContributor: true, programStatus: 'active' as const};
    expect(checkRunSelection(base).length, 0);
    expect(checkRunSelection({...base, isContributor: false})[0], 'Only the portfolio\'s contributors and approvers select runs');
    expect(checkRunSelection({...base, isContributor: false, isApprover: true}).length, 0);
    expect(checkRunSelection({...base, programStatus: 'terminated'})[0], 'The program is terminated');
    expect(checkRunSelection({...base, programStatus: 'on_hold'}).length, 0);
    expect(checkRunSelection({...base, rationale: '  '}).length, 1);
  });

  test('runs side by side: declared order, units in the label, what differs flagged', async () => {
    const rows = compareRuns([
      {title: 'A', inputs: {feedRate: 2, temperature: 36.5}, outputs: {titer: 3.14159, extra: 'x'}},
      {title: 'B', inputs: {feedRate: 4, temperature: 36.5}, outputs: {titer: 4.5}},
    ], {inputs: [{name: 'temperature', caption: 'Temperature', units: '°C'}, {name: 'feedRate', caption: 'Feed rate'}],
      outputs: [{name: 'titer'}]});
    expect(rows.map((r) => `${r.section}:${r.name}:${r.differs ? '≠' : '='}`).join(),
      'input:temperature:=,input:feedRate:≠,output:titer:≠,output:extra:≠');
    expect(rows[0].label, 'Temperature (°C)');
    expect(rows[2].values.join(' | '), '3.142 | 4.5');
    expect(rows[3].values.join(' | '), 'x | —');
    expect(formatRunValue(DG.DataFrame.fromColumns([DG.Column.fromList('int', 'a', [1, 2])])), 'table, 2 rows');
  });

  test('the history: the choice now, and the ones it superseded, newest first', async () => {
    const e = (runTitle: string, role: any, selectedOn: string | null) =>
      ({runTitle, role, selectedOn, selectedBy: null, rationale: null});
    const h = selectionHistory([e('A', 'superseded', '2026-01-01'), e('C', 'selected', '2026-03-01'),
      e('B', 'superseded', '2026-02-01'), e('D', 'candidate', null)]);
    expect(`${h.selected?.runTitle} | ${h.superseded.map((x) => x.runTitle).join()}`, 'C | B,A');
  });
});

/** Everything the platform tests made, removed in reverse in `after`. */
const undo: (() => Promise<unknown>)[] = [];
let program: ProgramRow;
let d5: StageRow;
let criterionId: string;
const saved: SavedModelRun[] = [];

async function saveRun(feedRate: number, unitOperation: string, title: string): Promise<SavedModelRun> {
  const f = DG.Func.find({package: 'ProjectPortalCore', name: 'FedBatchTiter'})[0];
  const call = await f.prepare({feedRate, temperature: 36.5, days: 14, seedingDensity: 0.5}).call();
  const r = await saveModelRun(call, 'EX-BIOC', program.code, {unitOperation, runType: 'EX-SIM', title, link: false});
  for (const id of r.tableIds)
    undo.push(async () => grok.dapi.tables.delete(await grok.dapi.tables.find(id)));
  undo.push(async () => grok.dapi.functions.calls.delete(await grok.dapi.functions.calls.find(r.runId)));
  undo.push(() => portalDb.modelRuns.delete(r.modelRun.id));
  return r;
}

category('Portal: run selection', () => {
  before(async () => {
    await loadAllPortfolios(_package);
    const bioc = await portfolioByCode('EX-BIOC');
    const d4 = (await portalDb.stages.getByKey({stage_model_id: bioc.stage_model_id, code: 'D4'}))!;
    d5 = (await portalDb.stages.getByKey({stage_model_id: bioc.stage_model_id, code: 'D5'}))!;
    criterionId = (await portalDb.gateCriterions.getByKey({stage_id: d5.id, code: 'process_description'}))!.id;
    const [ins] = await portalDb.programs.insert({portfolio_id: bioc.id, code: `RS-${Date.now() % 1e9}`,
      name: 'Run selection test', stage_model_id: bioc.stage_model_id, stage_id: d4.id, source: 'manual', modality: 'mAb'});
    undo.push(() => portalDb.programs.delete(ins.id));
    program = await portalDb.programs.get(ins.id);
    saved.push(await saveRun(2, 'EX-PROD', 'Feed 2'), await saveRun(4, 'EX-PROD', 'Feed 4'),
      await saveRun(3, 'EX-SEED', 'Seed train'));
  });

  after(async () => {
    for (const u of undo.splice(0).reverse())
      await u().catch(() => {});
    saved.splice(0);
    document.querySelectorAll('[data-u2-name="choose-run-dialog"]').forEach((d) => d.remove());
  });

  test('the criterion offers the program\'s production bioreactor runs only', async () => {
    const c = await portalDb.gateCriterions.get(criterionId);
    expect(c.artifact, 'model_run', 'the criterion is loaded with what it rests on');
    const r = await criterionRuns(program.id, criterionId);
    expect(r.candidates.map((i) => i.run.title).sort().join(), 'Feed 2,Feed 4');
    expect(r.excluded, 1, 'the seed-train run is left out');
    expect(r.refused.join('; '), '', 'the test user may select');
    expect(r.selected, null);
  });

  test('two runs side by side: feed rate and titer differ, temperature does not', async () => {
    const r = await criterionRuns(program.id, criterionId);
    const rows = await runComparison(r.candidates);
    const row = (n: string) => rows.find((x) => x.name === n)!;
    expect(`${row('feedRate').differs} ${row('temperature').differs} ${row('titer').differs}`, 'true false true');
    expect(row('feedRate').values.slice().sort().join(), '2,4');
  });

  test('a run selected, then another: the first is superseded and the criterion is met', async () => {
    await selectRunForCriterion(program.id, criterionId, saved[0].modelRun.id, 'lower feed first');
    await selectRunForCriterion(program.id, criterionId, saved[1].modelRun.id, 'higher titer at feed 4');
    const r = await criterionRuns(program.id, criterionId);
    expect(r.selectedRunId, saved[1].modelRun.id);
    expect(`${r.selected?.runTitle} · ${r.selected?.rationale}`, 'Feed 4 · higher titer at feed 4');
    expect(r.superseded.map((e) => e.runTitle).join(), 'Feed 2');
    const ready = await programReadiness('EX-BIOC', program.code, 'D5');
    expect(ready.readiness.criteria.find((c) => c.code === 'process_description')?.state, 'met');
    let refused = '';
    await selectRunForCriterion(program.id, criterionId, saved[0].modelRun.id, ' ').catch((e) => refused = e.message);
    expect(refused.includes('needs a rationale'), true, refused || 'a blank rationale was accepted');
  });

  test('the dialog lists the candidates and compares the ticked ones', async () => {
    const r = await criterionRuns(program.id, criterionId);
    const dlg = await openChooseRunDialog({programId: program.id, runs: r,
      criterion: {id: criterionId, code: 'process_description', name: 'Locked process description',
        runFilter: {unitOperations: ['EX-PROD']}}});
    const root = dlg.root;
    const boxes = Array.from(root.querySelectorAll('[data-u2-name^="candidate-feed-"] input')) as HTMLInputElement[];
    expect(boxes.length, 2);
    expect(root.querySelector('[data-u2-name="run-filter"]')?.textContent,
      'Runs on Production bioreactor · 1 other run of the program do not match');
    for (const b of boxes) {
      b.checked = true;
      b.dispatchEvent(new Event('change'));
    }
    let table: Element | null = null;
    for (let i = 0; i < 40 && !table; i++) {
      await DG.delay(250);
      table = root.querySelector('[data-u2-name="run-comparison"] table');
    }
    expect(!!table, true, 'the comparison is drawn');
    expect(table!.querySelector('[data-u2-name="compare-feedrate"]')?.classList.contains('pp-differs'), true);
    expect(root.querySelector('[data-u2-name="run-selection-history"]')?.textContent?.includes('Feed 2'), true);
    dlg.close();
  });
});
