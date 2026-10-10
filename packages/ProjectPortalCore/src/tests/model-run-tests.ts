import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {portalDb, ProgramRow, StageRow} from '../generated/db';
import {loadAllPortfolios, portfolioByCode} from '../portal';
import {buildPortfolioWorkspace} from '../workspace';
import {deleteRunTable, modelRunsOf, saveModelRun, SavedModelRun, selectRunForCriterion} from '../model-runs';
import {readRunAnnotations, RUN_OPTIONS, runFolderPath, runOptions} from '../model/model-runs';
import {_package} from '../package-test';

const MODEL = 'ProjectPortalCore:FedBatchTiter';

category('Portal: model run rules', () => {
  test('annotations round-trip through call options, and keep the tags already there', async () => {
    const o = runOptions({portfolio: 'EX-BIOC', program: 'EXB-1', stage: 'D5', unitOperation: 'PB',
      runType: 'SIM', purpose: 'characterization', title: 'T'}, ['mine', 'EXB-1']);
    expect(o.tags.join(), 'mine,EXB-1,portal,D5');
    const back = readRunAnnotations(o)!;
    expect([back.portfolio, back.program, back.stage, back.unitOperation, back.runType, back.purpose, back.title].join(),
      'EX-BIOC,EXB-1,D5,PB,SIM,characterization,T');
    expect(back.description, null);
    expect(readRunAnnotations({title: 'not a portal run'}), null);
  });

  test('a run goes to the folder asked for, else the current stage\'s first folder, else the program folder', async () => {
    expect(runFolderPath('Drug Substance/Process Lock|Regulatory', null).join('>'), 'Drug Substance>Process Lock');
    expect(runFolderPath('Drug Substance/Process Lock', ' Analytical / Methods ').join('>'), 'Analytical>Methods');
    expect(runFolderPath(null, null).length, 0);
    expect(runFolderPath('', '/').length, 0);
  });

  test('the example model gives a titer and a daily time course', async () => {
    const f = DG.Func.find({package: 'ProjectPortalCore', name: 'FedBatchTiter'})[0];
    expect(!!f, true, `${MODEL} is not published`);
    expect(f.options['role'], 'model', 'annotated as a model (meta.role)');
    const fed = await f.prepare({feedRate: 3, temperature: 36.5, days: 14, seedingDensity: 0.5}).call();
    const batch = await f.prepare({feedRate: 0, temperature: 36.5, days: 14, seedingDensity: 0.5}).call();
    expect(fed.outputs['titer'] > 2 && fed.outputs['titer'] < 6, true, `titer ${fed.outputs['titer']}`);
    expect(fed.outputs['titer'] > batch.outputs['titer'], true, 'feeding raises the titer');
    const tc = fed.outputs['timeCourse'] as DG.DataFrame;
    expect(tc.rowCount, 15);
    expect(tc.columns.names().join(), 'Day,VCD,Viability,Titer');
  });
});

/** Everything a test of this category made, removed in reverse in `after`. */
const undo: (() => Promise<unknown>)[] = [];
let program: ProgramRow;
let processLock: StageRow;
let tag: string;
let first: SavedModelRun;

async function runModel(feedRate: number): Promise<DG.FuncCall> {
  const f = DG.Func.find({package: 'ProjectPortalCore', name: 'FedBatchTiter'})[0];
  return await f.prepare({feedRate, temperature: 36.5, days: 14, seedingDensity: 0.5}).call();
}

/** Undoes a saved run: the model_run row, the link, the call and its uploaded tables. */
function forget(r: SavedModelRun): void {
  for (const id of r.tableIds)
    undo.push(() => deleteRunTable(id));
  undo.push(async () => grok.dapi.functions.calls.delete(await grok.dapi.functions.calls.find(r.runId)));
  if (r.folder)
    undo.push(() => grok.dapi.spaces.id(r.folder!.id).removeEntity(r.runId));
  undo.push(() => portalDb.modelRuns.delete(r.modelRun.id));
}

const subspaces = async (id: string) =>
  await grok.dapi.spaces.id(id).children.ofTypes('Project', false).list({pageSize: 1000}) as DG.Project[];

async function removeFolder(id: string): Promise<void> {
  for (const c of await subspaces(id))
    await removeFolder(c.id);
  const space = await grok.dapi.spaces.find(id).catch(() => null);
  if (space)
    await grok.dapi.spaces.delete(space);
}

category('Portal: model runs', () => {
  before(async () => {
    await loadAllPortfolios(_package);
    const cmc = await portfolioByCode('EX-BIOC');
    processLock = (await portalDb.stages.getByKey({stage_model_id: cmc.stage_model_id, code: 'D5'}))!;
    tag = `MR-${Date.now() % 1e9}`;
    const [ins] = await portalDb.programs.insert({portfolio_id: cmc.id, code: tag, name: 'Model run test',
      stage_model_id: cmc.stage_model_id, stage_id: processLock.id, source: 'manual', modality: 'mAb'});
    undo.push(() => portalDb.programs.delete(ins.id));
    await buildPortfolioWorkspace('EX-BIOC', false, ins.id);
    program = await portalDb.programs.get(ins.id);
    undo.push(async () => program.space_id && removeFolder(program.space_id));
    const [op] = await portalDb.unitOperations.insert({code: `${tag}-PB`, name: 'Production bioreactor'});
    undo.push(() => portalDb.unitOperations.delete(op.id));
    const [rt] = await portalDb.runTypes.insert({code: `${tag}-SIM`, name: 'Simulated'});
    undo.push(() => portalDb.runTypes.delete(rt.id));
    if (!await portalDb.models.getByKey({code: MODEL})) {
      const [m] = await portalDb.models.insert({code: MODEL, name: 'Fed-batch titer (example)',
        package: 'ProjectPortalCore', category: 'Bioprocess', status: 'experimental'});
      undo.push(() => portalDb.models.delete(m.id));
    }
    first = await saveModelRun(await runModel(3), 'EX-BIOC', tag, {unitOperation: `${tag}-PB`,
      runType: `${tag}-SIM`, purpose: 'characterization', description: 'baseline feed'});
    forget(first);
  });

  after(async () => {
    for (const u of undo.splice(0).reverse())
      await u().catch(() => {});
  });

  test('a run saved at process lock is linked into the stage\'s folder', async () => {
    expect(program.space_id != null, true, 'the test program has a workspace folder');
    expect(first.folder?.path.endsWith('/Drug Substance/Process Lock'), true, first.folder?.path);
    expect(first.folder!.id !== program.space_id, true);
    // listing a space's linked children fails once it links a FuncCall (platform); the project's links do not
    const space = await grok.dapi.projects.find(first.folder!.id);
    expect(space.links.some((e) => e.id === first.runId), true, 'the run is linked into the folder');
  });

  test('the annotations survive a reload of the call from the server', async () => {
    const back = await grok.dapi.functions.calls.include('inputs, outputs').find(first.runId);
    const a = readRunAnnotations(back.options as any)!;
    expect(a != null, true, 'the reloaded call carries the portal annotations');
    expect([a.portfolio, a.program, a.stage, a.unitOperation, a.runType, a.purpose].join(),
      `EX-BIOC,${tag},D5,${tag}-PB,${tag}-SIM,characterization`);
    expect(a.title, `Fed-batch titer (example) · ${tag}`);
    expect(a.description, 'baseline feed');
    expect((back.options['tags'] as string[]).includes(tag), true);
    expect(back.inputs['feedRate'], 3);
    expect(back.outputs['titer'], first.call.outputs['titer']);
    const found = await grok.dapi.functions.calls.filter(`options.${RUN_OPTIONS.program} = "${tag}"`).list();
    expect(found.map((c) => c.id).join(), first.runId, 'the program\'s runs are found by option');
  });

  test('the run is recorded as a model_run with the stage at run; by run id while the soft reference is refused', async () => {
    const runs = await modelRunsOf(program.id);
    expect(runs.length, 1);
    const r = runs[0];
    expect(r.run_id, first.runId);
    expect(r.stage_id, processLock.id);
    expect(r.model_name, MODEL);
    expect(r.purpose, 'characterization');
    expect(r.model_id != null, true, 'the catalogued model is linked');
    expect(r.unit_operation_id != null && r.run_type_id != null, true);
    if (first.recordedBy === 'func_call')
      expect(r.func_call, first.runId, 'the soft reference is written');
    else {
      expect(r.func_call ?? null, null);
      expect(!!first.softRefError, true, 'the refusal is reported');
    }
    const again = await saveModelRun(first.call, 'EX-BIOC', tag, {link: false});
    expect(again.modelRun.id, r.id, 'saving the same call again records no second run');
  });

  test('a gate document criterion selects the run, and a later run supersedes it', async () => {
    const criterion = (await portalDb.gateCriterions.getByKey({stage_id: processLock.id, code: 'process_description'}))!;
    expect(criterion.kind, 'document');
    const pick = await selectRunForCriterion(program.id, criterion.id, first.modelRun.id, 'baseline feed meets target');
    expect(pick.role, 'selected');
    const status = (await portalDb.criterionStatuses.query(
      {filter: `program_id = "${program.id}" and gate_criterion_id = "${criterion.id}"`}))[0];
    expect(status.status, 'met');
    expect(status.evidence?.includes(first.runId), true, status.evidence);

    const second = await saveModelRun(await runModel(5), 'EX-BIOC', tag, {purpose: 'verification'});
    forget(second);
    await selectRunForCriterion(program.id, criterion.id, second.modelRun.id, 'higher feed, confirmed');
    const artifacts = await portalDb.criterionArtifacts.query({filter: `criterion_status_id = "${status.id}"`});
    const roleOf = (runId: string) => artifacts.find((x) => x.model_run_id === runId)?.role;
    expect(`${roleOf(first.modelRun.id)},${roleOf(second.modelRun.id)}`, 'superseded,selected');
    expect((await modelRunsOf(program.id)).length, 2);
  });
});
