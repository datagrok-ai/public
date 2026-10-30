import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {awaitCheck, category, expect, test} from '@datagrok-libraries/test/src/test';
import {historyUtils} from '@datagrok-libraries/compute-utils/history-utils/src/history-utils';
import {
  getRunAuthorId, getRunCreated, getRunName, CREATED_BY_OPTION, CREATED_AT_OPTION, NESTED_ENTITY_OPTION,
} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {getModelFilter} from '@datagrok-libraries/compute-utils/model-catalog/src/model-handler';
import {deleteSavedRun, deleteSavedWorkflow} from './utils';
import {
  RunsChange, addRunsToSpace, addToSpace, bySavedDesc, canEditSpace, defaultRunName, deleteOrRemoveRuns,
  listMyModelRuns,
  loadRunAuthors,
  modelLabel, parseTags, placeRunFilter, planSpaceAdd, runVersion, runsChanged$,
} from '../components/RunStorage/run-storage';
import {RunDetailsForm} from '../components/RunStorage/RunDetailsForm';
import {listSpace} from '../components/RunStorage/space-index';
import {loadInstanceState} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/funccall-utils';
import {isFuncCallSerializedState, PipelineSerializedState} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {
  runWorkflow, saveWorkflowRun,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/headless';

export const DAY_CALCULATION = 'Bioreactorworkflowdemo:DayCalculation';
const BIOREACTOR_WORKFLOW = 'Bioreactorworkflowdemo:BioreactorWorkflow';

let counter = 0;
export const uniqueName = (prefix: string) => `${prefix} ${Date.now()}-${counter++}`;
const ids = (runs: DG.FuncCall[]) => runs.map((run) => run.id);

/** The runs of a space as the run browser lists them; with models, only theirs, newest first. */
export async function spaceRuns(space: DG.Project, models?: DG.Func | DG.Func[]): Promise<DG.FuncCall[]> {
  const {runs} = await listSpace(space).toPromise();
  return models ? runs.filter(placeRunFilter(models)).sort(bySavedDesc) : runs;
}
const SOMEONE_ELSE = '00000000-0000-0000-0000-000000000001';
const LONG_AGO = '2020-01-02T03:04:05.000Z';

const DAY_INPUTS = {
  solutionAdded: 0.5, substrateAdded: 0, incomingVolume: 1, incomingBiomass: 0.1, incomingSubstrate: 10,
  maximumGrowthRate: 0.4, halfSaturationConstant: 0.1, yieldCoefficient: 0.5, decayConstant: 0.01, dayDuration: 24,
};

/** Tracks what a test creates on the server and deletes it afterwards. */
export class RunStorageFixture {
  private runIds = new Set<string>();
  private spaces: DG.Project[] = [];
  private subspaces: DG.Project[] = [];
  private favorites: DG.FuncCall[] = [];

  async space(): Promise<DG.Project> {
    const space = await grok.dapi.spaces.createRootSpace(uniqueName('Run storage test'));
    this.spaces.push(space);
    return space;
  }

  async subspace(parent: DG.Project): Promise<DG.Project> {
    const space = await grok.dapi.spaces.id(parent.id).addSubspace(uniqueName('Run storage subspace'));
    this.subspaces.push(space);
    return space;
  }

  /** A DayCalculation call with its inputs; only a called one has the output table, which a save uploads and
   * shares, so tests that don't read outputs save it uncalled. */
  async dayRun({called = false} = {}): Promise<DG.FuncCall> {
    const call = DG.Func.byName(DAY_CALCULATION).prepare(DAY_INPUTS);
    if (called)
      await call.call();
    return call;
  }

  async save(call: DG.FuncCall, options: Parameters<typeof historyUtils.saveRun>[1] = {}): Promise<DG.FuncCall> {
    const saved = await historyUtils.saveRun(call, {newId: true, ...options});
    this.track(saved.id);
    return saved;
  }

  async saveDayRun(name: string, spaceId?: string, replace?: DG.FuncCall[]): Promise<DG.FuncCall> {
    return this.save(await this.dayRun(), {metadata: {name}, target: {spaceId, replace}});
  }

  /** A run saved straight through the platform, as runs were before the author option existed. */
  async saveLegacyRun(name: string): Promise<DG.FuncCall> {
    const call = await this.dayRun({called: true});
    call.newId();
    call.friendlyName = name;
    const saved = await grok.dapi.functions.calls.allPackageVersions().save(call);
    this.track(saved.id);
    return saved;
  }

  track(id: string) {
    this.runIds.add(id);
  }

  favorite(run: DG.FuncCall) {
    this.favorites.push(run);
  }

  async dispose() {
    const errors: string[] = [];
    const attempt = (cleanup: Promise<unknown>) => cleanup.catch((e) => errors.push(String(e?.message ?? e)));
    const all = (cleanups: Promise<unknown>[]) => Promise.all(cleanups.map(attempt));
    await all(this.favorites.map((run) => DG.Favorites.remove(run)));
    // one at a time, subspaces first: deleting spaces at once can deadlock on the server
    for (const space of [...this.subspaces, ...this.spaces])
      await attempt(grok.dapi.spaces.delete(space));
    await all([...this.runIds].map(deleteSavedRun));
    if (errors.length > 0)
      throw new Error(`Cleanup failed: ${errors.join('; ')}`);
  }
}

export async function withRunStorage(fn: (fixture: RunStorageFixture) => Promise<void>) {
  const fixture = new RunStorageFixture();
  try {
    await fn(fixture);
  } finally {
    await fixture.dispose();
  }
}

const dayFunc = () => DG.Func.byName(DAY_CALCULATION);

category('Run storage', () => {
  test('Saving a run as listed keeps its inputs and outputs', () => withRunStorage(async (f) => {
    const name = uniqueName('Listed');
    const deleted = await f.save(await f.dayRun({called: true}), {metadata: {name: `${name} a`}});
    const edited = await f.save(await f.dayRun({called: true}), {metadata: {name: `${name} b`}});
    const listed = await historyUtils.listMyRuns(dayFunc());
    const listedDeleted = listed.find((run) => run.id === deleted.id)!;
    const listedEdited = listed.find((run) => run.id === edited.id)!;

    await historyUtils.deleteRun(listedDeleted);
    historyUtils.setRunMetadata(listedEdited, {name: `${name} renamed`});
    await historyUtils.updateRunMeta(listedEdited);

    for (const id of [deleted.id, edited.id]) {
      const loaded = await historyUtils.loadRun(id);
      expect(loaded.inputs['solutionAdded'], 0.5);
      expect((loaded.outputs['dailyProfile'] as DG.DataFrame).rowCount > 0, true);
    }
    expect((await historyUtils.loadRun(deleted.id)).options['isDeleted'], true);
    expect(getRunName(await historyUtils.loadRun(edited.id)), `${name} renamed`);
  }));

  test('Saving a run loaded for its details keeps its inputs and outputs', () => withRunStorage(async (f) => {
    const name = uniqueName('Details');
    const edited = await f.save(await f.dayRun({called: true}), {metadata: {name: `${name} a`}});
    const deleted = await f.save(await f.dayRun({called: true}), {metadata: {name: `${name} b`}});
    const forEdit = await historyUtils.shallowLoadRun(edited.id);
    historyUtils.setRunMetadata(forEdit, {name: `${name} renamed`, description: 'Edited', tags: ['edited']});
    await historyUtils.updateRunMeta(forEdit);
    await historyUtils.deleteRun(await historyUtils.shallowLoadRun(deleted.id));

    for (const id of [edited.id, deleted.id]) {
      const loaded = await historyUtils.loadRun(id);
      expect(loaded.inputs['solutionAdded'], 0.5);
      expect((loaded.outputs['dailyProfile'] as DG.DataFrame).rowCount > 0, true);
    }
    const reloaded = await historyUtils.loadRun(edited.id);
    expect(getRunName(reloaded), `${name} renamed`);
    expect(reloaded.options['description'], 'Edited');
    expect((await historyUtils.loadRun(deleted.id)).options['isDeleted'], true);
  }));

  test('Saving to My runs stores the name, description and tags', () => withRunStorage(async (f) => {
    const name = uniqueName('My run');
    const saved = await f.save(await f.dayRun(),
      {metadata: {name, description: 'Reference run', tags: ['reference', 'q4']}, target: {}});
    const loaded = await historyUtils.loadRun(saved.id);
    expect(loaded.friendlyName, name);
    expect(loaded.options['description'], 'Reference run');
    expect(JSON.stringify(loaded.options['tags']), JSON.stringify(['reference', 'q4']));
    expect(ids(await historyUtils.listMyRuns(dayFunc())).includes(saved.id), true);
  }));

  test('Saving to a space links the run, and it shows in My runs too', () => withRunStorage(async (f) => {
    const space = await f.space();
    const saved = await f.save(await f.dayRun({called: true}),
      {metadata: {name: uniqueName('Space run')}, target: {spaceId: space.id}});
    expect(ids(await spaceRuns(space, dayFunc())).includes(saved.id), true);
    expect(ids(await historyUtils.listMyRuns(dayFunc())).includes(saved.id), true);
  }));

  test('Replacing in a space swaps the link and keeps the old run', () => withRunStorage(async (f) => {
    const space = await f.space();
    const name = uniqueName('Baseline');
    const old = await f.saveDayRun(name, space.id);
    const conflict = await historyUtils.findNameConflict(` ${name.toUpperCase()} `, space.id);
    expect(JSON.stringify(ids(conflict.runs)), JSON.stringify([old.id]));

    const replacement = await f.saveDayRun(name, space.id, conflict.runs);
    const linked = ids(await spaceRuns(space, dayFunc()));
    expect(linked.includes(replacement.id), true);
    expect(linked.includes(old.id), false);
    const kept = await historyUtils.loadRun(old.id);
    expect(kept.friendlyName, name);
    expect(kept.options['isDeleted'] == null, true);
    expect(ids(await historyUtils.listMyRuns(dayFunc())).includes(old.id), true);
  }));

  test('My runs keeps runs with the same name', () => withRunStorage(async (f) => {
    const name = uniqueName('Same name');
    const first = await f.saveDayRun(name);
    const second = await f.saveDayRun(name);
    const listed = ids(await historyUtils.listMyRuns(dayFunc()));
    expect(listed.includes(first.id) && listed.includes(second.id), true);
  }));

  test('A legacy run saved under a new name gets it as its title too', () => withRunStorage(async (f) => {
    const legacyCall = await f.dayRun();
    legacyCall.options['title'] = uniqueName('Legacy');
    const legacy = await f.save(legacyCall);
    const loaded = await historyUtils.loadRun(legacy.id);
    expect(getRunName(loaded), legacyCall.options['title']);

    const name = uniqueName('Renamed');
    const saved = await f.save(loaded, {metadata: {name}, target: {}});
    const reloaded = await historyUtils.loadRun(saved.id);
    expect(reloaded.options['title'], name);
    expect(reloaded.friendlyName, name);
    expect(getRunName(reloaded), name);
  }));

  test('A soft delete keeps the space links', () => withRunStorage(async (f) => {
    const space = await f.space();
    const saved = await f.saveDayRun(uniqueName('Soft delete'));
    await historyUtils.linkRun(saved.id, space.id);
    await historyUtils.deleteRun(await historyUtils.loadRun(saved.id));
    expect(ids(await historyUtils.listMyRuns(dayFunc())).includes(saved.id), false);
    expect(ids(await spaceRuns(space, dayFunc())).includes(saved.id), true);
  }));

  test('Saving a run deleted from My runs creates a live run', () => withRunStorage(async (f) => {
    const space = await f.space();
    const deleted = await f.saveDayRun(uniqueName('Deleted'), space.id);
    await historyUtils.deleteRun(await historyUtils.loadRun(deleted.id));
    const fromSpace = (await spaceRuns(space, dayFunc())).find((run) => run.id === deleted.id);
    const saved = await f.save(await historyUtils.loadRun(fromSpace!.id), {metadata: {name: uniqueName('Resaved')}});
    expect((await historyUtils.loadRun(saved.id)).options['isDeleted'] == null, true);
    expect(ids(await historyUtils.listMyRuns(dayFunc())).includes(saved.id), true);
  }));

  test('Removing a link keeps the run', () => withRunStorage(async (f) => {
    const space = await f.space();
    const saved = await f.saveDayRun(uniqueName('Unlink'), space.id);
    await historyUtils.unlinkRun(saved.id, space.id);
    expect(ids(await spaceRuns(space, dayFunc())).includes(saved.id), false);
    const loaded = await historyUtils.loadRun(saved.id);
    expect(loaded.options['isDeleted'] == null, true);
    expect(loaded.inputs['solutionAdded'], 0.5);
    expect(ids(await historyUtils.listMyRuns(dayFunc())).includes(saved.id), true);
  }));

  test('Add to space links the same run and keeps it in My runs', () => withRunStorage(async (f) => {
    const space = await f.space();
    const saved = await f.saveDayRun(uniqueName('Add to space'));
    await historyUtils.linkRun(saved.id, space.id);
    expect(JSON.stringify(ids(await spaceRuns(space, dayFunc()))),
      JSON.stringify([saved.id]));
    expect(ids(await historyUtils.listMyRuns(dayFunc())).includes(saved.id), true);
  }));

  test('Add to space can rename the run, keeping its author and time', () => withRunStorage(async (f) => {
    const space = await f.space();
    const saved = await f.saveDayRun(uniqueName('Before'));
    const run = await historyUtils.shallowLoadRun(saved.id);
    const created = run.options[CREATED_AT_OPTION];
    const name = uniqueName('After');
    await addToSpace(run, space, name);

    const linked = await spaceRuns(space, dayFunc());
    expect(JSON.stringify(ids(linked)), JSON.stringify([saved.id]));
    const loaded = await historyUtils.loadRun(saved.id);
    expect(getRunName(loaded), name);
    expect(loaded.options[CREATED_BY_OPTION], DG.User.current().id);
    expect(loaded.options[CREATED_AT_OPTION], created);
    expect(loaded.inputs['solutionAdded'], 0.5);
    expect(ids(await historyUtils.listMyRuns(dayFunc())).includes(saved.id), true);
  }));

  test('The name check finds runs and other items in a space', () => withRunStorage(async (f) => {
    const space = await f.space();
    const name = uniqueName('Feed schedule');
    const first = await f.saveDayRun(name, space.id);
    const second = await f.saveDayRun(name.toUpperCase(), space.id);
    await grok.dapi.spaces.id(space.id).addSubspace('Sensor exports');

    const conflict = await historyUtils.findNameConflict(`  ${name.toLowerCase()} `, space.id);
    expect(JSON.stringify(ids(conflict.runs).sort()), JSON.stringify([first.id, second.id].sort()));
    expect(conflict.other == null, true);
    const folder = await historyUtils.findNameConflict('SENSOR EXPORTS', space.id);
    expect(folder.runs.length, 0);
    expect(folder.other?.entityType, 'Project');

    const replacement = await f.saveDayRun(name, space.id, conflict.runs);
    expect(JSON.stringify(ids(await spaceRuns(space, dayFunc()))),
      JSON.stringify([replacement.id]));
  }));

  test('The name check skips the space\'s own storage connection', () => withRunStorage(async (f) => {
    const space = await f.space();
    const connection = (await historyUtils.listSpace(space.id)).items
      .find((item) => item.entityType === 'DataConnection');
    expect(connection != null, true, 'The space has no storage connection');
    const conflict = await historyUtils.findNameConflict(connection!.friendlyName || connection!.name, space.id);
    expect(conflict.other == null, true);
  }));

  test('Adding several runs to a space adds, skips, replaces and blocks each by its name', () =>
    withRunStorage(async (f) => {
      const space = await f.space();
      const call = await f.dayRun();
      const save = (name: string, spaceId?: string) => f.save(call, {metadata: {name}, target: {spaceId}});
      const name = uniqueName('Several');
      const taken = await save(`${name} taken`, space.id);
      const already = await save(`${name} already`, space.id);
      await grok.dapi.spaces.id(space.id).addSubspace(`${name} folder`);
      const fresh = await save(`${name} fresh`);
      const takenAgain = await save(`${name} TAKEN `);
      // blocked runs are never linked, so they need no saved copy
      const unsaved = (runName: string) => {
        const run = dayFunc().prepare({});
        run.newId();
        run.friendlyName = runName;
        return run;
      };
      const folderNamed = unsaved(`${name} folder`);
      const twins = [unsaved(`${name} twin`), unsaved(`${name} twin`)];

      const plan = planSpaceAdd([fresh, already, takenAgain, folderNamed, ...twins],
        await historyUtils.listSpace(space.id));
      const kind = (run: DG.FuncCall) => plan.find((outcome) => outcome.run.id === run.id)!.kind;
      expect(kind(fresh), 'add');
      expect(kind(already), 'skip');
      expect(kind(takenAgain), 'replace');
      expect(kind(folderNamed), 'blocked');
      expect(twins.map(kind).join(), 'blocked,blocked');

      const {added, failed} = await addRunsToSpace(space, plan);
      expect(failed.length, 0);
      expect(JSON.stringify(ids(added).sort()), JSON.stringify([fresh.id, takenAgain.id].sort()));
      const linked = ids((await historyUtils.listSpace(space.id)).runs);
      expect(JSON.stringify(linked.sort()), JSON.stringify([already.id, fresh.id, takenAgain.id].sort()));
      expect((await historyUtils.shallowLoadRun(taken.id)) != null, true);
    }));

  test('The name check covers only the chosen space', () => withRunStorage(async (f) => {
    const [spaceA, spaceB] = [await f.space(), await f.space()];
    const name = uniqueName('One space');
    await f.saveDayRun(name);
    const inA = await f.saveDayRun(name, spaceA.id);
    expect(JSON.stringify(ids((await historyUtils.findNameConflict(name, spaceA.id)).runs)), JSON.stringify([inA.id]));
    expect((await historyUtils.findNameConflict(name, spaceB.id)).runs.length, 0);
  }));

  test('The run name prefers the legacy title', async () => {
    const call = dayFunc().prepare({});
    call.friendlyName = 'Friendly';
    call.options['title'] = 'Title';
    expect(getRunName(call), 'Title');
    delete call.options['title'];
    expect(getRunName(call), 'Friendly');
  });

  test('Listings of several models keep runs of those models only', () => withRunStorage(async (f) => {
    const space = await f.space();
    await grok.dapi.spaces.id(space.id).addSubspace('Not a run');
    const workflowFunc = DG.Func.byName(BIOREACTOR_WORKFLOW);
    const day = await f.saveDayRun(uniqueName('Several models day'), space.id);
    const workflow = await f.save(workflowFunc.prepare({}),
      {metadata: {name: uniqueName('Several models workflow')}, target: {spaceId: space.id}});
    const both = [workflowFunc, dayFunc()];

    expect(JSON.stringify(ids(await spaceRuns(space, both)).sort()), JSON.stringify([day.id, workflow.id].sort()));
    expect(JSON.stringify(ids(await spaceRuns(space, [workflowFunc]))), JSON.stringify([workflow.id]));
    const mine = ids(await listMyModelRuns(both));
    expect(mine.includes(day.id) && mine.includes(workflow.id), true);
    expect(ids(await listMyModelRuns([workflowFunc])).includes(day.id), false);
  }));

  test('Listings keep runs of the model only', () => withRunStorage(async (f) => {
    const space = await f.space();
    const day = await f.saveDayRun(uniqueName('Day'), space.id);
    const workflowCall = DG.Func.byName(BIOREACTOR_WORKFLOW).prepare({});
    const workflow = await f.save(workflowCall,
      {metadata: {name: uniqueName('Workflow')}, target: {spaceId: space.id}});

    expect(JSON.stringify(ids(await spaceRuns(space, dayFunc()))), JSON.stringify([day.id]));
    expect(JSON.stringify(ids(await spaceRuns(space, DG.Func.byName(BIOREACTOR_WORKFLOW)))),
      JSON.stringify([workflow.id]));
    const mine = await historyUtils.listMyRuns(dayFunc());
    expect(mine.every((run) => run.func.name === dayFunc().name &&
      run.func.package?.name === dayFunc().package.name), true);
  }));

  test('Favorites use the platform favorites', () => withRunStorage(async (f) => {
    const saved = await f.saveDayRun(uniqueName('Favorite'));
    const run = await historyUtils.loadRun(saved.id);
    f.favorite(run);
    await DG.Favorites.add(run);
    expect(ids((await grok.dapi.entities.getFavorites()) as DG.FuncCall[]).includes(run.id), true);
    await DG.Favorites.remove(run);
    expect(ids((await grok.dapi.entities.getFavorites()) as DG.FuncCall[]).includes(run.id), false);
  }));

  test('A workflow saved to a space links only its own run', () => withRunStorage(async (f) => {
    const space = await f.space();
    const result = await runWorkflow(BIOREACTOR_WORKFLOW, '1.0');
    const name = uniqueName('Bioreactor workflow');
    const run = await saveWorkflowRun(BIOREACTOR_WORKFLOW, '1.0', result.state,
      {title: name, target: {spaceId: space.id}});
    try {
      const [state] = await loadInstanceState(run.id);
      const stepIds = (function collect(node: PipelineSerializedState): string[] {
        return isFuncCallSerializedState(node) ?
          (node.funcCallId ? [node.funcCallId] : []) : node.steps.flatMap(collect);
      })(state);
      expect(stepIds.length > 0, true);
      expect(run.friendlyName, name);
      const {runs} = await historyUtils.listSpace(space.id);
      expect(JSON.stringify(ids(runs)), JSON.stringify([run.id]));
    } finally {
      await deleteSavedWorkflow(run.id);
    }
  }));

  test('A save records who saved the run and when', () => withRunStorage(async (f) => {
    const before = Date.now();
    const saved = await f.saveDayRun(uniqueName('Author'));
    const loaded = await historyUtils.loadRun(saved.id);
    expect(loaded.options[CREATED_BY_OPTION], DG.User.current().id);
    expect(getRunCreated(loaded)!.valueOf() >= before - 1000, true);

    loaded.options[CREATED_BY_OPTION] = SOMEONE_ELSE;
    loaded.options[CREATED_AT_OPTION] = LONG_AGO;
    const copy = await historyUtils.loadRun((await f.save(loaded, {metadata: {name: uniqueName('Copy')}})).id);
    expect(copy.options[CREATED_BY_OPTION], DG.User.current().id);
    expect(copy.options[CREATED_AT_OPTION] !== LONG_AGO, true);
  }));

  test('Edit details and delete keep the original author and time', () => withRunStorage(async (f) => {
    const name = uniqueName('Keep author');
    const edited = await f.saveDayRun(`${name} a`);
    const deleted = await f.saveDayRun(`${name} b`);
    for (const id of [edited.id, deleted.id]) {
      const run = await historyUtils.shallowLoadRun(id);
      run.options[CREATED_BY_OPTION] = SOMEONE_ELSE;
      run.options[CREATED_AT_OPTION] = LONG_AGO;
      await historyUtils.updateRunMeta(run);
    }
    const forEdit = await historyUtils.shallowLoadRun(edited.id);
    historyUtils.setRunMetadata(forEdit, {name: `${name} renamed`});
    await historyUtils.updateRunMeta(forEdit);
    await historyUtils.deleteRun(await historyUtils.shallowLoadRun(deleted.id));

    for (const id of [edited.id, deleted.id]) {
      const loaded = await historyUtils.shallowLoadRun(id);
      expect(loaded.options[CREATED_BY_OPTION], SOMEONE_ELSE);
      expect(loaded.options[CREATED_AT_OPTION], LONG_AGO);
    }
    expect(getRunName(await historyUtils.shallowLoadRun(edited.id)), `${name} renamed`);
  }));

  test('Editing a run saved before the author option records its platform author', () => withRunStorage(async (f) => {
    const legacy = await f.saveLegacyRun(uniqueName('Legacy author'));
    const listed = (await historyUtils.listMyRuns(dayFunc())).find((run) => run.id === legacy.id)!;
    expect(listed.options[CREATED_BY_OPTION] == null, true);
    const started = listed.started.valueOf();
    historyUtils.setRunMetadata(listed, {name: uniqueName('Legacy renamed')});
    await historyUtils.updateRunMeta(listed);

    const loaded = await historyUtils.shallowLoadRun(legacy.id);
    expect(loaded.options[CREATED_BY_OPTION], DG.User.current().id);
    expect(getRunCreated(loaded)!.valueOf(), started);
  }));

  test('My runs lists runs by their original author', () => withRunStorage(async (f) => {
    const name = uniqueName('Original author');
    const mine = await f.saveDayRun(`${name} mine`);
    const legacy = await f.saveLegacyRun(`${name} legacy`);
    const theirs = await historyUtils.shallowLoadRun((await f.saveDayRun(`${name} theirs`)).id);
    theirs.options[CREATED_BY_OPTION] = SOMEONE_ELSE;
    await historyUtils.updateRunMeta(theirs);

    const listed = ids(await historyUtils.listMyRuns(dayFunc()));
    expect(listed.includes(mine.id), true);
    expect(listed.includes(legacy.id), true);
    expect(listed.includes(theirs.id), false);
    expect(listed.length, new Set(listed).size);
  }));

  test('A run by someone else in an editable space keeps its author when renamed', () => withRunStorage(async (f) => {
    const space = await f.space();
    const saved = await f.saveDayRun(uniqueName('Theirs'), space.id);
    const run = await historyUtils.shallowLoadRun(saved.id);
    run.options[CREATED_BY_OPTION] = SOMEONE_ELSE;
    await historyUtils.updateRunMeta(run);

    const forEdit = await historyUtils.shallowLoadRun(saved.id);
    expect(historyUtils.canEditRunDetails(getRunAuthorId(forEdit), await canEditSpace(space)), true);
    const name = uniqueName('Renamed theirs');
    historyUtils.setRunMetadata(forEdit, {name});
    await historyUtils.updateRunMeta(forEdit);

    const loaded = await historyUtils.shallowLoadRun(saved.id);
    expect(getRunName(loaded), name);
    expect(loaded.options[CREATED_BY_OPTION], SOMEONE_ELSE);
    expect(ids(await spaceRuns(space, dayFunc())).includes(saved.id), true);
  }));

  test('My runs lists runs of a model taken from the model catalog', () => withRunStorage(async (f) => {
    const workflow = DG.Func.byName(BIOREACTOR_WORKFLOW);
    const fromCatalog = (await grok.dapi.functions.filter(getModelFilter()).list())
      .find((func) => func.nqName === workflow.nqName);
    expect(fromCatalog != null, true, `${BIOREACTOR_WORKFLOW} is not in the model catalog`);
    const saved = await f.save(workflow.prepare({}), {metadata: {name: uniqueName('Catalog model')}});
    expect(ids(await historyUtils.listMyRuns(fromCatalog!)).includes(saved.id), true);
  }));

  test('Only saves through the run storage options tell open run browsers', () => withRunStorage(async (f) => {
    const changes: RunsChange[] = [];
    const sub = runsChanged$.subscribe((change) => changes.push(change));
    try {
      const old = await f.dayRun();
      old.newId();
      f.track((await historyUtils.saveRun(old)).id);
      await f.save(await f.dayRun(), {nested: true});
      await f.saveDayRun(uniqueName('Told'));
      await awaitCheck(() => changes.length > 0, 'the last save should reach run browsers', 3000);
      expect(JSON.stringify(changes), JSON.stringify([{spaceIds: []}]));
    } finally {
      sub.unsubscribe();
    }
  }));

  test('Workflow step runs keep the default name and stay out of My runs unless saved on their own', () =>
    withRunStorage(async (f) => {
      const result = await runWorkflow(BIOREACTOR_WORKFLOW, '1.0');
      const name = uniqueName('Bioreactor workflow steps');
      const changes: RunsChange[] = [];
      const sub = runsChanged$.subscribe((change) => changes.push(change));
      const run = await saveWorkflowRun(BIOREACTOR_WORKFLOW, '1.0', result.state, {title: name});
      sub.unsubscribe();
      try {
        const [state] = await loadInstanceState(run.id);
        const stepId = (function first(node: PipelineSerializedState): string | undefined {
          return isFuncCallSerializedState(node) ? node.funcCallId : node.steps.map(first).find((id) => !!id);
        })(state);
        expect(stepId != null, true);
        const step = await historyUtils.shallowLoadRun(stepId!);
        expect(getRunName(step) !== '' && getRunName(step) !== name, true);
        expect(step.options[CREATED_BY_OPTION], DG.User.current().id);
        expect(step.options[NESTED_ENTITY_OPTION], 'true');
        expect(JSON.stringify(changes), JSON.stringify([{spaceIds: []}]));
        const stepFunc = DG.Func.byName(step.func.nqName);
        expect(ids(await historyUtils.listMyRuns(stepFunc)).includes(step.id), false);
        expect(placeRunFilter(stepFunc)(step), false);
        expect(ids(await historyUtils.listMyRuns(DG.Func.byName(BIOREACTOR_WORKFLOW))).includes(run.id), true);

        const own = await f.save(await historyUtils.loadRun(step.id),
          {metadata: {name: uniqueName('Step on its own')}});
        expect((await historyUtils.shallowLoadRun(own.id)).options[NESTED_ENTITY_OPTION] == null, true);
        expect(ids(await historyUtils.listMyRuns(stepFunc)).includes(own.id), true);
      } finally {
        await deleteSavedWorkflow(run.id);
      }
    }));

  test('Space runs get their authors, older runs their platform author', () => withRunStorage(async (f) => {
    const space = await f.space();
    const name = uniqueName('Space author');
    const mine = await f.saveDayRun(`${name} mine`, space.id);
    const legacy = await f.saveLegacyRun(`${name} legacy`);
    await historyUtils.linkRun(legacy.id, space.id);
    const unknown = await historyUtils.shallowLoadRun((await f.saveDayRun(`${name} unknown`, space.id)).id);
    unknown.options[CREATED_BY_OPTION] = SOMEONE_ELSE;
    await historyUtils.updateRunMeta(unknown);

    const listed = await spaceRuns(space, dayFunc());
    expect(listed.length, 3);
    const authors = await loadRunAuthors(listed);
    const me = DG.User.current();
    expect(authors.get(mine.id)?.id, me.id);
    expect(authors.get(mine.id)?.name, me.friendlyName);
    expect(authors.get(legacy.id)?.id, me.id);
    expect(authors.has(unknown.id), false);
  }));

  test('The run version is the workflow version, else the package version', () => withRunStorage(async (f) => {
    const result = await runWorkflow(BIOREACTOR_WORKFLOW, '1.0');
    const workflow = await saveWorkflowRun(BIOREACTOR_WORKFLOW, '1.0', result.state,
      {title: uniqueName('Workflow version')});
    try {
      expect(runVersion(await historyUtils.shallowLoadRun(workflow.id)), '1.0');
    } finally {
      await deleteSavedWorkflow(workflow.id);
    }
    const day = await f.saveDayRun(uniqueName('Package version'));
    expect(runVersion(await historyUtils.shallowLoadRun(day.id)), dayFunc().package.version);
  }));

  test('Edit details is allowed for the author or with Edit on the space', async () => {
    const me = DG.User.current().id;
    const someoneElse = '00000000-0000-0000-0000-000000000001';
    expect(historyUtils.canEditRunDetails(me, false), true);
    expect(historyUtils.canEditRunDetails(someoneElse, true), true);
    expect(historyUtils.canEditRunDetails(someoneElse, false), false);
    expect(historyUtils.canEditRunDetails(undefined, false), false);
  });

  test('Add to space under a new name renames the saved run, not the one it was given', () =>
    withRunStorage(async (f) => {
      const space = await f.space();
      const name = uniqueName('Listed name');
      const given = await historyUtils.shallowLoadRun((await f.saveDayRun(name)).id);
      const renamed = uniqueName('Name in space');
      await addToSpace(given, space, renamed);
      expect(getRunName(given), name);
      expect(getRunName(await historyUtils.shallowLoadRun(given.id)), renamed);
      expect(ids(await spaceRuns(space)).includes(given.id), true);
    }));

  test('Add to space replaces runs of the same name, never the run itself', () => withRunStorage(async (f) => {
    const space = await f.space();
    const name = uniqueName('Replaced by add');
    const old = await f.saveDayRun(name, space.id);
    const run = await f.saveDayRun(name);
    await addToSpace(run, space, name, [old, run]);
    expect(JSON.stringify(ids(await spaceRuns(space))), JSON.stringify([run.id]));
  }));

  test('Space listings are newest first', () => withRunStorage(async (f) => {
    const space = await f.space();
    const first = await f.saveDayRun(uniqueName('Older'), space.id);
    const second = await f.saveDayRun(uniqueName('Newer'), space.id);
    expect(JSON.stringify(ids(await spaceRuns(space, dayFunc()))), JSON.stringify([second.id, first.id]));
  }));

  test('Deleting a run saved before the author option records its platform author', () => withRunStorage(async (f) => {
    const space = await f.space();
    const legacy = await f.saveLegacyRun(uniqueName('Legacy delete'));
    await historyUtils.linkRun(legacy.id, space.id);
    // space listings come without the platform author, so the delete has to look it up
    const listed = (await spaceRuns(space)).find((run) => run.id === legacy.id)!;
    const started = listed.started.valueOf();
    await historyUtils.deleteRun(listed);

    const loaded = await historyUtils.shallowLoadRun(legacy.id);
    expect(loaded.options['isDeleted'], true);
    expect(loaded.options[CREATED_BY_OPTION], DG.User.current().id);
    expect(getRunCreated(loaded)!.valueOf(), started);
  }));

  test('Delete and remove try every run and report the ones that failed', () => withRunStorage(async (f) => {
    const runs = [await f.saveDayRun(uniqueName('Deleted A')), await f.saveDayRun(uniqueName('Deleted B'))];
    const deleted = await deleteOrRemoveRuns(runs, null);
    expect(JSON.stringify([ids(deleted.done), deleted.failed.length]), JSON.stringify([ids(runs), 0]));
    for (const run of runs)
      expect((await historyUtils.shallowLoadRun(run.id)).options['isDeleted'], true);

    // a space deleted while it is browsed can't be changed any more
    const space = await grok.dapi.spaces.createRootSpace(uniqueName('Gone'));
    for (const run of runs)
      await historyUtils.linkRun(run.id, space.id);
    await grok.dapi.spaces.delete(space);
    const removed = await deleteOrRemoveRuns(runs, space);
    expect(removed.done.length, 0);
    expect(JSON.stringify(removed.failed.map(({run}) => run.id)), JSON.stringify(ids(runs)));
    expect(removed.failed.every(({error}) => error.length > 0), true);
  }));

  test('Edit details refuses a name taken in the space and saves another', () => withRunStorage(async (f) => {
    const space = await f.space();
    const taken = uniqueName('Taken');
    await f.saveDayRun(taken, space.id);
    const saved = await f.saveDayRun(uniqueName('Edited'), space.id);
    let done = 0;
    const root = ui.div();
    document.body.append(root);
    const given = await historyUtils.shallowLoadRun(saved.id);
    const givenName = getRunName(given);
    let savedRun: DG.FuncCall | undefined;
    const app = Vue.createApp(RunDetailsForm, {run: given, place: space, onDone: (run?: DG.FuncCall) => {
      done++;
      savedRun = run;
    }});
    app.mount(root);
    try {
      const nameInput = root.querySelector('input')!;
      const save = (name: string) => {
        nameInput.value = name;
        nameInput.dispatchEvent(new Event('input'));
        Array.from(root.querySelectorAll('button')).find((button) => button.textContent?.trim() === 'Save details')!
          .click();
      };
      save(taken.toUpperCase());
      await awaitCheck(() => !!root.querySelector('.c2-runs-error')?.textContent?.includes('is already in'),
        'a taken name should be refused', 10000);
      expect(done, 0);
      const renamed = uniqueName('Renamed');
      save(renamed);
      await awaitCheck(() => done === 1, 'a free name should save', 10000);
      expect(getRunName(await historyUtils.shallowLoadRun(saved.id)), renamed);
      expect(getRunName(given), givenName);
      expect(getRunName(savedRun!), renamed);
    } finally {
      app.unmount();
      root.remove();
    }
  }));

  test('Edit details saves once on a double click', () => withRunStorage(async (f) => {
    const saved = await f.saveDayRun(uniqueName('Double click'));
    let done = 0;
    const root = ui.div();
    document.body.append(root);
    const app = Vue.createApp(RunDetailsForm, {run: await historyUtils.shallowLoadRun(saved.id), onDone: () => done++});
    app.mount(root);
    try {
      const button = Array.from(root.querySelectorAll('button'))
        .find((item) => item.textContent?.trim() === 'Save details')!;
      button.click();
      button.click();
      await awaitCheck(() => done > 0, 'the details should save', 10000);
      // a request sent now comes back after a second save that started with the first one
      await historyUtils.shallowLoadRun(saved.id);
      expect(done, 1);
    } finally {
      app.unmount();
      root.remove();
    }
  }));

  test('Run names, details and tags default sensibly', async () => {
    const func = dayFunc();
    expect(defaultRunName(func, 'Loaded', 'Title link'), 'Loaded');
    expect(defaultRunName(func, '', 'Title link'), 'Title link');
    expect(defaultRunName(func).startsWith(`${modelLabel(func)} `), true);
    expect(JSON.stringify(parseTags(' a, b,, a ,')), JSON.stringify(['a', 'b']));
    const call = func.prepare({});
    call.options['tags'] = 'not a list';
    const {description, tags} = historyUtils.getRunMetadata(call);
    expect(JSON.stringify({description, tags}), JSON.stringify({description: '', tags: []}));
    expect((await historyUtils.listMyRuns([])).length, 0);
    expect(await canEditSpace(null), false);
  });
});
