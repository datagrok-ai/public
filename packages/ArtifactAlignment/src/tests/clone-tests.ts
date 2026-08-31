import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, category, expect, test} from '@datagrok-libraries/test/src/test';
import {historyUtils} from '@datagrok-libraries/compute-utils';
import {
  CONFIG_PATH, loadInstanceState,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/funccall-utils';
import {cloneRun} from '../service/clone-service';
import {OPT_FROZEN, OPT_PUBLICATION_ID} from '../domain/constants';
import {cleanupTestPrograms, makeSavedRun, makeTestProgram, trackCalls} from './fixtures';

async function trackedClone(...args: Parameters<typeof cloneRun>): ReturnType<typeof cloneRun> {
  const clone = await cloneRun(...args);
  trackCalls(clone.metaCall.id);
  return clone;
}

category('ArtifactAlignment: clone', () => {
  test('workflow freeze saves a stamped meta call over the shared config', async () => {
    const run = await makeSavedRun(4);
    const publicationId = crypto.randomUUID();
    const clone = await trackedClone(run.metaCallId, {publicationId, title: 'Frozen v1'});

    expect(clone.metaCall.id === run.metaCallId, false);
    const [config, metaCall] = await loadInstanceState(clone.metaCall.id);
    expect(metaCall.options[OPT_FROZEN], 'true');
    expect(metaCall.options[OPT_PUBLICATION_ID], publicationId);
    const state = config as any;
    expect(state.steps[0].funcCallId, run.stepCallId, 'step calls are shared, not copied');
  });

  test('the shared step run stays loadable, complete, and unstamped', async () => {
    const run = await makeSavedRun(5);
    await trackedClone(run.metaCallId, {publicationId: crypto.randomUUID()});
    const step = await historyUtils.loadRun(run.stepCallId);
    expect((step.outputs['result'] as DG.DataFrame).rowCount, 5);
    expect(step.options[OPT_FROZEN] == null, true, 'the source step must not be stamped');
    const source = await historyUtils.loadRun(run.metaCallId);
    expect(source.options[OPT_FROZEN] == null, true, 'the source meta call must not be stamped');
  });

  test('function-run freeze grants the program viewers group on its dataframes', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun(2);
    const clone = await trackedClone(run.stepCallId,
      {publicationId: crypto.randomUUID(), audience: program.viewers});
    const frozen = await historyUtils.loadRun(clone.metaCall.id);
    const tableInfo = (frozen.outputs['result'] as DG.DataFrame).getTableInfo();
    const perms: any = await grok.dapi.permissions.get(tableInfo);
    const grantedIds: string[] = [];
    for (const key of Object.keys(perms ?? {})) {
      for (const g of perms[key] ?? [])
        grantedIds.push(g?.id);
    }
    expect(grantedIds.includes(program.viewers.id), true,
      `viewers group not among grantees: ${grantedIds.join(', ')}`);
  });

  test('two publications of the same source produce independent meta calls', async () => {
    const run = await makeSavedRun();
    const c1 = await trackedClone(run.metaCallId, {publicationId: crypto.randomUUID()});
    const c2 = await trackedClone(run.metaCallId, {publicationId: crypto.randomUUID()});
    expect(c1.metaCall.id === c2.metaCall.id, false);
  });

  test('frozen copy content matches the source config structurally', async () => {
    const run = await makeSavedRun();
    const clone = await trackedClone(run.metaCallId, {publicationId: crypto.randomUUID()});
    const source = await historyUtils.loadRun(run.metaCallId);
    const frozen = await historyUtils.loadRun(clone.metaCall.id);
    const normalize = (s: string) => JSON.stringify(JSON.parse(s));
    expect(normalize(frozen.options[CONFIG_PATH]), normalize(source.options[CONFIG_PATH]),
      'configs must match, including the shared step ids');
  });

  after(async () => {
    await cleanupTestPrograms();
  });
}, {timeout: 120000});
