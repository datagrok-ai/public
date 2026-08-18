import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, category, expect, test} from '@datagrok-libraries/test/src/test';
import {historyUtils} from '@datagrok-libraries/compute-utils';
import {
  CONFIG_PATH, loadInstanceState,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/funccall-utils';
import {cloneRun} from '../service/clone-service';
import {OPT_FROZEN, OPT_PARENT_CALL_ID, OPT_PUBLICATION_ID} from '../domain/constants';
import {cleanupTestPrograms, makeSavedRun, makeTestProgram} from './fixtures';

category('ArtifactAlignment: clone', () => {
  test('deep clone rewrites child ids and freezes the whole tree readonly', async () => {
    const run = await makeSavedRun(4);
    const publicationId = crypto.randomUUID();
    const clone = await cloneRun(run.metaCallId, {publicationId, title: 'Frozen v1'});

    expect(clone.metaCall.id === run.metaCallId, false);
    const newChildId = clone.childIdMap[run.stepCallId];
    expect(newChildId != null, true, 'the step call must be cloned');
    expect(newChildId === run.stepCallId, false);

    const [config, metaCall] = await loadInstanceState(clone.metaCall.id);
    expect(metaCall.options[OPT_FROZEN], 'true');
    expect(metaCall.options[OPT_PUBLICATION_ID], publicationId);
    const state = config as any;
    expect(state.isReadonly, true, 'the root node must be readonly');
    expect(state.steps[0].isReadonly, true, 'step nodes must be readonly');
    expect(state.steps[0].funcCallId, newChildId);
  });

  test('the frozen child run is loadable and complete, and its source is untouched', async () => {
    const run = await makeSavedRun(5);
    const clone = await cloneRun(run.metaCallId, {publicationId: crypto.randomUUID()});
    const frozenChild = await historyUtils.loadRun(clone.childIdMap[run.stepCallId]);
    expect(frozenChild.options[OPT_PARENT_CALL_ID], clone.metaCall.id);
    const df: DG.DataFrame = frozenChild.outputs['result'];
    expect(df.rowCount, 5, 'the cloned run must carry the output dataframe');
    const source = await historyUtils.loadRun(run.stepCallId);
    expect((source.outputs['result'] as DG.DataFrame).rowCount, 5);
    expect(source.options[OPT_FROZEN] == null, true, 'the source run must not be stamped');
  });

  test('publish clone grants the program viewers group on the frozen dataframes', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun(2);
    const clone = await cloneRun(run.metaCallId,
      {publicationId: crypto.randomUUID(), audience: program.viewers});
    const frozenChild = await historyUtils.loadRun(clone.childIdMap[run.stepCallId]);
    const tableInfo = (frozenChild.outputs['result'] as DG.DataFrame).getTableInfo();
    const perms: any = await grok.dapi.permissions.get(tableInfo);
    const grantedIds: string[] = [];
    for (const key of Object.keys(perms ?? {})) {
      for (const g of perms[key] ?? [])
        grantedIds.push(g?.id);
    }
    expect(grantedIds.includes(program.viewers.id), true,
      `viewers group not among grantees: ${grantedIds.join(', ')}`);
  });

  test('two publications of the same source produce independent frozen copies', async () => {
    const run = await makeSavedRun();
    const c1 = await cloneRun(run.metaCallId, {publicationId: crypto.randomUUID()});
    const c2 = await cloneRun(run.metaCallId, {publicationId: crypto.randomUUID()});
    expect(c1.metaCall.id === c2.metaCall.id, false);
    expect(c1.childIdMap[run.stepCallId] === c2.childIdMap[run.stepCallId], false);
  });

  test('frozen copy content matches the source config structurally', async () => {
    const run = await makeSavedRun();
    const clone = await cloneRun(run.metaCallId, {publicationId: crypto.randomUUID()});
    const source = await historyUtils.loadRun(run.metaCallId);
    const frozen = await historyUtils.loadRun(clone.metaCall.id);
    const strip = (s: string) => JSON.stringify(JSON.parse(s), (k, v) =>
      (k === 'funcCallId' || k === 'uuid' || k === 'isReadonly') ? undefined : v);
    expect(strip(frozen.options[CONFIG_PATH]), strip(source.options[CONFIG_PATH]),
      'configs must match after masking ids and readonly flags');
  });

  after(async () => {
    await cleanupTestPrograms();
  });
}, {timeout: 120000});
