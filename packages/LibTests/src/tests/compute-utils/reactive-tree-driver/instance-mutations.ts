import * as DG from 'datagrok-api/dg';
import {category, test} from '@datagrok-libraries/utils/src/test';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {callHandler} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/utils';
import {snapshotCompare} from '../../../test-utils';

category('ComputeUtils: Driver state tree mutations', async () => {
  test('Pipeline append step subtree', async () => {
    const config = await callHandler<PipelineConfiguration>('LibTests:MockProvider2', {}).toPromise();
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const root = tree.nodeTree.getItem([]);
    const puuid = root.uuid;
    await tree.addSubTree(puuid, 'stepAdd', 2).toPromise();
    const state = tree.toSerializedState({disableNodesUUID: true});
    await snapshotCompare(state, 'Pipeline append step subtree');
  });

  test('Pipeline append static pipeline subtree', async () => {
    const config = await callHandler<PipelineConfiguration>('LibTests:MockProvider2', {}).toPromise();
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const root = tree.nodeTree.getItem([]);
    const puuid = root.uuid;
    await tree.addSubTree(puuid, 'pipeline1', 2).toPromise();
    const state = tree.toSerializedState({disableNodesUUID: true});
    await snapshotCompare(state, 'Pipeline append static pipeline subtree');
  });

  test('Pipeline append dynamic pipeline subtree', async () => {
    const config = await callHandler<PipelineConfiguration>('LibTests:MockProvider3', {}).toPromise();
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const root = tree.nodeTree.getItem([]);
    const puuid = root.uuid;
    await tree.addSubTree(puuid, 'pipelinePar', 0).toPromise();
    const state = tree.toSerializedState({disableNodesUUID: true});
    await snapshotCompare(state, 'Pipeline append dynamic pipeline subtree');
  });

  test('Pipeline insert subtree', async () => {
    const config = await callHandler<PipelineConfiguration>('LibTests:MockProvider2', {}).toPromise();
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const root = tree.nodeTree.getItem([]);
    const puuid = root.uuid;
    await tree.addSubTree(puuid, 'pipeline1', 0).toPromise();
    const state = tree.toSerializedState({disableNodesUUID: true});
    await snapshotCompare(state, 'Pipeline insert subtree');
  });

  test('Pipeline remove subtree', async () => {
    const config = await callHandler<PipelineConfiguration>('LibTests:MockProvider2', {}).toPromise();
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const pnode = tree.nodeTree.getItem([{idx: 0}]);
    await tree.removeSubtree(pnode.uuid).toPromise();
    const state = tree.toSerializedState({disableNodesUUID: true});
    await snapshotCompare(state, 'Pipeline remove subtree');
  });

  test('Pipeline move subtree', async () => {
    const config = await callHandler<PipelineConfiguration>('LibTests:MockProvider2', {}).toPromise();
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const pnode = tree.nodeTree.getItem([{idx: 0}]);
    await tree.moveSubtree(pnode.uuid, 1).toPromise();
    const state = tree.toSerializedState({disableNodesUUID: true});
    await snapshotCompare(state, 'Pipeline move subtree');
  });
});
