import * as DG from 'datagrok-api/dg';
import {category, test} from '@datagrok-libraries/utils/src/test';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {callHandler} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/utils';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';

category('ComputeUtils: Driver state tree persistence', async () => {
  test('Load and save simple config', async () => {
    const config = await callHandler<PipelineConfiguration>('LibTests:MockProvider1', {}).toPromise();
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromConfig(pconf);
    await tree.init().toPromise();
    const sc = tree.toSerializedState({ disableNodesUUID: true, disableCallsUUID: true });
    const id = await tree.save().toPromise();
    const loadedTree = await StateTree.load(id, pconf).toPromise();
    await loadedTree.init().toPromise();
    const lc = loadedTree.toSerializedState({ disableNodesUUID: true, disableCallsUUID: true });
    expectDeepEqual(lc, sc);
  });

  test('Load nested pipeline', async () => {
    // create and save nested pipeline
    const config = await callHandler<PipelineConfiguration>('LibTests:MockProvider2', {version: '1.0'}).toPromise();
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromConfig(pconf);
    await tree.init().toPromise();
    const sc = tree.toSerializedState({ disableNodesUUID: true, disableCallsUUID: true });
    const dbId = await tree.save().toPromise();
    // create outer pipeline
    const outerConfig = await callHandler<PipelineConfiguration>('LibTests:MockProvider3', {version: '1.0'}).toPromise();
    const outerPconf = await getProcessedConfig(outerConfig);
    const outerTree = StateTree.fromConfig(outerPconf);
    await outerTree.init().toPromise();
    // load nested tree into outer
    const root = outerTree.getItem([]);
    await outerTree.loadSubTree(root.uuid, dbId, 'pipelinePar', 1).toPromise();
    const loadedTree = outerTree.getNode([{idx: 1}]);
    const lc = StateTree.toStateRec(loadedTree, true, { disableNodesUUID: true, disableCallsUUID: true });
    expectDeepEqual(lc, sc);
  });

  test('Save nested pipeline', async () => {
    // create outer pipeline
    const outerConfig = await callHandler<PipelineConfiguration>('LibTests:MockProvider3', {version: '1.0'}).toPromise();
    const outerPconf = await getProcessedConfig(outerConfig);
    const outerTree = StateTree.fromConfig(outerPconf);
    await outerTree.init().toPromise();
    const nestedRoot = outerTree.getNode([{idx: 0}]);
    const sc = StateTree.toStateRec(nestedRoot, true, { disableNodesUUID: true, disableCallsUUID: true });
    // save nested pipeline
    const dbId = await outerTree.save(nestedRoot.getItem().uuid).toPromise();
    const config = await callHandler<PipelineConfiguration>('LibTests:MockProvider2', {version: '1.0'}).toPromise();
    const pconf = await getProcessedConfig(config);
    const loadedTree = await StateTree.load(dbId, pconf).toPromise();
    await loadedTree.init().toPromise();
    const lc = loadedTree.toSerializedState({ disableNodesUUID: true, disableCallsUUID: true });
    expectDeepEqual(lc, sc);
  });
});
