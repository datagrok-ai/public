import * as DG from 'datagrok-api/dg';
import {category, test} from '@datagrok-libraries/utils/src/test';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {callHandler} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/utils';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {PipelineStateStatic, StepFunCallStateBase} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';

category('ComputeUtils: Driver state tree readonly', async () => {
  test('Pipeline load readonly', async () => {
    const config = await callHandler<PipelineConfiguration>('LibTests:MockProvider1', {}).toPromise();
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const sc = tree.toSerializedState({disableNodesUUID: true, disableCallsUUID: true});
    const metaCall = await tree.save().toPromise();
    const loadedTree = await StateTree.load({dbId: metaCall!.id, config: pconf, isReadonly: true}).toPromise();
    await loadedTree.init().toPromise();
    const lc = loadedTree.toSerializedState({disableNodesUUID: true, disableCallsUUID: true});
    sc.isReadonly = true;
    (sc as PipelineStateStatic<StepFunCallStateBase, {}>).steps[0].isReadonly = true;
    (sc as PipelineStateStatic<StepFunCallStateBase, {}>).steps[1].isReadonly = true;
    expectDeepEqual(lc, sc);
  });

  test('Pipeline load readonly subtree', async () => {
    // create and save nested pipeline
    const config = await callHandler<PipelineConfiguration>('LibTests:MockProvider2', {version: '1.0'}).toPromise();
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const sc = tree.toSerializedState({disableNodesUUID: true, disableCallsUUID: true});
    const metaCall = await tree.save().toPromise();
    // create outer pipeline
    const outerConfig = await callHandler<PipelineConfiguration>('LibTests:MockProvider3', {version: '1.0'}).toPromise();
    const outerPconf = await getProcessedConfig(outerConfig);
    const outerTree = StateTree.fromPipelineConfig({config: outerPconf});
    await outerTree.init().toPromise();
    // load nested tree into outer
    const root = outerTree.nodeTree.getItem([]);
    await outerTree.loadSubTree(root.uuid, metaCall!.id, 'pipelinePar', 1, true).toPromise();
    const loadedTree = outerTree.nodeTree.getNode([{idx: 1}]);
    const lc = StateTree.toStateRec(loadedTree, true, {disableNodesUUID: true, disableCallsUUID: true});
    sc.isReadonly = true;
    (sc as PipelineStateStatic<StepFunCallStateBase, {}>).steps[0].isReadonly = true;
    (sc as PipelineStateStatic<StepFunCallStateBase, {}>).steps[1].isReadonly = true;
    ((sc as PipelineStateStatic<StepFunCallStateBase, {}>).steps[1] as PipelineStateStatic<StepFunCallStateBase, {}>).steps[0].isReadonly = true;
    ((sc as PipelineStateStatic<StepFunCallStateBase, {}>).steps[1] as PipelineStateStatic<StepFunCallStateBase, {}>).steps[1].isReadonly = true;
    expectDeepEqual(lc, sc);
  });

  test('Pipeline preserve readonly subtree', async () => {
    // create and save nested pipeline
    const config = await callHandler<PipelineConfiguration>('LibTests:MockProvider2', {version: '1.0'}).toPromise();
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const metaCall = await tree.save().toPromise();
    // create outer pipeline
    const outerConfig = await callHandler<PipelineConfiguration>('LibTests:MockProvider3', {version: '1.0'}).toPromise();
    const outerPconf = await getProcessedConfig(outerConfig);
    const outerTree = StateTree.fromPipelineConfig({config: outerPconf});
    await outerTree.init().toPromise();
    // load nested tree into outer
    const root = outerTree.nodeTree.getItem([]);
    await outerTree.loadSubTree(root.uuid, metaCall!.id, 'pipelinePar', 1, true).toPromise();
    const sc = outerTree.toSerializedState({disableNodesUUID: true, disableCallsUUID: true});
    // save outer with inner
    const metaCallOuter = await outerTree.save().toPromise();
    // load outer with inner
    const loadedTree = await StateTree.load({dbId: metaCallOuter!.id, config: outerPconf}).toPromise();
    const lc = loadedTree.toSerializedState({disableNodesUUID: true, disableCallsUUID: true});
    expectDeepEqual(lc, sc);
  });

});
