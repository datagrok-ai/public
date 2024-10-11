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
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const sc = tree.toSerializedState({disableNodesUUID: true, disableCallsUUID: true});
    const metaCall = await tree.save().toPromise();
    const loadedTree = await StateTree.load({dbId: metaCall!.id, config: pconf}).toPromise();
    await loadedTree.init().toPromise();
    const lc = loadedTree.toSerializedState({disableNodesUUID: true, disableCallsUUID: true});
    expectDeepEqual(lc, sc);
  });

  test('Save creates new uuids', async () => {
    const config = await callHandler<PipelineConfiguration>('LibTests:MockProvider1', {}).toPromise();
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const sc = tree.toSerializedState();
    const metaCall = await tree.save().toPromise();
    const loadedTree = await StateTree.load({dbId: metaCall!.id, config: pconf}).toPromise();
    await loadedTree.init().toPromise();
    const lc = loadedTree.toSerializedState();

    expectDeepEqual(!!sc.uuid, true, {prefix: 'pipeline created uuid'});
    expectDeepEqual(!!lc.uuid, true, {prefix: 'pipeline saved uuid'});

    expectDeepEqual(!!(sc as any).steps[0].uuid, true, {prefix: 'pipeline step1 created uuid'});
    expectDeepEqual(!!(lc as any).steps[0].uuid, true, {prefix: 'pipeline step1 saved uuid'});
    expectDeepEqual(!!(sc as any).steps[1].uuid, true, {prefix: 'pipeline step2 created uuid'});
    expectDeepEqual(!!(lc as any).steps[1].uuid, true, {prefix: 'pipeline step2 saved uuid'});
    expectDeepEqual(!!(sc as any).steps[0].funcCallId, true, {prefix: 'pipeline step1 created funcCall uuid'});
    expectDeepEqual(!!(lc as any).steps[0].funcCallId, true, {prefix: 'pipeline step1 saved funcCall uuid'});
    expectDeepEqual(!!(sc as any).steps[1].funcCallId, true, {prefix: 'pipeline step2 created funcCall uuid'});
    expectDeepEqual(!!(lc as any).steps[1].funcCallId, true, {prefix: 'pipeline step2 saved funcCall uuid'});

    expectDeepEqual(sc.uuid !== lc.uuid, true, {prefix: 'pipeline different node uuid'});
    expectDeepEqual((sc as any).steps[0].uuid !== (lc as any).steps[0].uuid, true, {prefix: 'pipeline step1 different node uuid'});
    expectDeepEqual((sc as any).steps[1].uuid !== (lc as any).steps[1].uuid, true, {prefix: 'pipeline step2 different node uuid'});
    expectDeepEqual((sc as any).steps[0].funcCallId !== (lc as any).steps[0].funcCallId, true, {prefix: 'pipeline step1 different funcCall uuid'});
    expectDeepEqual((sc as any).steps[1].funcCallId !== (lc as any).steps[1].funcCallId, true, {prefix: 'pipeline step2 different funcCall uuid'});
  });

  test('Load nested pipeline', async () => {
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
    await outerTree.loadSubTree(root.uuid, metaCall!.id, 'pipelinePar', 1, false).toPromise();
    const loadedTree = outerTree.nodeTree.getNode([{idx: 1}]);
    const lc = StateTree.toStateRec(loadedTree, true, {disableNodesUUID: true, disableCallsUUID: true});
    expectDeepEqual(lc, sc);
  });

  test('Save nested pipeline', async () => {
    // create outer pipeline
    const outerConfig = await callHandler<PipelineConfiguration>('LibTests:MockProvider3', {version: '1.0'}).toPromise();
    const outerPconf = await getProcessedConfig(outerConfig);
    const outerTree = StateTree.fromPipelineConfig({config: outerPconf});
    await outerTree.init().toPromise();
    const nestedRoot = outerTree.nodeTree.getNode([{idx: 0}]);
    const sc = StateTree.toStateRec(nestedRoot, true, {disableNodesUUID: true, disableCallsUUID: true});
    // save nested pipeline
    const metaCall = await outerTree.save(nestedRoot.getItem().uuid).toPromise();
    const config = await callHandler<PipelineConfiguration>('LibTests:MockProvider2', {version: '1.0'}).toPromise();
    const pconf = await getProcessedConfig(config);
    const loadedTree = await StateTree.load({dbId: metaCall!.id, config: pconf}).toPromise();
    await loadedTree.init().toPromise();
    const lc = loadedTree.toSerializedState({disableNodesUUID: true, disableCallsUUID: true});
    expectDeepEqual(lc, sc);
  });

});
