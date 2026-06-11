import * as DG from 'datagrok-api/dg';
import {category, test} from '@datagrok-libraries/test/src/test';
import {PipelineConfiguration, historyUtils} from '@datagrok-libraries/compute-utils';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {deserializeRestrictions} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/funccall-utils';
import {callHandler} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/utils';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {serialize} from '@datagrok-libraries/utils/src/json-serialization';

category('ComputeUtils: Driver state tree persistence', async () => {
  test('Load and save simple config', async () => {
    const config = await callHandler<PipelineConfiguration>('LibTests:MockWrapper1', {}).toPromise();
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
    const config = await callHandler<PipelineConfiguration>('LibTests:MockWrapper1', {}).toPromise();
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
    const config = await callHandler<PipelineConfiguration>('LibTests:MockWrapper2', {version: '1.0'}).toPromise();
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const sc = tree.toSerializedState({disableNodesUUID: true, disableCallsUUID: true});
    const metaCall = await tree.save().toPromise();
    // create outer pipeline
    const outerConfig = await callHandler<PipelineConfiguration>('LibTests:MockWrapper3', {version: '1.0'}).toPromise();
    const outerPconf = await getProcessedConfig(outerConfig);
    const outerTree = StateTree.fromPipelineConfig({config: outerPconf});
    await outerTree.init().toPromise();
    // load nested tree into outer
    const root = outerTree.nodeTree.getItem([]);
    await outerTree.loadSubTree(root.uuid, metaCall!.id, 'pipelinePar', 1, false).toPromise();
    const loadedTree = outerTree.nodeTree.getNode([{idx: 1}]);
    const lc = StateTree.toSerializedStateRec(loadedTree, {disableNodesUUID: true, disableCallsUUID: true});
    expectDeepEqual(lc, sc);
  });

  test('Save nested pipeline', async () => {
    // create outer pipeline
    const outerConfig = await callHandler<PipelineConfiguration>('LibTests:MockWrapper3', {version: '1.0'}).toPromise();
    const outerPconf = await getProcessedConfig(outerConfig);
    const outerTree = StateTree.fromPipelineConfig({config: outerPconf});
    await outerTree.init().toPromise();
    const nestedRoot = outerTree.nodeTree.getNode([{idx: 0}]);
    const sc = StateTree.toSerializedStateRec(nestedRoot, {disableNodesUUID: true, disableCallsUUID: true});
    // save nested pipeline
    const metaCall = await outerTree.save(nestedRoot.getItem().uuid).toPromise();
    const config = await callHandler<PipelineConfiguration>('LibTests:MockWrapper2', {version: '1.0'}).toPromise();
    const pconf = await getProcessedConfig(config);
    const loadedTree = await StateTree.load({dbId: metaCall!.id, config: pconf}).toPromise();
    await loadedTree.init().toPromise();
    const lc = loadedTree.toSerializedState({disableNodesUUID: true, disableCallsUUID: true});
    expectDeepEqual(lc, sc);
  });

  test('Save and load pipeline with action step', async () => {
    const config = await callHandler<PipelineConfiguration>('LibTests:MockWrapperAction', {version: '1.0'}).toPromise();
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const sc = tree.toSerializedState({disableNodesUUID: true, disableCallsUUID: true});
    const metaCall = await tree.save().toPromise();
    const loadedTree = await StateTree.load({dbId: metaCall!.id, config: pconf}).toPromise();
    await loadedTree.init().toPromise();
    const lc = loadedTree.toSerializedState({disableNodesUUID: true, disableCallsUUID: true});
    expectDeepEqual(lc, sc);
    // the action node survives the round-trip and is still marked as an action step in view state
    const actionNode = loadedTree.nodeTree.getNode([{idx: 1}]);
    const actionState = actionNode.getItem().toState({}) as any;
    expectDeepEqual(actionState.type, 'static');
    expectDeepEqual(actionState.isActionStep, true);
  });

  test('Save and load consistency dataframe', async () => {
    const conf = await callHandler<PipelineConfiguration>('LibTests:MockWrapperDF', {version: '1.0'}).toPromise();
    const pconf = await getProcessedConfig(conf);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const df = DG.DataFrame.fromColumns([DG.Column.fromInt32Array('x', new Int32Array([1, 2, 3]))]);
    const step2 = tree.nodeTree.getNode([{idx: 1}]);
    step2.getItem().getStateStore().setState('df', df, 'restricted');
    const metaCall = await tree.save().toPromise();
    // loading must not crash when a step carries a dataframe restriction
    const loadedTree = await StateTree.load({dbId: metaCall!.id, config: pconf}).toPromise();
    await loadedTree.init().toPromise();
    // verify the new on-disk format: an id reference, not an embedded dataframe blob
    const fcId = (loadedTree.toSerializedState() as any).steps[1].funcCallId;
    const fc = await historyUtils.loadRun(fcId);
    const blob: string = fc.options['INPUT_RESTRICTIONS'];
    expectDeepEqual(blob.includes('_DG_CONSISTENCY_DF_REF_'), true, {prefix: 'new format stores a df reference'});
    expectDeepEqual(blob.includes('"DataFrame"'), false, {prefix: 'new format does not embed the dataframe'});
    // the reference resolves back to the original dataframe
    const restored = await deserializeRestrictions(blob);
    expectDeepEqual(restored['df']?.type, 'restricted', {prefix: 'restored restriction type'});
    expectDeepEqual(restored['df']?.assignedValue instanceof DG.DataFrame, true, {prefix: 'restored assignedValue is a dataframe'});
    expectDeepEqual(restored['df']?.assignedValue?.rowCount, 3, {prefix: 'restored dataframe content'});
  });

  test('Load legacy embedded-DF consistency format', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromInt32Array('x', new Int32Array([1, 2, 3]))]);
    // old format embedded the dataframe binary directly in the restrictions blob
    const legacy = serialize({df: {type: 'restricted', assignedValue: df}}, {useJsonDF: false});
    expectDeepEqual(legacy.includes('"DataFrame"'), true, {prefix: 'sanity: legacy blob embeds the dataframe'});
    const restored = await deserializeRestrictions(legacy);
    expectDeepEqual(restored['df']?.assignedValue instanceof DG.DataFrame, true, {prefix: 'legacy dataframe restored'});
    expectDeepEqual(restored['df']?.assignedValue?.get('x', 1), 2, {prefix: 'legacy dataframe content'});
  });
});
