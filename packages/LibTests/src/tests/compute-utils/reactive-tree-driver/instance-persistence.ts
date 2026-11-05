import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test} from '@datagrok-libraries/test/src/test';
import {PipelineConfiguration, historyUtils} from '@datagrok-libraries/compute-utils';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {Driver} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/Driver';
import {FuncCallNode, PipelineNodeBase} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {deserializeRestrictions} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/funccall-utils';
import {callHandler} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/utils';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {serialize} from '@datagrok-libraries/utils/src/json-serialization';
import {PipelineStateStatic, StepFunCallState} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {deleteRuns, savedCallIds} from '../../../test-utils';

category('ComputeUtils: Driver state tree persistence', async () => {
  test('Load and save simple config', async () => {
    const config = await callHandler<PipelineConfiguration>('LibTests:MockWrapper1', {}).toPromise();
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const sc = tree.toSerializedState({disableNodesUUID: true, disableCallsUUID: true});
    const metaCall = await tree.save().toPromise();
    try {
      const loadedTree = await StateTree.load({dbId: metaCall!.id, config: pconf}).toPromise();
      await loadedTree.init().toPromise();
      const lc = loadedTree.toSerializedState({disableNodesUUID: true, disableCallsUUID: true});
      expectDeepEqual(lc, sc);
    } finally {
      await deleteRuns([metaCall?.id, ...savedCallIds(tree.toSerializedState())]);
    }
  });

  test('Save creates new uuids', async () => {
    const config = await callHandler<PipelineConfiguration>('LibTests:MockWrapper1', {}).toPromise();
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const sc = tree.toSerializedState();
    const metaCall = await tree.save().toPromise();
    try {
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
    } finally {
      await deleteRuns([metaCall?.id, ...savedCallIds(tree.toSerializedState())]);
    }
  });

  test('Load nested pipeline', async () => {
    // create and save nested pipeline
    const config = await callHandler<PipelineConfiguration>('LibTests:MockWrapper2', {version: '1.0'}).toPromise();
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const sc = tree.toSerializedState({disableNodesUUID: true, disableCallsUUID: true});
    const metaCall = await tree.save().toPromise();
    try {
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
      expectDeepEqual(lc, {...sc, sourceRunId: metaCall!.id});
    } finally {
      await deleteRuns([metaCall?.id, ...savedCallIds(tree.toSerializedState())]);
    }
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
    try {
      const config = await callHandler<PipelineConfiguration>('LibTests:MockWrapper2', {version: '1.0'}).toPromise();
      const pconf = await getProcessedConfig(config);
      const loadedTree = await StateTree.load({dbId: metaCall!.id, config: pconf}).toPromise();
      await loadedTree.init().toPromise();
      const lc = loadedTree.toSerializedState({disableNodesUUID: true, disableCallsUUID: true});
      expectDeepEqual(lc, sc);
    } finally {
      await deleteRuns([metaCall?.id, ...savedCallIds(StateTree.toSerializedStateRec(nestedRoot))]);
    }
  });

  test('Save and load pipeline with action step', async () => {
    const config = await callHandler<PipelineConfiguration>('LibTests:MockWrapperAction', {version: '1.0'}).toPromise();
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const sc = tree.toSerializedState({disableNodesUUID: true, disableCallsUUID: true});
    const metaCall = await tree.save().toPromise();
    try {
      const loadedTree = await StateTree.load({dbId: metaCall!.id, config: pconf}).toPromise();
      await loadedTree.init().toPromise();
      const lc = loadedTree.toSerializedState({disableNodesUUID: true, disableCallsUUID: true});
      expectDeepEqual(lc, sc);
      // the action node survives the round-trip and is still marked as an action step in view state
      const actionNode = loadedTree.nodeTree.getNode([{idx: 1}]);
      const actionState = actionNode.getItem().toState({}) as any;
      expectDeepEqual(actionState.type, 'static');
      expectDeepEqual(actionState.isActionStep, true);
    } finally {
      await deleteRuns([metaCall?.id, ...savedCallIds(tree.toSerializedState())]);
    }
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
    try {
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
    } finally {
      await deleteRuns([metaCall?.id, ...savedCallIds(tree.toSerializedState())]);
    }
  });

  test('Driver loadPipeline with explicit config keeps run metadata', async () => {
    const config = await callHandler<PipelineConfiguration>('LibTests:MockWrapper1', {}).toPromise();
    const pconf = await getProcessedConfig(config);
    const driver = new Driver();
    try {
      await driver.sendCommand({event: 'initPipeline', provider: '', config: pconf});
      await driver.sendCommand({event: 'savePipeline', title: 'meta title', description: 'meta desc'});
      const savedId = driver.currentMetaCallData$.value?.id;
      expectDeepEqual(!!savedId, true, {prefix: 'Saved metaCall id'});
      await driver.sendCommand({event: 'loadPipeline', funcCallId: savedId!, config: pconf});
      const meta = driver.currentMetaCallData$.value;
      expectDeepEqual(meta?.id, savedId, {prefix: 'Loaded metaCall id'});
      expectDeepEqual(meta?.title, 'meta title', {prefix: 'Loaded title'});
      expectDeepEqual(meta?.description, 'meta desc', {prefix: 'Loaded description'});
    } finally {
      const steps = (driver.currentState$.value as PipelineStateStatic<StepFunCallState, {}> | undefined)?.steps ?? [];
      const stepIds = steps.map((step) => (step as StepFunCallState).funcCall?.id);
      await deleteRuns([driver.currentMetaCallData$.value?.id, ...stepIds]);
      driver.close();
    }
  });

  test('Driver savePipeline with a space target links only the workflow run there', async () => {
    const config = await callHandler<PipelineConfiguration>('LibTests:MockWrapper1', {}).toPromise();
    const pconf = await getProcessedConfig(config);
    const driver = new Driver();
    const space = await grok.dapi.spaces.createRootSpace(`Driver save target ${Date.now()}`);
    try {
      await driver.sendCommand({event: 'initPipeline', provider: '', config: pconf});
      await driver.sendCommand({event: 'savePipeline', title: 'space title', target: {spaceId: space.id}});
      const meta = driver.currentMetaCallData$.value;
      expectDeepEqual(meta?.title, 'space title', {prefix: 'Saved title'});
      const linked = await historyUtils.listSpaceRuns(space.id);
      expectDeepEqual(linked.map((run) => run.id), [meta?.id], {prefix: 'Linked runs'});
    } finally {
      const steps = (driver.currentState$.value as PipelineStateStatic<StepFunCallState, {}> | undefined)?.steps ?? [];
      const stepIds = steps.map((step) => (step as StepFunCallState).funcCall?.id);
      await grok.dapi.spaces.delete(space);
      await deleteRuns([driver.currentMetaCallData$.value?.id, ...stepIds]);
      driver.close();
    }
  });

  test('Driver savePipeline that fails only to link into the space keeps the run', async () => {
    const config = await callHandler<PipelineConfiguration>('LibTests:MockWrapper1', {}).toPromise();
    const pconf = await getProcessedConfig(config);
    const driver = new Driver();
    const space = await grok.dapi.spaces.createRootSpace(`Driver deleted space ${Date.now()}`);
    await grok.dapi.spaces.delete(space);
    try {
      await driver.sendCommand({event: 'initPipeline', provider: '', config: pconf});
      const saved = await driver.sendCommand({event: 'savePipeline', title: 'unlinked', target: {spaceId: space.id}});
      expectDeepEqual(saved != null, true, {prefix: 'Save resolves with the run'});
      const meta = driver.currentMetaCallData$.value;
      expectDeepEqual(meta?.title, 'unlinked', {prefix: 'Saved title'});
      expectDeepEqual((await historyUtils.loadRun(meta!.id!)).id, meta?.id, {prefix: 'Saved run'});
    } finally {
      const steps = (driver.currentState$.value as PipelineStateStatic<StepFunCallState, {}> | undefined)?.steps ?? [];
      const stepIds = steps.map((step) => (step as StepFunCallState).funcCall?.id);
      await deleteRuns([driver.currentMetaCallData$.value?.id, ...stepIds]);
      driver.close();
    }
  });

  test('Read-only nested load records its source and stays read-only after save', async () => {
    const config = await callHandler<PipelineConfiguration>('LibTests:MockWrapper2', {version: '1.0'}).toPromise();
    const tree = StateTree.fromPipelineConfig({config: await getProcessedConfig(config)});
    await tree.init().toPromise();
    const metaCall = await tree.save().toPromise();
    const outerConfig = await callHandler<PipelineConfiguration>('LibTests:MockWrapper3', {version: '1.0'}).toPromise();
    const outerPconf = await getProcessedConfig(outerConfig);
    const outerTree = StateTree.fromPipelineConfig({config: outerPconf});
    let outerMetaCall: DG.FuncCall | undefined;
    try {
      await outerTree.init().toPromise();
      const root = outerTree.nodeTree.getItem([]);
      await outerTree.loadSubTree(root.uuid, metaCall!.id, 'pipelinePar', 1, true).toPromise();
      const step = outerTree.nodeTree.getItem([{idx: 1}, {idx: 0}]);
      let error: string | undefined;
      try {
        await outerTree.runStep(step.uuid).toPromise();
      } catch (e: any) {
        error = e.message;
      }
      expectDeepEqual(error?.includes('read-only'), true, {prefix: 'Run'});
      outerMetaCall = await outerTree.save().toPromise();
      const loadedTree = await StateTree.load({dbId: outerMetaCall!.id, config: outerPconf}).toPromise();
      await loadedTree.init().toPromise();
      const loaded = loadedTree.nodeTree.getItem([{idx: 1}]) as PipelineNodeBase;
      const loadedStep = loadedTree.nodeTree.getItem([{idx: 1}, {idx: 0}]) as FuncCallNode;
      expectDeepEqual(loaded.sourceRunId, metaCall!.id, {prefix: 'Source'});
      expectDeepEqual(loaded.isReadonly, true, {prefix: 'Locked item'});
      expectDeepEqual(loadedStep.instancesWrapper.getInstance()?.isReadonly, true, {prefix: 'Locked call'});
      expectDeepEqual(loadedTree.nodeTree.getItem([{idx: 0}]).isReadonly, false, {prefix: 'Editable item'});
    } finally {
      await deleteRuns([metaCall?.id, outerMetaCall?.id, ...savedCallIds(tree.toSerializedState()),
        ...savedCallIds(outerTree.toSerializedState())]);
    }
  });

  test('Duplicate keeps a nested read-only load read-only', async () => {
    const config = await callHandler<PipelineConfiguration>('LibTests:MockWrapper1', {version: '1.0'}).toPromise();
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const metaCall = await tree.save().toPromise();
    try {
      const outerConfig =
        await callHandler<PipelineConfiguration>('LibTests:MockWrapper3', {version: '1.0'}).toPromise();
      const outerTree = StateTree.fromPipelineConfig({config: await getProcessedConfig(outerConfig)});
      await outerTree.init().toPromise();
      const item = outerTree.nodeTree.getItem([{idx: 0}]);
      await outerTree.loadSubTree(item.uuid, metaCall!.id, 'pipeline1', 1, true, true).toPromise();
      await outerTree.duplicateSubtree(item.uuid).toPromise();
      const copy = outerTree.nodeTree.getItem([{idx: 1}, {idx: 1}]);
      const copyStep = outerTree.nodeTree.getItem([{idx: 1}, {idx: 1}, {idx: 0}]) as FuncCallNode;
      expectDeepEqual(copy.isReadonly, true, {prefix: 'Locked item'});
      expectDeepEqual((copy as PipelineNodeBase).sourceRunId, metaCall!.id, {prefix: 'Source'});
      expectDeepEqual(copyStep.isReadonly, true, {prefix: 'Locked step'});
      expectDeepEqual(copyStep.instancesWrapper.getInstance()?.isReadonly, true, {prefix: 'Locked call'});
      expectDeepEqual(outerTree.nodeTree.getItem([{idx: 1}, {idx: 0}]).isReadonly, false, {prefix: 'Editable step'});
    } finally {
      await deleteRuns([metaCall?.id, ...savedCallIds(tree.toSerializedState())]);
    }
  });

  test('Loading into a read-only workflow fails', async () => {
    const config = await callHandler<PipelineConfiguration>('LibTests:MockWrapper1', {version: '1.0'}).toPromise();
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const metaCall = await tree.save().toPromise();
    try {
      const outerConfig =
        await callHandler<PipelineConfiguration>('LibTests:MockWrapper2', {version: '1.0'}).toPromise();
      const outerPconf = await getProcessedConfig(outerConfig);
      const outerTree = StateTree.fromPipelineConfig({config: outerPconf, isReadonly: true});
      await outerTree.init().toPromise();
      const root = outerTree.nodeTree.getItem([]);
      let error: string | undefined;
      try {
        await outerTree.loadSubTree(root.uuid, metaCall!.id, 'pipeline1', 1, false, true).toPromise();
      } catch (e: any) {
        error = e.message;
      }
      expectDeepEqual(error?.includes('read-only'), true, {prefix: 'Error'});
    } finally {
      await deleteRuns([metaCall?.id, ...savedCallIds(tree.toSerializedState())]);
    }
  });

  test('Read-only step load records its source and feeds later steps', async () => {
    const config = await callHandler<PipelineConfiguration>('LibTests:MockWrapper1', {version: '1.0'}).toPromise();
    const pconf = await getProcessedConfig(config);
    const source = DG.Func.byName('LibTests:TestAdd2').prepare({a: 2, b: 3});
    await source.call();
    source.options['title'] = 'Read-only source';
    await historyUtils.saveRun(source, {newId: true});
    const tree = StateTree.fromPipelineConfig({config: pconf});
    let metaCall: DG.FuncCall | undefined;
    try {
      await tree.init().toPromise();
      await tree.loadSubTree(tree.nodeTree.getItem([]).uuid, source.id, 'step1', 0, true, true).toPromise();
      const loaded = tree.nodeTree.getItem([{idx: 0}]) as FuncCallNode;
      expectDeepEqual(loaded.isReadonly, true, {prefix: 'Locked step'});
      expectDeepEqual(loaded.sourceRunId, source.id, {prefix: 'Source'});
      expectDeepEqual(loaded.sourceRunTitle, 'Read-only source', {prefix: 'Source title'});
      expectDeepEqual(loaded.getStateStore().getState('res'), 5, {prefix: 'Loaded output'});
      expectDeepEqual(tree.nodeTree.getItem([{idx: 1}]).getStateStore().getState('a'), 5, {prefix: 'Later step'});
      let error: string | undefined;
      try {
        await tree.runStep(loaded.uuid).toPromise();
      } catch (e: any) {
        error = e.message;
      }
      expectDeepEqual(error?.includes('read-only'), true, {prefix: 'Run'});
      metaCall = await tree.save().toPromise();
      const loadedTree = await StateTree.load({dbId: metaCall!.id, config: pconf}).toPromise();
      await loadedTree.init().toPromise();
      const reloaded = loadedTree.nodeTree.getItem([{idx: 0}]) as FuncCallNode;
      expectDeepEqual(reloaded.isReadonly, true, {prefix: 'Reloaded locked step'});
      expectDeepEqual(reloaded.sourceRunId, source.id, {prefix: 'Reloaded source'});
      expectDeepEqual(reloaded.sourceRunTitle, 'Read-only source', {prefix: 'Reloaded source title'});
      expectDeepEqual(JSON.stringify(tree.toSerializedState()).includes('sourceRunTitle'), false,
        {prefix: 'Not saved'});
      expectDeepEqual(reloaded.instancesWrapper.getInstance()?.isReadonly, true, {prefix: 'Reloaded locked call'});
    } finally {
      await deleteRuns([source.id, metaCall?.id, ...savedCallIds(tree.toSerializedState())]);
    }
  });

  test('Editable step load records its source and stays runnable', async () => {
    const config = await callHandler<PipelineConfiguration>('LibTests:MockWrapper1', {version: '1.0'}).toPromise();
    const pconf = await getProcessedConfig(config);
    const source = DG.Func.byName('LibTests:TestMul2').prepare({a: 7, b: 2});
    await source.call();
    source.options['title'] = 'Editable source';
    await historyUtils.saveRun(source, {newId: true});
    const tree = StateTree.fromPipelineConfig({config: pconf});
    try {
      await tree.init().toPromise();
      await tree.loadSubTree(tree.nodeTree.getItem([]).uuid, source.id, 'step2', 1, false, true).toPromise();
      const loaded = tree.nodeTree.getItem([{idx: 1}]) as FuncCallNode;
      expectDeepEqual(loaded.isReadonly, false, {prefix: 'Editable step'});
      expectDeepEqual(loaded.sourceRunId, source.id, {prefix: 'Source'});
      expectDeepEqual(loaded.sourceRunTitle, 'Editable source', {prefix: 'Source title'});
      expectDeepEqual(loaded.getStateStore().getState('res'), 14, {prefix: 'Loaded output'});
      loaded.getStateStore().editState('b', 3);
      await tree.runStep(loaded.uuid).toPromise();
      expectDeepEqual(loaded.getStateStore().getState('res'), 21, {prefix: 'Rerun output'});
    } finally {
      await deleteRuns([source.id, ...savedCallIds(tree.toSerializedState())]);
    }
  });

  test('A deleted source run leaves the loaded step without a title', async () => {
    const config = await callHandler<PipelineConfiguration>('LibTests:MockWrapper1', {version: '1.0'}).toPromise();
    const pconf = await getProcessedConfig(config);
    const source = DG.Func.byName('LibTests:TestAdd2').prepare({a: 2, b: 3});
    await source.call();
    source.options['title'] = 'Deleted source';
    await historyUtils.saveRun(source, {newId: true});
    const tree = StateTree.fromPipelineConfig({config: pconf});
    let metaCall: DG.FuncCall | undefined;
    let sourceDeleted = false;
    try {
      await tree.init().toPromise();
      await tree.loadSubTree(tree.nodeTree.getItem([]).uuid, source.id, 'step1', 0, true, true).toPromise();
      metaCall = await tree.save().toPromise();
      await deleteRuns([source.id]);
      sourceDeleted = true;
      const loadedTree = await StateTree.load({dbId: metaCall!.id, config: pconf}).toPromise();
      await loadedTree.init().toPromise();
      const reloaded = loadedTree.nodeTree.getItem([{idx: 0}]) as FuncCallNode;
      expectDeepEqual(reloaded.sourceRunId, source.id, {prefix: 'Source'});
      expectDeepEqual(reloaded.sourceRunTitle, null, {prefix: 'No title'});
      expectDeepEqual(reloaded.getStateStore().getState('res'), 5, {prefix: 'Own copy loaded'});
    } finally {
      await deleteRuns(
        [sourceDeleted ? undefined : source.id, metaCall?.id, ...savedCallIds(tree.toSerializedState())]);
    }
  });

  test('Duplicated run step saves and loads as run', async () => {
    const pconf = await getProcessedConfig({
      id: 'pipeline1',
      type: 'sequential',
      nqName: 'LibTests:MockWrapper5',
      version: '1.0',
      stepTypes: [{id: 'mul', nqName: 'LibTests:TestMul2'}],
      initialSteps: ['mul'],
    });
    const tree = StateTree.fromPipelineConfig({config: pconf});
    await tree.init().toPromise();
    const source = tree.nodeTree.getItem([{idx: 0}]) as FuncCallNode;
    source.getStateStore().editState('a', 2);
    source.getStateStore().editState('b', 3);
    await tree.runStep(source.uuid).toPromise();
    await tree.duplicateSubtree(source.uuid).toPromise();
    const metaCall = await tree.save().toPromise();
    try {
      const loadedTree = await StateTree.load({dbId: metaCall!.id, config: pconf}).toPromise();
      await loadedTree.init().toPromise();
      const loadedSource = loadedTree.nodeTree.getItem([{idx: 0}]) as FuncCallNode;
      const copy = loadedTree.nodeTree.getItem([{idx: 1}]) as FuncCallNode;
      expectDeepEqual(copy.funcCallState$.value?.isOutputOutdated, false, {prefix: 'Run copy'});
      expectDeepEqual(copy.getStateStore().getState('res'), 6, {prefix: 'Output'});
      expectDeepEqual(copy.getStateStore().getState('a'), 2, {prefix: 'Input'});
      expectDeepEqual(!!copy.instancesWrapper.getInstance()?.getFuncCall().started, true, {prefix: 'Started'});
      expectDeepEqual(copy.instancesWrapper.id !== loadedSource.instancesWrapper.id, true, {prefix: 'Own call'});
    } finally {
      await deleteRuns([metaCall?.id, ...savedCallIds(tree.toSerializedState())]);
    }
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
