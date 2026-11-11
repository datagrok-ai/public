import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import type * as Vue from 'vue';
import {
  isFuncCallState,
  PipelineState,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import type {MetaCallInfo} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {snapshotState} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/state-snapshot';
import {historyUtils} from '@datagrok-libraries/compute-utils/history-utils/src/history-utils';
import {STEP_HISTORY_OPTION, asItemMetadata, defaultRunName} from '../RunStorage/run-storage';
import {openSaveRunDialog} from '../RunStorage/SaveRunDialog';
import {entryFromSnapshot} from '../RunComparison/entry-extraction';
import {
  addRunToComparison, addToComparison, compareWithOriginal, snapshotName,
} from '../RunComparison/comparison-views';
import {findNodeWithPathByUuid, nodeMetadata, nodeName} from './tree-utils';
import type {DriverAPI} from '../../composables/use-reactive-tree-driver';

export function createNodeRuns(
  {treeState, states, currentMetaCallData, hasNotSavedEdits, savePipeline, saveDynamicItem}: DriverAPI,
  chosenStepState: Readonly<Vue.Ref<PipelineState | undefined>>,
) {
  const snapshotOf = (state: PipelineState) => entryFromSnapshot(
    snapshotState(state, states.calls, states.consistency),
    snapshotName(nodeName(states.descriptions, state)), {
      metadata: nodeMetadata(states.descriptions, state.uuid,
        state.uuid === treeState.value?.uuid ? currentMetaCallData.value : undefined),
      stepHistory: isFuncCallState(state) && !!state.enableHistory,
    });
  // the function whose history holds the node's saved runs
  const modelOf = (state: PipelineState) => isFuncCallState(state) ?
    state.funcCall?.func : (state.nqName ? DG.Func.byName(state.nqName) : undefined);
  const sourceTitle = (state: PipelineState) => isFuncCallState(state) ? state.sourceRunTitle : undefined;
  const compareNodeWithOriginal = (state: PipelineState) =>
    compareWithOriginal(snapshotOf(state), state.sourceRunId!, sourceTitle(state), modelOf(state));
  // an unedited saved model and read-only copies equal their saved runs; editable copies aren't tracked
  const savedRunOf = (state: PipelineState) => state.uuid === treeState.value?.uuid ?
    (hasNotSavedEdits.value ? undefined : currentMetaCallData.value?.id) :
    (state.isReadonly ? state.sourceRunId : undefined);
  const addNodeToComparison = (state: PipelineState) => {
    const runId = savedRunOf(state);
    if (runId) {
      addRunToComparison(runId, state.uuid === treeState.value?.uuid ? currentMetaCallData.value?.title :
        sourceTitle(state), modelOf(state));
    } else
      addToComparison(snapshotOf(state), modelOf(state));
  };
  const addModelToComparison = () => addNodeToComparison(treeState.value!);

  // the Save dialog starts from the loaded run's name, description and tags, then the node's links,
  // then the model name and date
  const askSaveRun = (state: PipelineState, loaded?: MetaCallInfo) => {
    const func = modelOf(state);
    if (!func)
      throw new Error(`Node ${state.uuid} has no model to save`);
    const {title, description, tags} = nodeMetadata(states.descriptions, state.uuid);
    return openSaveRunDialog({
      func,
      name: defaultRunName(func, loaded?.title, title),
      description: loaded?.description || description,
      tags: loaded?.tags?.length ? loaded.tags : tags,
    });
  };

  const saveSubTreeState = async (uuid: string) => {
    const state = findNodeWithPathByUuid(uuid, treeState.value!)?.state;
    const result = state ? await askSaveRun(state) : null;
    if (!result)
      return;
    await saveDynamicItem(uuid, asItemMetadata(result));
  };

  const saveEntireModelState = async () => {
    if (!treeState.value) return;
    const loaded = currentMetaCallData.value;
    const result = await askSaveRun(treeState.value, loaded?.id ? loaded : undefined);
    if (!result)
      return;
    await savePipeline(asItemMetadata(result));
  };

  const saveStepToHistory = async (fc: DG.FuncCall) => {
    const state = chosenStepState.value;
    const result = state ? await askSaveRun(state) : null;
    if (!result)
      return;
    fc.options[STEP_HISTORY_OPTION] = 'true';
    try {
      await historyUtils.saveRun(fc, {newId: true, metadata: result.metadata, target: result.target});
      grok.shell.info('Step saved');
    } catch (e: any) {
      grok.shell.error(e);
    }
  };

  return {
    compareNodeWithOriginal, addNodeToComparison, addModelToComparison,
    saveSubTreeState, saveEntireModelState, saveStepToHistory,
  };
}
