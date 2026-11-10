import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {useUrlSearchParams} from '@vueuse/core';
import {PipelineState} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {pinView} from '../../utils';
import {
  findNodeWithPathByUuid, findTreeNodeByPath, nodeMetadata, resolveChosenUuid, resolveSingleStep,
} from './tree-utils';
import type {DriverAPI} from '../../composables/use-reactive-tree-driver';

// The chosen step repair must run before the currentStep sync: it reads searchParams.currentStep as its
// fallback path, which the sync clears for a removed step. It does, as it subscribes to treeState first.
export function useTreeRouting(
  {treeState, states, currentMetaCallData, hasNotSavedEdits, loadPipeline}: DriverAPI,
  {currentView, singleStep, modelName, initialRunId}: {
    currentView: DG.ViewBase,
    singleStep: Readonly<Vue.Ref<ReturnType<typeof resolveSingleStep>>>,
    modelName: Readonly<Vue.Ref<string>>,
    initialRunId?: string,
  },
) {
  const chosenStepUuid = Vue.ref<string | undefined>();
  const searchParams = useUrlSearchParams<{id?: string, currentStep?: string}>('history');

  const setViewName = (name: string = '') => {
    currentView.name = name;
  };

  const setViewPath = (path: string = '') => {
    currentView.path = path;
  };

  Vue.watch(searchParams, (params) => {
    const paramsRaw = [];
    if (params.currentStep)
      paramsRaw.push(`currentStep=${params.currentStep.replace(/ /g, '+')}`);
    if (params.id)
      paramsRaw.push(`id=${params.id}`);
    setViewPath(paramsRaw.length ? `?${paramsRaw.join('&')}`: '?');
  });

  const rootTitle = Vue.computed(() => nodeMetadata(states.descriptions, treeState.value?.uuid).title);

  Vue.watch([currentMetaCallData, hasNotSavedEdits, rootTitle], ([metadata, hasNotSavedEdits, rootTitle]) => {
    if (hasNotSavedEdits)
      pinView(currentView);
    const {id, title} = metadata?.id ? metadata : {id: undefined, title: undefined};
    searchParams.id = id && !hasNotSavedEdits ? id : undefined;
    if (title)
      setViewName(`${modelName.value} - ${title}${hasNotSavedEdits ? ' (edited)' : ''}`);
    else if (!id && rootTitle)
      setViewName(`${modelName.value} - ${rootTitle}`);
    else
      setViewName(modelName.value);
  });


  const setCurrentStepByPath = (pendingStepPath: string, treeState: PipelineState) => {
    if (pendingStepPath) {
      const nodeWithPath = findTreeNodeByPath(
        pendingStepPath.split(' ').map((pathSegment) => Number.parseInt(pathSegment)),
        treeState,
      );
      if (nodeWithPath?.state)
        chosenStepUuid.value = nodeWithPath.state.uuid;
    }
  };

  let pendingStepPath: string | null = null;
  // Programmatic open (OpenWorkflowRun) passes the id via call.aux, deep links via the start URL
  let runIdToLoad = initialRunId;

  Vue.watch(treeState, (treeState) => {
    if (!treeState)
      return;

    if (pendingStepPath) {
      if (!singleStep.value)
        setCurrentStepByPath(pendingStepPath, treeState);
      pendingStepPath = null;
      return;
    }

    if (!runIdToLoad && globalThis.initialURLHandled)
      return;

    // Getting inital URL user entered with
    const startUrl = new URL(grok.shell.startUri);
    const loadingId = runIdToLoad ?? startUrl.searchParams.get('id');
    pendingStepPath = (runIdToLoad || singleStep.value) ? null : startUrl.searchParams.get('currentStep');
    runIdToLoad = undefined;
    globalThis.initialURLHandled = true;

    if (loadingId)
      loadPipeline(loadingId);
    else if (pendingStepPath) {
      setCurrentStepByPath(pendingStepPath, treeState);
      pendingStepPath = null;
    }
  }, {immediate: true});

  const chosenStep = Vue.computed(() => {
    if (!treeState.value || !chosenStepUuid.value)
      return null;
    return findNodeWithPathByUuid(chosenStepUuid.value, treeState.value) ?? null;
  });

  // Owns chosenStepUuid: defaults it on first load and repairs it when the
  // selected node disappears (e.g. the selected step was removed). Falls back
  // along the old positional path to the nearest surviving ancestor, so a
  // removed last child resolves to its parent instead of dangling.
  Vue.watch([treeState, chosenStepUuid], ([treeState]) => {
    if (!treeState)
      return;
    if (singleStep.value) {
      if (chosenStepUuid.value !== singleStep.value.step.uuid)
        chosenStepUuid.value = singleStep.value.step.uuid;
      return;
    }
    const fallbackPath = searchParams.currentStep ?
      searchParams.currentStep.split(' ').map((segment) => Number.parseInt(segment)) :
      undefined;
    const next = resolveChosenUuid(chosenStepUuid.value, treeState, fallbackPath);
    if (next !== chosenStepUuid.value)
      chosenStepUuid.value = next;
  }, {immediate: true});

  Vue.watch(chosenStep, (newStep) => {
    if (newStep && !singleStep.value)
      searchParams.currentStep = newStep.pathSegments.join(' ');
    else
      searchParams.currentStep = undefined;
  }, {immediate: true});

  const chosenStepState = Vue.computed(() => chosenStep.value?.state);

  const chooseStep = (uuid: string) => {
    chosenStepUuid.value = uuid;
  };

  return {chosenStepUuid: Vue.readonly(chosenStepUuid), chosenStepState, chooseStep};
}
