import * as grok from 'datagrok-api/grok';
import * as Vue from 'vue';
import {
  isFuncCallState,
  PipelineState,
  ViewAction,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import * as Utils from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {busyReason, Control, READONLY_REASON, runControl} from '../../utils';
import {getRelevantGlobalActions, hasInconsistencies, hasSubtreeAnyInconsistencies} from './tree-utils';
import type {DriverAPI} from '../../composables/use-reactive-tree-driver';

export function useTreeControls(
  {treeMutationsLocked, isGlobalLocked, treeState, states, hiddenActions}: DriverAPI,
  chosenStepUuid: Readonly<Vue.Ref<string | undefined>>,
  chosenStepState: Readonly<Vue.Ref<PipelineState | undefined>>,
) {
  // Ribbon controls key their visibility on tree presence only, not on the transient locks.
  // treeMutationsLocked/isGlobalLocked cycle on every link-propagation burst and every run,
  // which would otherwise hide/re-show the controls (flicker). treeState is buffered during
  // the global lock, so its last-good value keeps the controls' content stable while busy.
  const isTreeLoaded = Vue.computed(() => !!treeState.value);
  const treeBusy = Vue.computed(() => treeMutationsLocked.value || isGlobalLocked.value);

  // The ribbon lives in the view chrome, outside the compositor overlay that covers the
  // component body while busy, so the mutating controls stay visually unchanged during a lock
  // and surface a warning on click (the driver otherwise no-ops them with a console warning).
  const busyGuard = (action: string) => ({
    disabled: treeBusy.value,
    disabledReason: busyReason(action),
    disabledStyle: 'none' as const,
    disabledReasonMode: 'popup' as const,
  });

  const menuActions = Vue.computed(() => {
    if (!treeState.value || !chosenStepUuid.value)
      return {};
    const globalActions = getRelevantGlobalActions(treeState.value, chosenStepUuid.value);
    const currentStepActions = chosenStepState.value?.actions?.filter((action) => action.position === 'menu') ?? [];
    const actions = [...globalActions, ...currentStepActions]
      .filter((action) => !hiddenActions.value?.has(action.uuid));
    return actions.reduce((acc, action) => {
      const menuCategory = action.menuCategory ?? 'Actions';
      if (action.position === 'menu' || action.position === 'globalmenu') {
        if (acc[menuCategory])
          acc[menuCategory].push(action);
        else
          acc[menuCategory] = [action];
      }
      return acc;
    }, {} as Record<string, ViewAction[]>);
  });

  const buttonActions = Vue.computed(() => {
    return chosenStepState.value?.actions?.reduce((acc, action) => {
      if (action.position === 'buttons' && !hiddenActions.value?.has(action.uuid))
        acc.push(action);

      return acc;
    }, [] as ViewAction[]);
  });

  // actions on the whole workflow warn while it is busy, like its ribbon controls
  const whenIdle = (action: string, run: () => void) =>
    treeBusy.value ? grok.shell.warning(busyReason(action)) : run();

  // a step's Run and Update buttons; a workflow has no Run, and its Update is the row's substeps update
  const runControlOf = (state: PipelineState): Control => {
    if (!isFuncCallState(state) || !state.funcCall)
      return {shown: false};
    const callState = states.calls[state.uuid];
    return runControl({
      isOutputOutdated: callState?.isOutputOutdated,
      allowRerun: Utils.getFeature(Utils.getFeatures(state.funcCall.func), 'rerun', false),
      autorun: false,
      isRunnable: callState?.isRunnable,
      isRunning: callState?.isRunning,
      isReadonly: state.isReadonly,
      isBlocked: treeBusy.value,
    });
  };
  const updateControlOf = (state: PipelineState): Control => ({
    shown: isFuncCallState(state) ? !!hasInconsistencies(states.consistency[state.uuid]) :
      !!hasSubtreeAnyInconsistencies(state, states.calls, states.consistency),
    disabledReason: state.isReadonly ? READONLY_REASON : treeBusy.value ? busyReason('updating') : undefined,
  });
  const updateControl = Vue.computed(() =>
    chosenStepState.value ? updateControlOf(chosenStepState.value) : {shown: false});
  const runControlOfChosen = Vue.computed(() =>
    chosenStepState.value ? runControlOf(chosenStepState.value) : {shown: false} as Control);
  const runLabel = Vue.computed(() => {
    const state = chosenStepState.value;
    const func = state && isFuncCallState(state) ? state.funcCall?.func : undefined;
    return (func ? Utils.getRunLabel(func) : null) ?? 'Run';
  });

  return {
    isTreeLoaded, treeBusy, busyGuard, menuActions, buttonActions, whenIdle,
    runControlOf, updateControlOf, updateControl, runControlOfChosen, runLabel,
  };
}
