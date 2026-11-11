import * as Vue from 'vue';
import {
  isFuncCallState,
  PipelineState,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';

type PanelRef = Readonly<Vue.Ref<{$el: HTMLElement} | null>>;

export function useDockPanels(
  chosenStepUuid: Readonly<Vue.Ref<string | undefined>>,
  chosenStepState: Readonly<Vue.Ref<PipelineState | undefined>>,
  chooseStep: (uuid: string) => void,
  panels: {
    stepsTree: PanelRef,
    inspector: PanelRef,
    rfv: PanelRef,
    pipelineView: Readonly<Vue.Ref<{$el: HTMLElement, showStepsToAdd: () => void} | null>>,
  },
) {
  const treeHidden = Vue.ref(false);
  const inspectorHidden = Vue.ref(true);
  const rfvHidden = Vue.ref(true);
  const pipelineViewHidden = Vue.ref(true);

  const toggleTree = () => {
    treeHidden.value = !treeHidden.value;
  };
  const toggleInspector = () => {
    inspectorHidden.value = !inspectorHidden.value;
  };

  const handlePanelClose = (el: HTMLElement) => {
    if (el === panels.stepsTree.value?.$el) treeHidden.value = true;
    if (el === panels.inspector.value?.$el) inspectorHidden.value = true;
    if (el === panels.rfv.value?.$el) rfvHidden.value = true;
    if (el === panels.pipelineView.value?.$el) pipelineViewHidden.value = true;
  };

  // a step opened by Back/Next takes the focus once the dock has placed its view, so the keyboard stays in the
  // workflow instead of falling back to the page
  let focusOpenedStep = false;
  const navigateTo = (uuid: string) => {
    focusOpenedStep = uuid !== chosenStepUuid.value;
    chooseStep(uuid);
  };

  // chooses a step and shows its panel again if it was closed
  const showStep = (state: PipelineState) => {
    chooseStep(state.uuid);
    if (isFuncCallState(state))
      rfvHidden.value = false;
    else
      pipelineViewHidden.value = false;
  };

  // the keyboard way to a tree row's plus: the workflow's steps to add
  const chooseStepToAdd = async (state: PipelineState) => {
    showStep(state);
    await Vue.nextTick();
    panels.pipelineView.value?.showStepsToAdd();
  };

  const onPanelOpened = (el: HTMLElement) => {
    if (focusOpenedStep && (el === panels.rfv.value?.$el || el === panels.pipelineView.value?.$el)) {
      focusOpenedStep = false;
      el.focus({preventScroll: true});
    }
  };

  Vue.watch(chosenStepState, (state) => {
    if (!state)
      return;
    if (isFuncCallState(state)) rfvHidden.value = false;
    if (!isFuncCallState(state)) pipelineViewHidden.value = false;
  });

  return {
    treeHidden: Vue.readonly(treeHidden),
    inspectorHidden: Vue.readonly(inspectorHidden),
    rfvHidden: Vue.readonly(rfvHidden),
    pipelineViewHidden: Vue.readonly(pipelineViewHidden),
    toggleTree, toggleInspector, handlePanelClose, navigateTo, showStep, chooseStepToAdd, onPanelOpened,
  };
}
